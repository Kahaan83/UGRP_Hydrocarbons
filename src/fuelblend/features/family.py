"""T1: strict family-level (Murcko-scaffold) holdout helpers.

`scaffold_family_key` maps a SMILES to a homolog-family key: the generic Murcko
scaffold for anything with a ring framework, else an open-chain bucket keyed by
unsaturation. Acyclic n-alkanes of every chain length therefore share one key,
and toluene/ethylbenzene share the phenyl scaffold — the lesson 0.56 C was
molecule-group friendly.
"""
from __future__ import annotations

import numpy as np


def scaffold_family_key(smiles) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    m = Chem.MolFromSmiles(str(smiles)) if smiles else None
    if m is None:
        return "unreadable"
    try:
        sc = MurckoScaffold.GetScaffoldForMol(m)
        core = Chem.MolToSmiles(sc) if sc.GetNumAtoms() else ""
    except Exception:
        core = ""
    aromatic_bonds = sum(1 for b in m.GetBonds() if b.GetIsAromatic())
    unsat = sum(1 for b in m.GetBonds() if b.GetBondType() != Chem.BondType.SINGLE)
    if not core or core == "C":  # acyclic chain: bucket by unsaturation class
        return f"open_u{max(unsat, aromatic_bonds)}"
    return "sc:" + core


def leave_family_folds(records, n_folds: int = 5, seed: int = 20261030):
    """Deterministic leave-family-out folds: every family sits in exactly one test fold."""
    fams = sorted({r["family"] for r in records})
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(fams))
    fold_of = {fams[perm[i]]: i % n_folds for i in range(len(fams))}
    buckets = {i: [] for i in range(n_folds)}
    for r in records:
        buckets[fold_of[r["family"]]].append(r)
    folds = []
    for i in range(n_folds):
        folds.append(dict(
            test=buckets[i],
            train=[r for j in range(n_folds) if j != i for r in buckets[j]],
            families=[f for f, k in fold_of.items() if k == i],
        ))
    return folds


def plausibility_report(fuels: dict) -> dict:
    """Guard against input collusions: paraffin-heavy (tie-break heavier) fuels should flash higher."""
    items = [v for v in fuels.values() if "flash_C" in v]
    pairs, viol, smapes = 0, 0, []

    def paraffinity(u):
        d = u.get("dist", {})
        return d.get("n_paraffin", 0) + d.get("iso_paraffin", 0)

    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i], items[j]
            if a is b:
                continue
            if paraffinity(a) != paraffinity(b):
                hi, lo = (a, b) if paraffinity(a) > paraffinity(b) else (b, a)
            else:
                hi, lo = (a, b) if a.get("mw_mean", 0) >= b.get("mw_mean", 0) else (b, a)
            if hi["flash_C"] < lo["flash_C"]:
                viol += 1
            denom = max((abs(hi["flash_C"]) + abs(lo["flash_C"])) / 2.0, 1e-9)
            smapes.append(abs(hi["flash_C"] - lo["flash_C"]) / denom * 100.0)
            pairs += 1
    return {"pairs": pairs, "pairwise_violations": viol,
            "smape_pct": float(np.mean(smapes)) if smapes else float("nan")}
