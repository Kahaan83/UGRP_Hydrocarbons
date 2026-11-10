"""r2 reconstruction of CRC 647 compositions from DESCRIPTORS ONLY (aromatics, diaromatics, olefins, D86).

The first reconstruction took the first registry compound of every (class, carbon number) bin, which put cage
hydrocarbons, polycyclic aromatics and perdeuterated alkanes into jet fuels.  Here every bin is represented by a
structurally realistic fuel molecule (``is_realistic``); a bin without one sends its mass to a fallback lump
instead of to an arbitrary molecule.  Nothing in this module reads a target property.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from fuelblend.mixture.fuel_hindcast import RHO_CLASS, carbon_distribution

log = logging.getLogger(__name__)

KINDS = ("n_paraffin", "iso_paraffin", "cycloalkane", "aromatic_mono", "aromatic_di", "alkene_alkyne")


def _saturated(m) -> bool:
    return all(b.GetBondType() == Chem.BondType.SINGLE and not b.GetIsAromatic() for b in m.GetBonds())


def _aromatic_rings(m):
    """Bond-index sets of the rings whose bonds are all aromatic."""
    return [set(r) for r in m.GetRingInfo().BondRings() if all(m.GetBondWithIdx(b).GetIsAromatic() for b in r)]


def is_realistic(smiles: str, kind: str) -> bool:
    m = Chem.MolFromSmiles(smiles)
    if m is None or kind not in KINDS:
        return False
    if any(a.GetSymbol() not in ("C", "H") or a.GetIsotope() != 0 or a.GetFormalCharge() != 0 for a in m.GetAtoms()):
        return False
    ri = m.GetRingInfo()
    n_rings = ri.NumRings()
    carbons = [a for a in m.GetAtoms() if a.GetSymbol() == "C"]
    if kind in ("n_paraffin", "iso_paraffin"):
        if n_rings or not _saturated(m):
            return False
        deg = [a.GetDegree() for a in carbons]
        if kind == "n_paraffin":
            return max(deg) <= 2
        return max(deg) <= 3 and 1 <= sum(d == 3 for d in deg) <= 3
    if kind == "cycloalkane":
        return (_saturated(m) and n_rings in (1, 2) and all(len(r) in (5, 6) for r in ri.AtomRings())
                and rdMolDescriptors.CalcNumBridgeheadAtoms(m) == 0 and rdMolDescriptors.CalcNumSpiroAtoms(m) == 0)
    if kind == "alkene_alkyne":
        n_double = sum(b.GetBondType() == Chem.BondType.DOUBLE for b in m.GetBonds())
        n_triple = sum(b.GetBondType() == Chem.BondType.TRIPLE for b in m.GetBonds())
        return n_rings == 0 and n_double == 1 and n_triple == 0 and not any(b.GetIsAromatic() for b in m.GetBonds())
    arom = _aromatic_rings(m)
    if kind == "aromatic_di":
        return n_rings == 2 and len(arom) == 2 and len(arom[0] & arom[1]) == 1
    # aromatic_mono: one aromatic ring; a second ring must be saturated apart from the bonds it shares with it
    if len(arom) != 1 or n_rings > 2:
        return False
    for r in ri.BondRings():
        if set(r) == arom[0]:
            continue
        for b in r:
            if b in arom[0]:
                continue
            bond = m.GetBondWithIdx(b)
            if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                return False
    return True


class RealisticReps:
    """One realistic registry compound per (kind, carbon number) bin (lowest InChIKey).

    Candidates must be usable in the registry, have valid pure eta(T)/rho(T) in the provider, carry an Antoine
    set and an LFL value (``vp_lookup(smiles)`` is not None, so one composition serves every downstream
    property), and pass ``is_realistic``.
    """

    def __init__(self, reg, provider, vp_lookup, nc=None):
        self.ik, self.smiles = list(reg.ik), list(reg.smiles)
        has_extra = np.array([vp_lookup(s) is not None for s in self.smiles])
        if nc is None:
            nc = np.array([sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(s).GetAtoms()) for s in self.smiles])
        cand = [i for i in range(len(self.ik)) if reg.usable[i] and has_extra[i] and np.isfinite(nc[i])]
        provider.ensure(cand)
        self.by: dict[tuple[str, int], list[int]] = {}
        for i in (i for i in cand if provider.valid[i]):
            for kind in KINDS:
                if is_realistic(self.smiles[i], kind):
                    self.by.setdefault((kind, int(nc[i])), []).append(i)
        for k in self.by:
            self.by[k].sort(key=lambda i: self.ik[i])
        self.ns = {k: sorted(n for (kk, n) in self.by if kk == k) for k in KINDS}

    def get(self, kind: str, n_carbon: int, max_gap: int = 3):
        ns = self.ns.get(kind, [])
        if not ns:
            return None
        best = min(ns, key=lambda k: (abs(k - n_carbon), k))
        return self.by[(kind, best)][0] if abs(best - n_carbon) <= max_gap else None


def nominal_sat_split_detail(jet_class_csv):
    """((n, iso, cyclo), used_fuel_ids, [(skipped_id, reason)])."""
    from fuelblend.data.loaders_jet import included_jet_fuel_ids

    df = pd.read_csv(jet_class_csv)
    used, skipped, shares = [], [], []
    for fid in included_jet_fuel_ids(df):
        cf = json.loads(df.loc[df.fuel_id == fid, "class_fractions"].iloc[0])
        sat = [cf.get(k, np.nan) for k in ("n_paraffin", "iso_paraffin", "cycloalkane")]
        if abs(sum(cf.values()) - 1.0) > 0.02:
            why = "class shares do not sum to 1"
        elif not np.isfinite(sat).all() or sum(sat) <= 0:
            why = "no n/iso/cyclo split"
        else:
            shares.append(np.array(sat) / sum(sat))
            used.append(fid)
            continue
        log.warning("nominal split: skipping %s (%s)", fid, why)
        skipped.append((fid, why))
    if not shares:
        raise ValueError("no GC jet fuel with a complete n/iso/cyclo split")
    m = np.mean(shares, axis=0)
    return (float(m[0]), float(m[1]), float(m[2])), used, skipped


def nominal_sat_split(jet_class_csv: Path) -> tuple[float, float, float]:
    """Mean n / iso / cyclo share of the saturates over the GC jet fuels (a prior, not a fit to CRC)."""
    return nominal_sat_split_detail(jet_class_csv)[0]


def diaromatic_mass(row: Mapping) -> float:
    """Mass fraction of the diaromatic lump (aromatic mass x D/A), 0 when the row cannot be reconstructed."""
    A = float(row.get("aromatics_vol_pct", np.nan))
    if not np.isfinite(A) or A <= 0:
        return 0.0
    get = lambda k: float(row[k]) if k in row and np.isfinite(row[k]) else 0.0
    D, O = min(get("diaromatics_vol_pct"), A), get("olefins_vol_pct")
    ma, mo, ms = A * RHO_CLASS["aromatic"], O * RHO_CLASS["olefin"], max(100.0 - A - O, 0.0) * RHO_CLASS["sat"]
    return ma / (ma + mo + ms) * D / A


def reconstruct(row: Mapping, sat_split, reps) -> tuple[pd.DataFrame, float]:
    cols = ["inchikey", "smiles", "fraction", "kind"]
    empty = pd.DataFrame(columns=cols)
    cd = carbon_distribution(row)
    A = float(row.get("aromatics_vol_pct", np.nan))
    if cd is None or not np.isfinite(A):
        log.warning("reconstruct: skipped (%s)", "fewer than 4 distillation points" if cd is None else "no aromatics value")
        return empty, 1.0
    get = lambda k: float(row[k]) if k in row and np.isfinite(row[k]) else 0.0
    O, D = get("olefins_vol_pct"), min(get("diaromatics_vol_pct"), A)
    S = max(100.0 - A - O, 0.0)
    m = {"aromatic": A * RHO_CLASS["aromatic"], "alkene_alkyne": O * RHO_CLASS["olefin"], "sat": S * RHO_CLASS["sat"]}
    tot = sum(m.values())
    m = {k: v / tot for k, v in m.items()}

    items: dict[tuple[str, int], float] = {}

    def add(kind, n, mass):
        items[(kind, n)] = items.get((kind, n), 0.0) + mass

    for kind, share in zip(("n_paraffin", "iso_paraffin", "cycloalkane"), sat_split):
        for n, f in cd.items():
            add(kind, n, m["sat"] * share * f)
    for n, f in cd.items():
        add("alkene_alkyne", n, m["alkene_alkyne"] * f)
        add("aromatic_mono", n, m["aromatic"] * (1 - D / A if A > 0 else 1.0) * f)

    out: dict[tuple[int, str], float] = {}
    fallback = 0.0
    for (kind, n), mass in items.items():
        if mass <= 0:
            continue
        i = reps.get(kind, n)
        if i is None:
            fallback += mass
        else:
            out[(i, kind)] = out.get((i, kind), 0.0) + mass
    if D > 0:
        mass = m["aromatic"] * D / A
        i = reps.get("aromatic_di", 11, max_gap=99)
        if i is None:
            fallback += mass
        else:
            out[(i, "aromatic_di")] = out.get((i, "aromatic_di"), 0.0) + mass
    rows = [dict(inchikey=reps.ik[i], smiles=reps.smiles[i], fraction=v, kind=k) for (i, k), v in out.items()]
    return (pd.DataFrame(rows, columns=cols) if rows else empty), float(fallback)
