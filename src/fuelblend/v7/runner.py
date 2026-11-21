"""Shared v7 runner: Track B pool loader (folds computed before locked-label masking), result writer, paired bootstrap."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold

from fuelblend.v7 import guard

ROOT = guard.ROOT
POOL_CSV = ROOT / "results/v6/gnn/pool.csv"


def _group(smiles: str) -> str:
    """Generic Murcko scaffold; acyclic molecules grouped by carbon number (as scripts/v6/gnn/train_eval.py)."""
    m = Chem.MolFromSmiles(smiles)
    if m.GetRingInfo().NumRings() == 0:
        return f"acyc_C{sum(a.GetAtomicNum() == 6 for a in m.GetAtoms())}"
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))


def scaffold_folds(df: pd.DataFrame, k: int = 5) -> np.ndarray:
    """Fold id per row of a compound frame with columns smiles and grp (computed here if grp is absent)."""
    grp = df["grp"] if "grp" in df else df["smiles"].map(_group)
    out = np.full(len(df), -1)
    for f, (_, te) in enumerate(GroupKFold(k).split(df, groups=grp)):
        out[te] = f
    return out


def load_track_b_pool(k: int = 5) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    """Wide compound frame (inchikey, smiles, is_hc, grp, one column per task), fold ids, task names.

    Folds are computed on the full pool first (identical to v6 Track B), then YSI labels of the ysi_holdout lock are set to
    NaN, so folds are never recomputed on a reduced pool. df.attrs['n_locked_ysi'] records the count."""
    pool = pd.read_csv(POOL_CSV)
    comp = pool.drop_duplicates("inchikey")[["inchikey", "smiles", "is_hc"]].reset_index(drop=True)
    comp["grp"] = comp.smiles.map(_group)
    tasks = sorted(pool.target.unique())
    wide = pd.pivot_table(pool, index="inchikey", columns="target", values="value").reindex(comp.inchikey)[tasks]
    for t in tasks:
        comp[t] = wide[t].values
    folds = scaffold_folds(comp, k)
    comp, n = guard.mask_locked(comp, "inchikey", locks=("ysi_holdout",), label_cols=["ysi"])
    comp.attrs["n_locked_ysi"] = n
    return comp, folds, tasks


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return ""


def write_results(out_dir, *, summary: pd.DataFrame, oof: pd.DataFrame, meta: dict, drops: dict | None = None) -> None:
    """Write summary.csv/cv_summary.csv, oof.csv and meta.json (device info, git commit, guard drop counts)."""
    from fuelblend.ml.dmpnn import device_info
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / "cv_summary.csv", index=False)
    oof.to_csv(out / "oof.csv", index=False)
    m = {**meta, **device_info(), "git_commit": _git_commit(), "guard_drops": drops or {}}
    (out / "meta.json").write_text(json.dumps(m, indent=2, default=str), encoding="utf-8")


def paired_boot_gain(err_a, err_b, groups, n: int = 4000, seed: int = 0) -> tuple[float, float]:
    """90 % group-bootstrap CI of gain = mean|err_b| - mean|err_a| (positive: a, the candidate, is better).

    Same resampling scheme as octane_blend_v4.cluster_boot (sum over resampled groups / rows), 5th and 95th percentiles."""
    g = np.asarray(groups)
    d = np.abs(np.asarray(err_b, float)) - np.abs(np.asarray(err_a, float))
    ok = np.isfinite(d)
    g, d = g[ok], d[ok]
    ug = np.unique(g)
    S = np.array([d[g == u].sum() for u in ug])
    N = np.array([(g == u).sum() for u in ug], float)
    idx = np.random.default_rng(seed).integers(0, len(ug), (n, len(ug)))
    est = S[idx].sum(1) / N[idx].sum(1)
    return float(np.percentile(est, 5)), float(np.percentile(est, 95))


def fit_with_oom_retry(fit_fn, **kw):
    """Call fit_fn(**kw); on CUDA OOM halve kw['bs'] once and retry with oom_retry=1 recorded in the second call's kwargs.
    A second OOM is re-raised."""
    import torch
    try:
        return fit_fn(**kw)
    except torch.cuda.OutOfMemoryError:
        torch.cuda.empty_cache()
        kw = {**kw, "bs": max(1, int(kw.get("bs", 32)) // 2), "oom_retry": 1}
        return fit_fn(**kw)
