"""Learned pair excess-volume correction to ideal-volume mixture density (v5 delta_ml_density).

V_ML = V_ideal + sum_{i<j} x_i x_j g(d_i, d_j, T)   [cm3 per mol of mixture],   rho_ML = sum_i x_i M_i / V_ML   [kg/m3].
g is a LightGBM regressor on pair-symmetric features (mean and |difference| of RDKit descriptors and the pure molar volume, unordered
structural class pair, T).  No composition inside g: the correction is zero for a pure component, order invariant, defined for any n.

Public API: descriptors(smiles), excess_volume(smiles_list, x_mole, T_K, rho_i=None), density_ml(...), domain_status(...).
Pure densities rho_i (kg/m3) default to the frozen pure grid (read only); training used the tier-1 ThermoML fits.
"""
from __future__ import annotations

import json
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MODEL_DIR = ROOT / "models" / "v5" / "ml_density_pair"
CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]
CLASS_PAIRS = ["|".join(p) for p in combinations_with_replacement(sorted(CLASSES), 2)]
STATIC = ["n_carbon", "mw", "n_rings", "n_arom_rings", "n_ch3", "n_quat_c"]
BASE_FEATS = STATIC + ["v_molar"]
FEATURES = [f"{p}_{b}" for b in BASE_FEATS for p in ("mean", "absdiff")] + ["class_pair", "T"]
LGB_PARAMS = dict(n_estimators=400, learning_rate=0.05, num_leaves=15, min_child_samples=40, subsample=0.8, subsample_freq=1,
                  colsample_bytree=0.8, reg_lambda=1.0, random_state=42, deterministic=True, force_row_wise=True, n_jobs=1, verbose=-1)
MIN_SYSTEMS_PER_PAIR = 5
_MODEL = None


def descriptors_from_mol(m):
    """-> (static vector aligned to STATIC, structural class) for an RDKit molecule (hydrocarbon expected)."""
    from rdkit.Chem import Descriptors, rdMolDescriptors
    heavy = list(m.GetAtoms())
    ch3 = sum(1 for a in heavy if a.GetAtomicNum() == 6 and a.GetDegree() == 1 and a.GetTotalNumHs() == 3)
    quat = sum(1 for a in heavy if a.GetAtomicNum() == 6 and a.GetDegree() == 4 and all(n.GetAtomicNum() == 6 for n in a.GetNeighbors()))
    nr, na = rdMolDescriptors.CalcNumRings(m), rdMolDescriptors.CalcNumAromaticRings(m)
    if any(a.GetIsAromatic() for a in heavy):
        cls = "aromatic"
    elif any(b.GetBondTypeAsDouble() > 1.0 for b in m.GetBonds()):
        cls = "alkene_alkyne"
    elif nr > 0:
        cls = "cycloalkane"
    elif all(a.GetDegree() <= 2 for a in heavy):
        cls = "n_paraffin"
    else:
        cls = "iso_paraffin"
    nc = sum(1 for a in heavy if a.GetAtomicNum() == 6)
    return np.array([nc, Descriptors.MolWt(m), nr, na, ch3, quat], float), cls


def descriptors(smiles: str):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"bad SMILES {smiles}")
    return descriptors_from_mol(m)


def pair_frame(sa, sb, cls_a, cls_b, rho_a, rho_b, T) -> pd.DataFrame:
    """Vectorised: sa, sb (k, 6) static descriptors, cls_*, rho_* (kg/m3), T (K) arrays of length k."""
    sa, sb = np.atleast_2d(sa), np.atleast_2d(sb)
    va = sa[:, 1] / np.asarray(rho_a, float) * 1000.0
    vb = sb[:, 1] / np.asarray(rho_b, float) * 1000.0
    da, db = np.column_stack([sa, va]), np.column_stack([sb, vb])
    cols = {}
    for k, b in enumerate(BASE_FEATS):
        cols[f"mean_{b}"] = 0.5 * (da[:, k] + db[:, k])
        cols[f"absdiff_{b}"] = np.abs(da[:, k] - db[:, k])
    cp = ["|".join(sorted((a, b))) for a, b in zip(cls_a, cls_b)]
    cols["class_pair"] = np.array([CLASS_PAIRS.index(c) if c in CLASS_PAIRS else -1 for c in cp])
    cols["T"] = np.broadcast_to(np.asarray(T, float), (len(sa),)).copy()
    return pd.DataFrame(cols)[FEATURES]


def load_model():
    global _MODEL
    if _MODEL is None:
        import lightgbm as lgb
        _MODEL = lgb.Booster(model_file=str(MODEL_DIR / "model.txt"))
    return _MODEL


def _grid_rho(smiles_list, T_K):
    from fuelblend.ml.viscosity_pair import context
    ctx = context()
    idx = [ctx.index_of(s) for s in smiles_list]
    if not all(i is not None and bool(ctx.prov.valid[i]) for i in idx):
        raise ValueError("compound without frozen-grid density; pass rho_i explicitly")
    return np.asarray(ctx.prov.lookup(np.array(idx), np.full(len(idx), T_K))[1], float)


def _prep(smiles_list, x_mole, T_K, rho_i):
    d = [descriptors(s) for s in smiles_list]
    S, cls = np.array([a for a, _ in d]), [c for _, c in d]
    x = np.asarray(x_mole, float)
    x = x / x.sum()
    rho = _grid_rho(smiles_list, T_K) if rho_i is None else np.asarray(rho_i, float)
    return S, cls, x, rho


def interaction_matrix(S, cls, rho, T_K) -> np.ndarray:
    n = len(S)
    ii, jj = np.triu_indices(n, k=1)
    G = np.zeros((n, n))
    if len(ii):
        g = np.asarray(load_model().predict(pair_frame(S[ii], S[jj], [cls[i] for i in ii], [cls[j] for j in jj], rho[ii], rho[jj], T_K)), float)
        G[ii, jj] = g
        G[jj, ii] = g
    return G


def excess_volume(smiles_list, x_mole, T_K: float, rho_i=None) -> float:
    """sum_{i<j} x_i x_j g_ij in cm3/mol.  Zero for a pure component."""
    S, cls, x, rho = _prep(smiles_list, x_mole, T_K, rho_i)
    return float(0.5 * x @ interaction_matrix(S, cls, rho, T_K) @ x)


def density_ml(smiles_list, x_mole, T_K: float, rho_i=None) -> dict:
    S, cls, x, rho = _prep(smiles_list, x_mole, T_K, rho_i)
    M = S[:, 1]
    v_ideal = float(np.sum(x * M / rho) * 1000.0)
    ve = float(0.5 * x @ interaction_matrix(S, cls, rho, T_K) @ x)
    mm = float(np.sum(x * M))
    return dict(rho_ideal=mm / v_ideal * 1000.0, rho_ml=mm / (v_ideal + ve) * 1000.0, excess_volume_cm3_mol=ve)


def domain_status(smiles_list) -> dict:
    meta = json.loads((MODEL_DIR / "meta.json").read_text())
    try:
        d = [descriptors(s) for s in smiles_list]
    except ValueError:
        return dict(in_domain=False, reason="unparseable SMILES")
    cls = [c for _, c in d]
    counts = meta["class_pair_training_systems"]
    for a in range(len(cls)):
        for b in range(a + 1, len(cls)):
            if counts.get("|".join(sorted((cls[a], cls[b]))), 0) < MIN_SYSTEMS_PER_PAIR:
                return dict(in_domain=False, reason="class pair under-represented in training")
    return dict(in_domain=True, reason="ok")
