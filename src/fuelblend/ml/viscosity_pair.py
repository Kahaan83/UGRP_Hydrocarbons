"""Learned pair-interaction correction to the Refutas mixture viscosity (v5 delta_ml_viscosity).

ln eta_ML = ln eta_Refutas + sum_{i<j} x_i x_j g(d_i, d_j, T), x = mole fractions.  g is a LightGBM regressor on pair-symmetric features
(mean and |difference| of per-compound descriptors, unordered class pair, 1000/T).  There is no composition feature inside g, so the
correction is zero for a pure component, invariant to component order, and defined for any number of components.

Public API: pair_features(smiles_i, smiles_j, T_K), interaction_matrix(smiles_list, T_K), ln_correction(smiles_list, x_mole, T_K),
domain_status(smiles_list, T_K).  Heavy objects (registry, frozen pure grid, booster) load lazily on first use.
"""
from __future__ import annotations

import json
from functools import lru_cache
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MODEL_DIR = ROOT / "models" / "v5" / "ml_viscosity_pair"
COMPOUNDS = ROOT / "data" / "unified" / "compounds.csv"
GRID = ROOT / "data" / "cache" / "pure_grid.npz"
CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]
CLASS_PAIRS = ["|".join(p) for p in combinations_with_replacement(sorted(CLASSES), 2)]
STATIC = ["n_carbon", "mw", "n_rings", "n_arom_rings", "n_ch3", "n_quat_c"]
DYNAMIC = ["ln_eta", "v_molar"]
BASE_FEATS = STATIC + DYNAMIC
FEATURES = [f"{p}_{b}" for b in BASE_FEATS for p in ("mean", "absdiff")] + ["class_pair", "inv_T"]
LGB_PARAMS = dict(n_estimators=400, learning_rate=0.05, num_leaves=15, min_child_samples=40, subsample=0.8, subsample_freq=1,
                  colsample_bytree=0.8, reg_lambda=1.0, random_state=42, deterministic=True, force_row_wise=True, n_jobs=1, verbose=-1)
MIN_SYSTEMS_PER_PAIR = 5
COLD_GATE_K = 283.15           # correction defined for T >= this unless meta says the cold gate passed


class Context:
    """Registry + frozen pure-property grid + static descriptors for every registry compound."""

    def __init__(self):
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
        from fuelblend.mixture.cold_flow import compound_meta
        self.reg = CompoundRegistry(COMPOUNDS)
        _, self.cls, self.nc = compound_meta(COMPOUNDS)
        self.prov = PurePropertyProvider(self.reg, GRID, verbose=False)
        n = len(self.reg.ik)
        st = np.full((n, len(STATIC)), np.nan)
        for i, s in enumerate(self.reg.smiles):
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            ch3 = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6 and a.GetDegree() == 1 and a.GetTotalNumHs() == 3)
            quat = sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 6 and a.GetDegree() == 4
                       and all(nb.GetAtomicNum() == 6 for nb in a.GetNeighbors()))
            st[i] = [self.nc[i], self.reg.mw[i], rdMolDescriptors.CalcNumRings(m), rdMolDescriptors.CalcNumAromaticRings(m), ch3, quat]
        self.static = st
        self._Chem = Chem

    def index_of(self, smiles: str):
        m = self._Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        return self.reg.ik2i.get(self._Chem.MolToInchiKey(m))

    def usable(self, i) -> bool:
        return i is not None and bool(self.prov.valid[i]) and bool(self.reg.usable[i]) and self.cls[i] in CLASSES and np.isfinite(self.static[i]).all()


_CTX: Context | None = None
_MODEL = None


def context() -> Context:
    global _CTX
    if _CTX is None:
        _CTX = Context()
    return _CTX


def pair_feature_frame(ia, ib, T_K, ctx: Context | None = None) -> pd.DataFrame:
    """Vectorised features for compound-index pairs ia[k], ib[k] at T_K[k] (symmetric in a <-> b)."""
    ctx = ctx or context()
    ia, ib = np.asarray(ia, np.int64), np.asarray(ib, np.int64)
    T = np.broadcast_to(np.asarray(T_K, float), ia.shape)
    cols = {}
    sa, sb = ctx.static[ia], ctx.static[ib]
    le_a, rho_a = ctx.prov.lookup(ia, T)
    le_b, rho_b = ctx.prov.lookup(ib, T)
    da = np.column_stack([sa, le_a, ctx.reg.mw[ia] / rho_a * 1000.0])
    db = np.column_stack([sb, le_b, ctx.reg.mw[ib] / rho_b * 1000.0])
    for k, b in enumerate(BASE_FEATS):
        cols[f"mean_{b}"] = 0.5 * (da[:, k] + db[:, k])
        cols[f"absdiff_{b}"] = np.abs(da[:, k] - db[:, k])
    cp = ["|".join(sorted((ctx.cls[a], ctx.cls[b]))) for a, b in zip(ia, ib)]
    cols["class_pair"] = np.array([CLASS_PAIRS.index(c) if c in CLASS_PAIRS else -1 for c in cp])
    cols["inv_T"] = 1000.0 / T
    return pd.DataFrame(cols)[FEATURES]


def load_model():
    global _MODEL
    if _MODEL is None:
        import lightgbm as lgb
        _MODEL = lgb.Booster(model_file=str(MODEL_DIR / "model.txt"))
    return _MODEL


def predict_g(frame: pd.DataFrame) -> np.ndarray:
    return np.asarray(load_model().predict(frame[FEATURES]), float)


def pair_features(smiles_i: str, smiles_j: str, T_K: float) -> pd.DataFrame:
    ctx = context()
    ia, ib = ctx.index_of(smiles_i), ctx.index_of(smiles_j)
    if not (ctx.usable(ia) and ctx.usable(ib)):
        raise ValueError("compound without frozen-grid / registry data")
    return pair_feature_frame([ia], [ib], [T_K], ctx)


def interaction_matrix(smiles_list, T_K: float) -> np.ndarray:
    """Symmetric (n, n) matrix of g, zero diagonal.  Memoised on (members, T): g does not depend on the fractions, so an optimiser
    that re-evaluates many recipes of one pool at one temperature pays for the LightGBM call once (a copy is returned)."""
    return _interaction_matrix_cached(tuple(smiles_list), float(T_K)).copy()


@lru_cache(maxsize=4096)
def _interaction_matrix_cached(smiles_list: tuple, T_K: float) -> np.ndarray:
    ctx = context()
    idx = [ctx.index_of(s) for s in smiles_list]
    if not all(ctx.usable(i) for i in idx):
        raise ValueError("compound without frozen-grid / registry data")
    n = len(idx)
    ii, jj = np.triu_indices(n, k=1)
    G = np.zeros((n, n))
    if len(ii):
        g = predict_g(pair_feature_frame(np.array(idx)[ii], np.array(idx)[jj], np.full(len(ii), T_K), ctx))
        G[ii, jj] = g
        G[jj, ii] = g
    return G


def ln_correction(smiles_list, x_mole, T_K: float) -> float:
    """sum_{i<j} x_i x_j g_ij.  Zero for a pure component or when only one component has weight."""
    x = np.asarray(x_mole, float)
    x = x / x.sum()
    G = interaction_matrix(smiles_list, T_K)
    return float(0.5 * x @ G @ x)


@lru_cache(maxsize=1)
def _meta() -> dict:
    return json.loads((MODEL_DIR / "meta.json").read_text())


def domain_status(smiles_list, T_K: float) -> dict:
    """in_domain per the declared adoption rule: every component has own grid data, every class pair has >= 5 training systems,
    T inside the cold gate, known composition (the caller must pass False for reconstructed real fuels)."""
    ctx = context()
    meta = _meta()
    idx = [ctx.index_of(s) for s in smiles_list]
    if not all(ctx.usable(i) for i in idx):
        return dict(in_domain=False, reason="component without grid data")
    cls = [ctx.cls[i] for i in idx]
    counts = meta["class_pair_training_systems"]
    for a in range(len(cls)):
        for b in range(a + 1, len(cls)):
            if counts.get("|".join(sorted((cls[a], cls[b]))), 0) < MIN_SYSTEMS_PER_PAIR:
                return dict(in_domain=False, reason="class pair under-represented in training")
    t_min = meta.get("T_min_K", COLD_GATE_K)
    if T_K < t_min:
        return dict(in_domain=False, reason=f"T below the validated gate {t_min} K")
    return dict(in_domain=True, reason="ok")
