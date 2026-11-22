"""v7 Task 15: heavy-tail composition features, ARD Matern GP and a Mahalanobis novelty flag for fuel freezing point."""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem

FEATS = ["nC14plus_mass", "naphthalene_mass", "dicyclo_mass", "tm_max_present", "tm_p90_mass"]
_NAPH = Chem.MolFromSmarts("c1ccc2ccccc2c1")
X_MIN = 1e-4


def _species(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    nC = sum(a.GetAtomicNum() == 6 for a in m.GetAtoms())
    arom = any(a.GetIsAromatic() for a in m.GetAtoms())
    return nC, m.HasSubstructMatch(_NAPH), (not arom) and m.GetRingInfo().NumRings() >= 2


def heavy_tail_features(comp_df: pd.DataFrame) -> pd.DataFrame:
    """comp_df columns: fuel, smiles, frac (mass fraction, renormalised per fuel), tm_K (NaN if unknown). One row per fuel, columns FEATS."""
    out = {}
    for fuel, g in comp_df.groupby("fuel", sort=False):
        w = g.frac.to_numpy(float)
        w = w / w.sum() if w.sum() > 0 else w
        sp = [_species(s) for s in g.smiles]
        nC14 = float(sum(wi for wi, s in zip(w, sp) if s[0] >= 14))
        naph = float(sum(wi for wi, s in zip(w, sp) if s[1]))
        dic = float(sum(wi for wi, s in zip(w, sp) if s[2]))
        tm = g.tm_K.to_numpy(float)
        ok = np.isfinite(tm) & (w > X_MIN)
        if ok.any():
            tmax = float(tm[ok].max())
            o = np.argsort(-tm[ok])
            cum = np.cumsum(w[ok][o]) / w[ok].sum()
            p90 = float(tm[ok][o][np.searchsorted(cum, 0.10)])
        else:
            tmax = p90 = np.nan
        out[fuel] = [nC14, naph, dic, tmax, p90]
    return pd.DataFrame.from_dict(out, orient="index", columns=FEATS)


def standardise(X_tr):
    X_tr = np.asarray(X_tr, float)
    mu, sd = X_tr.mean(0), X_tr.std(0)
    return mu, np.where(sd < 1e-9, 1.0, sd)


def fit_gp(X, y, restarts: int = 5, seed: int = 0):
    """GP on already standardised X: ConstantKernel * Matern(nu=2.5, ARD) + WhiteKernel, normalised target."""
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    X = np.asarray(X, float)
    k = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(length_scale=np.ones(X.shape[1]), length_scale_bounds=(1e-1, 1e2), nu=2.5) \
        + WhiteKernel(0.1, (1e-4, 1.0))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=restarts, random_state=seed).fit(X, np.asarray(y, float))


def novelty_score(X_tr, X_te) -> np.ndarray:
    """Mahalanobis distance of each X_te row to the training cloud, in the standardised feature space of X_tr (ridge-regularised covariance)."""
    mu, sd = standardise(X_tr)
    Z = (np.asarray(X_tr, float) - mu) / sd
    P = np.linalg.inv(np.cov(Z.T, bias=True).reshape(Z.shape[1], Z.shape[1]) + 1e-3 * np.eye(Z.shape[1]))
    D = (np.asarray(X_te, float) - mu) / sd
    return np.sqrt(np.einsum("ij,jk,ik->i", D, P, D))


def novelty_threshold(X_tr, q: float = 95.0) -> float:
    """q-th percentile of the leave-one-out novelty scores of the training rows."""
    X_tr = np.asarray(X_tr, float)
    s = [novelty_score(np.delete(X_tr, i, 0), X_tr[i:i + 1])[0] for i in range(len(X_tr))]
    return float(np.percentile(s, q))
