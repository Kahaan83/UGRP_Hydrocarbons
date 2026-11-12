"""v5 Task 13: low-capacity learned residual on top of physics, from measured composition bins.

bin_features: one fuel's bins (hc_class, n_carbon, mass_pct; hc_class a SubReps kind or 'no_rep:<label>') -> fixed feature vector.
fit_predict_residual: ridge on features standardised with the TRAINING rows only; alpha by inner leave-one-group-out on the training rows only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

CLASSES = ("n_paraffin", "iso_paraffin", "cyclo_mono", "cyclo_di", "alkylbenzene", "indan_tetralin", "aromatic_di", "alkene_alkyne",
           "no_rep_aromatic", "no_rep_other")
GROUPS = {"n": ("n_paraffin",), "iso": ("iso_paraffin",), "cyclo": ("cyclo_mono", "cyclo_di"),
          "aro": ("alkylbenzene", "indan_tetralin", "aromatic_di", "no_rep_aromatic")}
FEATURE_NAMES: list[str] = ([f"w_{c}" for c in CLASSES] + ["cn_mean", "cn_sd"]
                            + [f"cn_{s}_{g}" for g in GROUPS for s in ("mean", "sd")] + ["w_below_C8", "w_above_C15"])
_ARO_WORDS = ("arom", "naph", "fluor", "inden", "benzo", "phenan", "acenaph", "tetralin", "indan")


def _cls(h: str) -> str:
    if not h.startswith("no_rep"):
        return h
    return "no_rep_aromatic" if any(w in h.lower() for w in _ARO_WORDS) else "no_rep_other"


def _ms(cn, w):
    if w.sum() <= 0:
        return 0.0, 0.0
    m = float(np.average(cn, weights=w))
    return m, float(np.sqrt(np.average((cn - m) ** 2, weights=w)))


def bin_features(bins: pd.DataFrame) -> np.ndarray:
    c = np.array([_cls(str(h)) for h in bins.hc_class])
    cn = bins.n_carbon.to_numpy(float)
    w = bins.mass_pct.to_numpy(float)
    w = w / w.sum()
    f = [float(w[c == k].sum()) for k in CLASSES]
    f += list(_ms(cn, w))
    for g in GROUPS.values():
        m = np.isin(c, g)
        f += list(_ms(cn[m], w[m]))
    f += [float(w[cn < 8].sum()), float(w[cn > 15].sum())]
    return np.array(f, float)


def _ridge(X, y, alpha):
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    ym = y.mean()
    beta = np.linalg.solve(Z.T @ Z + alpha * np.eye(Z.shape[1]), Z.T @ (y - ym))
    return lambda Xn: ym + ((Xn - mu) / sd) @ beta


def fit_predict_residual(X_train, r_train, g_train, X_test, alpha_grid=(0.1, 1, 10, 100)) -> np.ndarray:
    X_train, r_train, X_test = np.asarray(X_train, float), np.asarray(r_train, float), np.asarray(X_test, float)
    g = np.asarray(g_train).astype(str)
    best, best_mse = alpha_grid[0], np.inf
    for a in alpha_grid:
        err = []
        for gg in np.unique(g):
            tr, te = g != gg, g == gg
            if tr.sum() < 2:
                continue
            err.append((_ridge(X_train[tr], r_train[tr], a)(X_train[te]) - r_train[te]) ** 2)
        mse = float(np.concatenate(err).mean()) if err else np.inf
        if mse < best_mse:
            best, best_mse = a, mse
    return _ridge(X_train, r_train, best)(X_test)
