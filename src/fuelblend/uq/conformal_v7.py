"""v7 Task 19: weighted, family-conditional (Mondrian) one-sided conformal margins and covariate-shift weights.

Row-level calibration (one score per row), complementary to `conformal.py` (group-max scores). The finite-sample rank is
ceil((n + 1)(1 - alpha)); with weights, the test point carries weight 1 (weights are normalised to mean 1) at +inf.
"""
from __future__ import annotations

import warnings

import numpy as np

W_LO, W_HI = 0.05, 20.0


def one_sided_margin(resid, alpha: float = 0.10, weights=None) -> float:
    """Weighted (1 - alpha) quantile of the unsafe-side scores `resid` with the finite-sample correction; +inf if too few rows."""
    r = np.asarray(resid, float)
    r = r[np.isfinite(r)] if weights is None else r
    n = r.size
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    if n == 0:
        return float("inf")
    if weights is None:
        k = int(np.ceil((n + 1) * (1 - alpha)))
        if k > n:
            warnings.warn(f"margin undefined: n={n} too few for alpha={alpha}; returning inf", RuntimeWarning, stacklevel=2)
            return float("inf")
        return float(np.sort(r)[k - 1])
    w = np.asarray(weights, float)
    ok = np.isfinite(r) & np.isfinite(w)
    r, w = r[ok], w[ok]
    w = w / w.mean()
    o = np.argsort(r)
    r, w = r[o], w[o]
    cum = np.cumsum(w) / (w.sum() + 1.0)                 # test point: weight 1 at +inf
    j = np.searchsorted(cum, 1 - alpha - 1e-12)
    if j >= len(r):
        warnings.warn("weighted margin undefined (not enough effective rows); returning inf", RuntimeWarning, stacklevel=2)
        return float("inf")
    return float(r[j])


def group_conditional_margins(resid, groups, alpha: float = 0.10, min_n: int = 10) -> dict:
    """Mondrian margins per family; a family with fewer than `min_n` calibration rows gets the global margin. Key '__global__' holds the global one."""
    r, g = np.asarray(resid, float), np.asarray(groups)
    out = {"__global__": one_sided_margin(r, alpha)}
    for f in sorted(set(g)):
        m = g == f
        out[str(f)] = one_sided_margin(r[m], alpha) if m.sum() >= min_n else out["__global__"]
    return out


def shift_weights(X_cal, X_test) -> np.ndarray:
    """Density-ratio weights for the calibration rows (test over calibration, logistic regression), clipped to [0.05, 20]."""
    from sklearn.linear_model import LogisticRegression
    Xc, Xt = np.asarray(X_cal, float), np.asarray(X_test, float)
    mu, sd = Xc.mean(0), np.where(Xc.std(0) < 1e-9, 1.0, Xc.std(0))
    Z = np.vstack([(Xc - mu) / sd, (Xt - mu) / sd])
    lab = np.r_[np.zeros(len(Xc)), np.ones(len(Xt))]
    clf = LogisticRegression(C=1.0, max_iter=1000).fit(Z, lab)
    p = clf.predict_proba((Xc - mu) / sd)[:, 1]
    w = (p / np.clip(1 - p, 1e-9, None)) * (len(Xc) / len(Xt))
    return np.clip(w, W_LO, W_HI)
