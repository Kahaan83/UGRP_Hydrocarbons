"""v5 BASE: per-temperature cold-viscosity offset K1 (moved unchanged from scripts/v5/c1_cold_viscosity_correction.py).

Offset = mean ln(meas/pred) over training rows at each temperature; corrected prediction = pred * exp(offset)."""
from __future__ import annotations

import numpy as np


def fit_k1(T, resid) -> dict:
    """K1: per-temperature constant ln offset = mean(ln(meas/pred)) over training rows.  T keyed by round(T, 2)."""
    T, r = np.round(np.asarray(T, float), 2), np.asarray(resid, float)
    ok = np.isfinite(r) & np.isfinite(T)
    return {float(t): float(r[ok & (T == t)].mean()) for t in np.unique(T[ok])}


def apply_k1(k1: dict, T) -> np.ndarray:
    """Offset at each T; NaN where the temperature was not in the training data (K1 is not usable there)."""
    return np.array([k1.get(round(float(t), 2), np.nan) for t in np.atleast_1d(T)])
