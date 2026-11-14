"""Carbon-number interpolation of Antoine constants in a homologous series (v5 delta_ws3_pure_gaps, Part A).

Convention: log10(P/kPa) = A - B / (T/K + C).

* ``grid_interp`` (designated): for each T of a grid, log10 P of every knot is interpolated in carbon number with a monotone PCHIP, then
  the target's (A, B, C) is refit to the interpolated grid.
* ``param_interp`` (comparator): PCHIP of A, B and C separately (ill-conditioned: A, B, C are strongly correlated).
Knots are never modified; the returned constants carry ``antoine_source = "interpolated"`` when written by the delta script.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

T_GRID = np.arange(293.15, 453.15 + 1e-9, 5.0)
SOURCE = "interpolated"


def log10p(abc: Sequence[float], T) -> np.ndarray:
    A, B, C = (float(v) for v in abc)
    return A - B / (np.asarray(T, float) + C)


def refit_abc(T, logp, x0=(6.5, 2000.0, -80.0)) -> tuple[float, float, float]:
    """Least-squares (A, B, C) to log10 P(T); C bounded to keep T + C > 0 on the data."""
    T = np.asarray(T, float)
    logp = np.asarray(logp, float)
    cmin = -float(T.min()) + 20.0
    lo, hi = [-20.0, 1.0, cmin], [30.0, 2.0e4, 300.0]
    x0 = np.clip(np.asarray(x0, float), np.array(lo) + 1e-6, np.array(hi) - 1e-6)
    r = least_squares(lambda p: log10p(p, T) - logp, x0, bounds=(lo, hi), x_scale=[1.0, 1000.0, 50.0], max_nfev=20000)
    return tuple(float(v) for v in r.x)


def grid_interp(knots: Mapping[int, Sequence[float]], n_target: int, T_grid=T_GRID) -> tuple[tuple[float, float, float], np.ndarray]:
    """(A, B, C) of carbon number ``n_target`` from the knots {n: (A, B, C)}; also returns the interpolated log10 P grid."""
    ns = np.array(sorted(knots), float)
    if n_target in knots:
        raise ValueError("target is a knot; knots are never overwritten")
    if not (ns.min() < n_target < ns.max()):
        raise ValueError("target outside the knot range (no extrapolation)")
    M = np.array([log10p(knots[int(n)], T_grid) for n in ns])          # knots x T
    lp = np.array([PchipInterpolator(ns, M[:, j])(n_target) for j in range(len(T_grid))], float)
    i = int(np.argmin(np.abs(ns - n_target)))
    return refit_abc(T_grid, lp, x0=knots[int(ns[i])]), lp


def param_interp(knots: Mapping[int, Sequence[float]], n_target: int) -> tuple[float, float, float]:
    ns = np.array(sorted(knots), float)
    if n_target in knots:
        raise ValueError("target is a knot; knots are never overwritten")
    P = np.array([knots[int(n)] for n in ns], float)
    return tuple(float(PchipInterpolator(ns, P[:, k])(n_target)) for k in range(3))


def rms_log10(abc, T, P_kPa) -> float:
    T, P = np.asarray(T, float), np.asarray(P_kPa, float)
    ok = P > 0
    if not ok.any():
        return float("nan")
    return float(np.sqrt(np.mean((log10p(abc, T[ok]) - np.log10(P[ok])) ** 2)))
