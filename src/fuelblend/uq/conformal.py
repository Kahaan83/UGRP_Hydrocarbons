"""Group-level split-conformal margins computed from residuals (v4 Task 10).

Residual convention everywhere: ``residual = predicted - measured``.

* ``direction="over"``  : the unsafe side is over-prediction (score = residual), e.g. flash point, RON/MON,
  heat of combustion, smoke point against a minimum limit.
* ``direction="under"`` : the unsafe side is under-prediction (score = -residual), e.g. freezing point against a
  maximum limit.

The independent unit is the *group* (fuel / system / blendstock), never the row: the score of a group is the
maximum unsafe-side residual over its rows, and the margin is the finite-sample conformal quantile
(k-th smallest of n group scores, k = ceil((n + 1) * (1 - alpha))).  If k > n there is no valid margin at that
alpha: the function returns ``inf`` and emits a warning (never a small optimistic number).
"""
from __future__ import annotations

import math
import warnings

import numpy as np

_DIRECTIONS = ("over", "under")


def conformal_rank(n_groups: int, alpha: float) -> int:
    """k = ceil((n + 1)(1 - alpha)); rounding guards against float fuzz such as 9.000000000000002."""
    return int(math.ceil(round((n_groups + 1) * (1.0 - alpha), 9)))


def _group_max(scores: np.ndarray, groups: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=float)
    groups = np.asarray(groups)
    if scores.shape[0] != groups.shape[0]:
        raise ValueError("residuals and groups must have the same length")
    ok = np.isfinite(scores)
    scores, groups = scores[ok], groups[ok]
    if scores.size == 0:
        return np.array([], dtype=float)
    gs = groups.astype(str)
    order = np.argsort(gs, kind="stable")
    g = gs[order]
    s = scores[order]
    starts = np.r_[0, np.flatnonzero(g[1:] != g[:-1]) + 1]
    return np.maximum.reduceat(s, starts)


def group_scores(residuals, groups, direction: str = "over", two_sided: bool = False) -> np.ndarray:
    """Per-group maximum unsafe-side score (non-finite residuals are dropped)."""
    r = np.asarray(residuals, dtype=float)
    if two_sided:
        s = np.abs(r)
    elif direction == "over":
        s = r
    elif direction == "under":
        s = -r
    else:
        raise ValueError(f"direction must be one of {_DIRECTIONS}, got {direction!r}")
    return _group_max(s, groups)


def _quantile(scores: np.ndarray, alpha: float) -> float:
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    n = int(scores.size)
    k = conformal_rank(n, alpha)
    if n == 0 or k > n:
        warnings.warn(
            f"conformal margin undefined: n_groups={n} is too few for alpha={alpha} "
            f"(needs ceil((n+1)(1-alpha))={k} <= n); returning inf",
            RuntimeWarning, stacklevel=3)
        return float("inf")
    return float(np.sort(scores)[k - 1])


def one_sided_margin(residuals, groups, alpha: float = 0.10, direction: str = "over") -> float:
    """Conformal margin on the unsafe side only (may be negative: no group was ever unsafe at that level)."""
    return _quantile(group_scores(residuals, groups, direction), alpha)


def two_sided_margin(residuals, groups, alpha: float = 0.10) -> float:
    """Symmetric conformal margin on |residual| (use log residuals for density-like/viscosity quantities)."""
    return _quantile(group_scores(residuals, groups, two_sided=True), alpha)


def coverage(residuals, groups, margin: float, direction: str = "over", two_sided: bool = False) -> tuple[int, int]:
    """(groups whose max unsafe-side residual <= margin, n_groups).  Counts only; used for the locked check."""
    s = group_scores(residuals, groups, direction, two_sided)
    if s.size == 0:
        return 0, 0
    return int(np.sum(s <= margin)), int(s.size)
