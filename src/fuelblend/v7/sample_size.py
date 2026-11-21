"""Group-level subsampling for sample-size learning curves; nested across fractions for a fixed seed."""
from __future__ import annotations
import numpy as np

FRACS = (0.1, 0.2, 0.4, 0.6, 0.8, 1.0)


def group_subsample(groups: np.ndarray, frac: float, rng: np.random.Generator) -> np.ndarray:
    """Indices of rows whose group is in a random subset of round(frac * n_groups) groups (>=1).

    The group order is one permutation drawn from rng, and the first k groups are kept, so for the same rng seed the
    groups at a smaller fraction are a subset of those at a larger one."""
    g = np.asarray(groups)
    ug = np.unique(g)
    order = rng.permutation(len(ug))
    k = max(1, int(round(frac * len(ug))))
    keep = set(ug[order[:k]].tolist())
    return np.where([x in keep for x in g.tolist()])[0]


def slope_ci(n_train, mae_draws, n_boot: int = 2000, seed: int = 0) -> tuple[float, float, float]:
    """Slope of MAE vs log(n) over given n values. mae_draws: (n_units, len(n_train)) per-unit (repeat x fold) MAEs.

    Returns (lo, hi, point) with a 90 % bootstrap CI over units."""
    x = np.log(np.asarray(n_train, float)); d = np.asarray(mae_draws, float)
    rng = np.random.default_rng(seed)
    slope = lambda m: np.polyfit(x, m.mean(0), 1)[0]
    bs = [slope(d[rng.integers(0, len(d), len(d))]) for _ in range(n_boot)]
    return float(np.percentile(bs, 5)), float(np.percentile(bs, 95)), float(slope(d))
