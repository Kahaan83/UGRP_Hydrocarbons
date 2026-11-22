import numpy as np
from fuelblend.v7.sample_size import group_subsample, slope_ci


def test_subsample_is_group_level():
    g = np.repeat(np.arange(20), 5)
    idx = group_subsample(g, 0.4, np.random.default_rng(0))
    kept = np.unique(g[idx])
    assert len(kept) == 8
    assert all((g[idx] == k).sum() == 5 for k in kept)   # all rows of a kept group, no partial groups


def test_fracs_monotone_nesting():
    g = np.repeat(np.arange(30), 3)
    a = set(g[group_subsample(g, 0.2, np.random.default_rng(1))]); b = set(g[group_subsample(g, 0.4, np.random.default_rng(1))])
    c = set(g[group_subsample(g, 0.8, np.random.default_rng(1))])
    assert a <= b <= c


def test_slope_ci_flat_vs_falling():
    rng = np.random.default_rng(0)
    n = np.array([50, 100, 200, 400]); fl = slope_ci(n[-2:], rng.normal(5, .01, (30, 2)))
    assert fl[0] <= 0 <= fl[1] or abs(fl[2]) < .1
    dn = slope_ci(n[-2:], np.tile([4.0, 3.0], (30, 1)) + rng.normal(0, .01, (30, 2)))
    assert dn[1] < 0
