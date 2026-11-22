import numpy as np
from fuelblend.uq import conformal_v7 as C


def test_margin_exact_on_uniform():
    r = np.arange(1, 101, dtype=float)
    assert C.one_sided_margin(r, 0.10) == 91.0            # rank ceil(101 * 0.9) = 91
    # unit weights reproduce the unweighted margin
    assert C.one_sided_margin(r, 0.10, weights=np.ones(100)) == 91.0


def test_mondrian_fallback_small_group():
    rng = np.random.default_rng(0)
    r = np.r_[rng.normal(0, 1, 60), rng.normal(5, 1, 5)]
    g = np.array(["big"] * 60 + ["small"] * 5)
    m = C.group_conditional_margins(r, g, 0.10, min_n=10)
    assert m["small"] == m["__global__"]                   # fewer than 10 rows: global margin
    assert m["big"] == C.one_sided_margin(r[:60], 0.10)


def test_weights_clipped():
    rng = np.random.default_rng(1)
    Xc = rng.normal(0, 1, (200, 2))
    Xt = rng.normal(4, 0.3, (30, 2))                       # strong shift
    w = C.shift_weights(Xc, Xt)
    assert w.min() >= C.W_LO - 1e-12 and w.max() <= C.W_HI + 1e-12
    assert w.argmax() == np.argmax(Xc @ np.array([1.0, 1.0]))   # largest weight where calibration looks most like the test cloud
