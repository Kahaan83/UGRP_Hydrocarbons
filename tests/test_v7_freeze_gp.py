import numpy as np
import pandas as pd
from fuelblend.v7 import freeze_gp as F


def _comp():
    return pd.DataFrame(dict(
        fuel=["a"] * 4 + ["b"] * 2,
        smiles=["CCCCCCCCCC", "CCCCCCCCCCCCCCCC", "c1ccc2ccccc2c1", "C1CCC2CCCCC2C1", "CCCCCCCC", "Cc1ccccc1"],
        frac=[0.4, 0.3, 0.2, 0.1, 0.5, 0.5], tm_K=[243.5, 291.0, 353.4, 230.0, 216.4, np.nan]))


def test_heavy_tail_features_sum_bounds():
    f = F.heavy_tail_features(_comp())
    assert list(f.columns) == F.FEATS
    m = f[["nC14plus_mass", "naphthalene_mass", "dicyclo_mass"]]
    assert ((m >= 0) & (m <= 1)).all().all()
    assert (f.naphthalene_mass + f.dicyclo_mass <= 1 + 1e-12).all()      # aromatic vs non-aromatic: disjoint
    assert (f.tm_p90_mass <= f.tm_max_present + 1e-12).all()
    assert abs(f.loc["a", "nC14plus_mass"] - 0.3) < 1e-9 and abs(f.loc["a", "naphthalene_mass"] - 0.2) < 1e-9
    assert abs(f.loc["a", "dicyclo_mass"] - 0.1) < 1e-9 and f.loc["a", "tm_max_present"] == 353.4


def test_novelty_flags_far_point():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 3))
    thr = F.novelty_threshold(X)
    near, far = F.novelty_score(X, np.zeros((1, 3)))[0], F.novelty_score(X, np.full((1, 3), 10.0))[0]
    assert near < thr < far
