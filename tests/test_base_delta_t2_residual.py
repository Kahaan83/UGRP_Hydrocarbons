"""v5 BASE Task 13: T2 residual model."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.mixture.t2_residual import FEATURE_NAMES, bin_features, fit_predict_residual  # noqa: E402


def test_constant_residual_recovered():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(30, 5))
    p = fit_predict_residual(X, np.full(30, 0.37), np.arange(30) % 10, rng.normal(size=(7, 5)))
    assert np.allclose(p, 0.37, atol=1e-6)


def test_test_rows_never_enter_scaler_or_alpha():
    rng = np.random.default_rng(1)
    X, r, g = rng.normal(size=(40, 4)), rng.normal(size=40), np.arange(40) % 8
    Xt = rng.normal(size=(5, 4))
    p1 = fit_predict_residual(X, r, g, Xt)
    Xt2 = Xt.copy()
    Xt2[1:] *= 1e6                                       # other test rows changed wildly: row 0 prediction must not move
    p2 = fit_predict_residual(X, r, g, Xt2)
    assert np.isclose(p1[0], p2[0])


def test_feature_vector_shape():
    b = pd.DataFrame(dict(hc_class=["n_paraffin", "alkylbenzene", "no_rep:diaromatic"], n_carbon=[10, 9, 12], mass_pct=[50, 30, 20]))
    f = bin_features(b)
    assert f.shape == (len(FEATURE_NAMES),) and np.isclose(f[:10].sum(), 1.0)


def test_train_and_test_groups_disjoint():
    s = importlib.util.spec_from_file_location("dt2r", ROOT / "scripts/v5/base/delta_t2_residual.py")
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    X = m.features()
    d = m.targets()["flash_point_C"]
    for test_src, train in m.DIRECTIONS:
        te = m.run_direction(d, X, test_src, train)          # asserts disjointness internally
        assert te is not None and te.source.eq(test_src).all()
