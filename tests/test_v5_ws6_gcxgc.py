"""WS6 helper checks: feature matrix closure, Holm, baselines use training rows only."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "v5" / "base"))
sys.path.insert(0, str(ROOT / "src"))
import delta_ws6_gcxgc as W  # noqa: E402


def test_matrix_rows_sum_to_one_and_map_classes():
    d = pd.DataFrame({"group": ["n-alkane", "iso-alkane", "alkylbenzene", "mystery"], "cn": [12, 25, 7, 10], "w": [50.0, 30.0, 15.0, 5.0]})
    v, tot, frac = W.matrix(d, "group", "cn", "w", W.AV_MAP)
    assert abs(v.sum() - 1.0) < 1e-12 and abs(tot - 100.0) < 1e-12 and abs(frac - 0.95) < 1e-12
    assert v[W.CLASSES.index("iso_paraffin") * len(W.CN_BINS) + W.CN_BINS.index(17)] > 0   # C25 -> >=17 bin
    assert v[W.CLASSES.index("alkylbenzene") * len(W.CN_BINS) + W.CN_BINS.index(8)] > 0    # C7 -> <=8 bin


def test_holm_orders_and_stops():
    r = W.holm({"a": 0.001, "b": 0.02, "c": 0.5})
    assert r == {"a": True, "b": True, "c": False}


def test_baselines_use_training_rows_only():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, len(W.FEATS)))
    ytr, ptr = rng.normal(size=12), rng.normal(size=12)
    out = W.fit_predict(X[:12], ytr, ptr, X[12:], np.zeros(8), False)
    assert np.allclose(out["C0"], ytr.mean()) and np.allclose(out["C2"], (ytr - ptr).mean())
