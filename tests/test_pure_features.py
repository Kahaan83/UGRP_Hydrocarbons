"""Pure-compound feature tests (WP4). Run: python tests/test_pure_features.py"""
import sys, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "src"))
from fuelblend.features import pure_features as pf

def _rows():
    return [
        dict(inchikey="TESTKEY0001", smiles="CCCC", T_K=300.0, cls="n_paraffin"),
        dict(inchikey="TESTKEY0002", smiles="CCCC", T_K=350.0, cls="aromatic"),
    ]

def test_feature_matrix_finite():
    X = pf.pure_feature_frame(_rows())
    assert X.shape == (2, pf.N_FEATURES)
    assert np.isfinite(X).all()
    assert pf.N_FEATURES <= 25

def test_logT_and_class_onehot():
    X = pf.pure_feature_frame(_rows())
    assert abs(X[0, pf.FEATURE_NAMES.index("logT")] - math.log(300.0)) < 1e-9
    assert abs(X[0, pf.FEATURE_NAMES.index("T_K")] - 300.0) < 1e-9
    assert abs(X[0, pf.FEATURE_NAMES.index("inv_T_K")] - 1000.0 / 300.0) < 1e-6
    i_np = pf.FEATURE_NAMES.index("cls_n_paraffin")
    i_ar = pf.FEATURE_NAMES.index("cls_aromatic")
    assert X[0, i_np] == 1.0 and X[0, i_ar] == 0.0
    assert X[1, i_np] == 0.0 and X[1, i_ar] == 1.0

if __name__ == "__main__":
    test_feature_matrix_finite()
    test_logT_and_class_onehot()
    print("ALL PURE FEATURE TESTS PASSED")
