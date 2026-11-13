import numpy as np
import pandas as pd

from fuelblend.thermo import fusion_estimate as FE


def test_family_of():
    cases = {"CCCCCCCCCC": "n_alkane", "CC(C)CC(C)(C)C": "iso_alkane", "C=CCCCC": "acyclic_unsaturated",
             "CCCCC1CCCCC1": "cycloalkane_mono", "C1CCC2CCCCC2C1": "cycloalkane_poly", "C1=CCCCC1": "cycloalkene",
             "Cc1ccccc1": "benzene_mono", "c1ccc2c(c1)CCCC2": "aromatic_naphthenic", "Cc1cccc2ccccc12": "aromatic_poly",
             "c1ccc(cc1)-c1ccccc1": "aromatic_poly"}
    for s, f in cases.items():
        assert FE.family_of(s) == f, s
    assert set(cases.values()) == set(FE.FAMILIES)


def test_features_finite_core_columns():
    X = FE.features(["CCCCCCCCCC", "Cc1ccccc1", "C1CCC2CCCCC2C1"])
    assert len(X) == 3
    assert np.isfinite(X[["n_carbon", "mw", "wiener", "symmetry"]].to_numpy()).all()
    assert X.n_carbon.tolist() == [10, 7, 10]


def test_measured_never_replaced():
    measured = {"A": {"tm_K": 250.0, "source": "verified"}}
    est = {"A": {"tm_K": 999.0, "source": FE.SOURCE}, "B": {"tm_K": 200.0, "source": FE.SOURCE}}
    out = FE.fill(measured, est)
    assert out["A"]["tm_K"] == 250.0 and out["A"]["source"] == "verified"
    assert out["B"]["source"] == FE.SOURCE


def test_oof_is_leave_group_out():
    X = pd.DataFrame({"a": np.arange(20, dtype=float)})
    y = np.arange(20, dtype=float)
    g = np.repeat(["x", "y"], 10)
    p = FE.oof_predict(X, y, g)
    assert np.isfinite(p).all()
    # a model trained only on group y (a >= 10) cannot reach the low targets of group x
    assert p[:10].min() >= 9.0
