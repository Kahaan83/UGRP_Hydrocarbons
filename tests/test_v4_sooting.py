"""Task 14 tests: YSI scale anchors, linear mole mixing, D1322 censoring."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines import sooting as S  # noqa: E402

pytestmark = pytest.mark.skipif(not (S.DATA / "ysi_hc.csv").exists(), reason="sooting inputs not built")


def test_ysi_unified_scale_anchor():
    t = S.ysi_table()
    assert t[S._canon("c1ccccc1")] == pytest.approx(100.3, abs=0.05)
    assert t[S._canon("CCCCCCC")] == pytest.approx(36.0, abs=0.05)
    assert t[S._canon("Cc1ccccc1")] == pytest.approx(170.9, abs=0.05)


def test_ysi_mix_linear():
    y = S.ysi_mix(["c1ccccc1", "CCCCCCC"], [0.5, 0.5])
    assert y == pytest.approx((100.3 + 36.0) / 2, abs=0.05)
    # unnormalised input is renormalised; unknown component is an error, not a silent fill
    assert S.ysi_mix(["c1ccccc1", "CCCCCCC"], [1, 1]) == pytest.approx(y)
    with pytest.raises(KeyError):
        S.ysi_mix(["C"], [1.0])


def test_smoke_point_censoring_respected():
    df = pd.read_csv(S.DATA / "smoke_point_pure_citedtyped.csv")
    assert df.censored.sum() >= 3
    ex = S.exact_rows(df)
    assert not ex.censored.any() and (ex.smoke_point_mm < S.SP_CAP_MM).all()
    c_with = S.fit_ysi_sp(df)
    c_without = S.fit_ysi_sp(df[~df.censored])
    assert c_with == pytest.approx(c_without)
    # a censored row with an absurd value must not move the fit
    bad = pd.concat([df, df.iloc[[0]].assign(smoke_point_mm=500.0, censored=True)])
    assert S.fit_ysi_sp(bad) == pytest.approx(c_with)
    # predictions never exceed the lamp cap, bulk fit ignores censored rows too
    assert S.smoke_point_from_ysi(1.0, 100.0) <= S.SP_CAP_MM
    a = pd.read_csv(S.DATA / "smoke_point_ajftd.csv")
    b = pd.concat([a, a.iloc[[0]].assign(value=50.0, censored=True)])
    assert np.allclose(S.fit_bulk(a)[0], S.fit_bulk(b)[0])
    with pytest.raises(ValueError):
        S.fit_ysi_sp(df[df.censored])


def test_bulk_model_direction():
    assert S.smoke_point_from_bulk(10.0, 1.0) > S.smoke_point_from_bulk(22.0, 1.0)
