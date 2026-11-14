"""AJFTD smoke point + naphthalenes loader (v5 delta_smoke_cert 1.1). Fixture = real rows of ajftd_all_metrics.csv."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.data.loaders_ajftd_v5 import load_ajftd_smoke_naph

FIX = Path(__file__).parent / "data" / "ajftd_smoke_naph_sample.csv"
COLS = ["fuel_id", "source_group", "smoke_point_mm", "smoke_qualifier", "naphthalenes_vol_pct",
        "aromatics_vol_pct", "method_smoke"]


@pytest.fixture()
def df(tmp_path):
    arom = tmp_path / "long.csv"
    pd.DataFrame({"fuel_id": [156, 156, 10], "property": ["aromatics_vol_pct"] * 3,
                  "value": [18.0, 20.0, 12.0]}).to_csv(arom, index=False)
    return load_ajftd_smoke_naph(raw_path=FIX, aromatics_path=arom).set_index("fuel_id")


def test_columns(df):
    assert list(df.reset_index().columns) == COLS


def test_no_test_dropped(df):
    assert 1003 not in df.index and 1015 not in df.index  # only 'No Test' / 'Not Found' smoke rows


def test_censored_smoke_kept_as_qualifier(df):
    assert np.isnan(df.loc[165, "smoke_point_mm"]) and df.loc[165, "smoke_qualifier"] == ">"
    # a fuel with censored AND uncensored rows keeps the median of the uncensored ones
    assert df.loc[161, "smoke_point_mm"] == pytest.approx(42.0)
    assert df.loc[161, "smoke_qualifier"] == ""


def test_median_of_repeats(df):
    assert df.loc[156, "smoke_point_mm"] == pytest.approx(22.5)
    assert df.loc[156, "method_smoke"] == "D1322"
    assert df.loc[156, "aromatics_vol_pct"] == pytest.approx(19.0)


def test_naphthalenes_vol_only(df):
    # fuel 24444 has mass% 1.6 (D2425) and vol% 1.1 (D1840); only vol% is used, no conversion
    assert 24444 in df.index
    assert df.loc[24444, "naphthalenes_vol_pct"] == pytest.approx(1.1)
