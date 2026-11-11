"""v5 D1: AJFTD bulk-descriptor loader.  Fixture = ~44 rows sliced from the real scrape (tests/data/)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.data import loaders_ajftd as L

FIX = Path(__file__).parent / "data" / "ajftd_all_metrics_sample.csv"


@pytest.fixture(scope="module")
def sample():
    raw = pd.read_csv(FIX, encoding="utf-8-sig", low_memory=False)
    log = {}
    return raw, L.build_long(raw, log), log


def test_numeric_only_and_provenance(sample):
    _, lg, _ = sample
    assert list(lg.columns) == L.LONG_COLS
    assert lg["value"].notna().all() and np.isfinite(lg["value"]).all()
    assert (lg["provenance_class"] == "real_fuel_bulk_descriptors").all()
    assert not lg.duplicated().any()


def test_duplicate_and_nonvalues_dropped(sample):
    raw, lg, log = sample
    assert log["dropped"]["exact_duplicate_rows"] >= 1
    assert "sulfur" not in set(lg["property"])
    assert len(lg) < len(raw)


def test_viscosity_temperature_encoding(sample):
    _, lg, _ = sample
    v = lg[lg["property"] == "viscosity_kin_cSt"]
    assert set(v["T_source"]) <= {"reported", "assumed_spec_-20C", "reported_parsed_text"}
    rep = v[v["T_source"] == "reported"]
    assert set(np.round(rep["T_K"], 2)) <= {253.15, 313.15, 233.15, 293.15}
    ass = v[v["T_source"] == "assumed_spec_-20C"]
    assert len(ass) >= 1 and np.allclose(ass["T_K"], 253.15)


def test_density_api_and_mislabelled_heat(sample):
    _, lg, _ = sample
    d = lg[(lg["property"] == "density_kg_m3") & (lg["derived"] == "from API gravity")]
    assert len(d) and d["value"].between(780, 880).all()
    assert abs(L.api_to_density_kg_m3(44.0) - 141.5 / 175.5 * 999.016) < 1e-9
    h = lg[lg["property"] == "heat_of_combustion_net_MJ_kg"]
    assert h["value"].between(40, 46).all()


def test_distillation_without_method_kept(sample):
    raw, lg, _ = sample
    nm = raw[(raw["metric"] == "distillation") & raw["Test Method"].isna() & pd.to_numeric(raw["Value"], errors="coerce").notna()]
    assert len(nm) and (lg["property"].str.startswith("distillation_T")).sum() >= len(nm)


def test_source_group_and_year():
    assert L.source_group("JF KENNEDY INT'L AIRPORT_1194") == "JFK"
    assert L.source_group("JAA20110073") == "JAA" and L.source_group("JA120130024") == "JA"
    assert L.source_group("HR820100004") == "HR" and L.source_group("(02239-070416)T108") == "PAREN_BATCH"
    assert L.source_group("10245") == "NUMERIC_ID" and L.source_group("1600776-002(37176) x") == "LAB_REPORT"
    assert L._test_year("7/5/2016 00:00") == 2016.0


def test_qualifier_and_temp_parsing():
    assert L.parse_value("<1") == (1.0, "<")
    assert np.isnan(L.parse_value("No Test")[0])
    assert L.parse_temp_c("-20°C (-4°F),,,-40") == -20.0
    assert np.isnan(L.parse_temp_c("No Test"))


def test_wide_table(sample):
    _, lg, _ = sample
    w = L.build_wide(lg)
    assert w["fuel_id"].is_unique
    assert "viscosity_kin_cSt_m20C" in w.columns and "density_kg_m3_15C" in w.columns


# ---- spot check of 10 values against the raw scrape (skipped if extras/ is absent) ----
SPOT = [  # (fuel_id, property, T_K or None, value present among raw replicates)
    ("2", "flash_point_C", None, 48.0), ("2", "flash_point_C", None, 49.0),
    ("2", "freezing_point_C", None, -52.0), ("2", "aromatics_vol_pct", None, 17.4),
    ("2", "density_kg_m3", 288.15, 804.0), ("2", "distillation_T100_C", None, 270.0),
    ("2", "distillation_T0_C", None, 159.0), ("2", "viscosity_kin_cSt", 253.15, 5.2),
    ("2", "viscosity_kin_cSt", 233.15, 9.318), ("2", "heat_of_combustion_net_MJ_kg", None, 43.06),
]


@pytest.mark.skipif(not L.RAW.exists(), reason="raw AJFTD scrape not on disk")
def test_spot_check_against_raw():
    raw = pd.read_csv(L.RAW, encoding="utf-8-sig", low_memory=False)
    lg = L.build_long(raw[raw["fuel_id"].astype(str) == "2"])
    for fid, prop, T, val in SPOT:
        s = lg[(lg["fuel_id"].astype(str) == fid) & (lg["property"] == prop)]
        if T is not None:
            s = s[np.isclose(s["T_K"], T, atol=0.01)]
        assert np.isclose(s["value"], val, atol=1e-6).any(), (prop, T, val, sorted(s["value"].unique()))
