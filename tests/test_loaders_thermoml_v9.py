"""Tests for fuelblend.data.loaders_thermoml_v9 and ThermoML parser constraints.

Tests:
- Real post-2020 fixture (tests/data/thermoml_v9_tca_2020_178536.xml) parses to expected T, phase, and value
- je8006138.xml regression test stays green (constraint T/P under ConstraintType preserved)
- No excess properties extracted (excess viscosity, excess molar volume)
- No vapour-phase composition variables accepted
- P <= 1000 kPa filter enforced
- loader outputs conform to measurements_dedup schema with canonical units
"""
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import pytest

from fuelblend.data import thermoml_extract as TE
from fuelblend.data.loaders_boehm2022 import COLS

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_2020 = ROOT / "tests" / "data" / "thermoml_v9_tca_2020_178536.xml"
FIXTURE_JE8006138 = ROOT / "tests" / "data" / "je8006138.xml"
EXTRACTED_CSV = ROOT / "data" / "raw" / "thermoml_v9" / "extracted.csv"


def test_real_new_record_fixture_parses_expected():
    """Real post-2020 fixture parses expected T, phase, and density/viscosity values."""
    assert FIXTURE_2020.exists(), f"Fixture missing: {FIXTURE_2020}"
    root = ET.parse(FIXTURE_2020).getroot()
    stats = Counter()
    comps, rows = TE.parse_file(root, "10.1016/j.tca.2020.178536.xml", TE.WANTED, max_p=1000.0, stats=stats)
    
    # Check endo-dicyclopentadiene pure density at 288.15 K
    pure_dens = [r for r in rows if r["prop"] == "Mass density, kg/m3" and r["comp_type"] == "pure" and abs(r["T_K"] - 288.15) < 0.01]
    endo_dens = [r for r in pure_dens if comps[r["ids"][0]]["name"] == "endo-dicyclopentadiene"]
    assert len(endo_dens) == 1
    d_row = endo_dens[0]
    assert d_row["prop_phase"] == "Liquid"
    assert pytest.approx(d_row["value"], rel=1e-4) == 985.8
    assert pytest.approx(d_row["P_kPa"], rel=1e-3) == 101.2

    # Check viscosity (Pa*s) exists
    visc_rows = [r for r in rows if r["prop"] == "Viscosity, Pa*s"]
    assert len(visc_rows) > 0
    for vr in visc_rows:
        assert vr["prop_phase"] == "Liquid"
        assert vr["value"] > 0


def test_je8006138_regression_stays_green():
    """Historic trap: constraint T and P live under Constraint > ConstraintID > ConstraintType."""
    assert FIXTURE_JE8006138.exists(), f"Fixture missing: {FIXTURE_JE8006138}"
    root = ET.parse(FIXTURE_JE8006138).getroot()
    stats = Counter()
    comps, rows = TE.parse_file(root, "je8006138.xml", TE.WANTED, max_p=1000.0, stats=stats)
    assert len(rows) > 0
    # Every row must have non-NaN T_K and P_kPa
    for r in rows:
        assert r["T_K"] is not None and not np.isnan(r["T_K"]), f"T_K lost in {r}"
        assert r["P_kPa"] is not None and not np.isnan(r["P_kPa"]), f"P_kPa lost in {r}"


def test_no_excess_properties_extracted():
    """ThermoML trap: whitelist property names, drop excess properties (e.g. excess viscosity)."""
    assert "Excess molar volume" not in TE.WANTED
    assert "Excess viscosity" not in TE.WANTED
    for prop in TE.WANTED:
        assert not prop.lower().startswith("excess")


def test_no_vapour_phase_composition():
    """ThermoML trap: ignore vapour-phase composition variables (only Liquid phase composition)."""
    assert FIXTURE_2020.exists()
    root = ET.parse(FIXTURE_2020).getroot()
    stats = Counter()
    comps, rows = TE.parse_file(root, "fixture.xml", TE.WANTED, max_p=1000.0, stats=stats)
    for r in rows:
        assert r["prop_phase"] == "Liquid"


def test_pressure_filter_drops_above_1000kpa():
    """ThermoML trap: drop P > 1000 kPa."""
    assert FIXTURE_2020.exists()
    root = ET.parse(FIXTURE_2020).getroot()
    stats = Counter()
    comps, rows = TE.parse_file(root, "fixture.xml", TE.WANTED, max_p=1000.0, stats=stats)
    for r in rows:
        if r["P_kPa"] is not None and not np.isnan(r["P_kPa"]):
            assert r["P_kPa"] <= 1000.0


def test_loaders_thermoml_v9_output_schema():
    """Test loader produces valid measurements_dedup schema with canonical units."""
    from fuelblend.data.loaders_thermoml_v9 import load_thermoml_v9
    
    assert EXTRACTED_CSV.exists(), f"Extracted CSV missing: {EXTRACTED_CSV}"
    df = load_thermoml_v9(EXTRACTED_CSV)
    
    # 1. Output columns match COLS
    assert list(df.columns) == list(COLS)
    assert len(df) > 0

    # 2. Canonical units and ranges
    assert set(df["property"].unique()).issubset({"density_kg_m3", "viscosity_dyn_mPas", "vapor_pressure_kPa", "viscosity_kin_cSt"})
    
    # Density in kg/m3 (e.g., 600 - 1200)
    rho = df[df["property"] == "density_kg_m3"]
    assert (rho["value"] >= 500).all() and (rho["value"] <= 1500).all()

    # Dynamic viscosity in mPa*s (not Pa*s; values > 0.05 mPa*s for liquids)
    mu = df[df["property"] == "viscosity_dyn_mPas"]
    assert (mu["value"] >= 0.01).all() and (mu["value"] <= 1000).all()

    # 3. System IDs valid and hydrocarbon-only
    assert df["system_id"].notna().all()
    assert (df["is_hc"] == True).all()

    # 4. Spot check: endo-dicyclopentadiene density at 298.15 K
    spot = rho[(rho["ref"].str.contains("10.1016/j.tca.2020.178536")) & (rho["kind"] == "pure") & (abs(rho["T_K"] - 298.15) < 0.1)]
    assert len(spot) >= 1
    assert pytest.approx(spot.iloc[0]["value"], rel=1e-3) == 977.2
