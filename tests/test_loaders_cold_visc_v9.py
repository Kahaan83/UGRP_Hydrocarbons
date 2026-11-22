"""Tests for v9 real-fuel cold viscosity loader."""
from __future__ import annotations

import pytest
from fuelblend.data import loaders_cold_visc_v9 as L


def test_row_counts():
    df = L.load_cold_visc_v9()
    assert len(df) > 50
    assert "viscosity_kin_cSt" in df["property"].values
    assert (df["property"] == "viscosity_kin_cSt").sum() >= 50


def test_spot_values():
    df = L.load_cold_visc_v9()
    # Edwards 2017 Table A-3 spot checks:
    ed = df[df["source_id"] == "edwards2017"].set_index(["sample_id", "property", "T_K"])
    assert ed.loc[("POSF-10264", "viscosity_kin_cSt", 253.15), "value"] == 3.5
    assert ed.loc[("POSF-10264", "viscosity_kin_cSt", 233.15), "value"] == 6.6
    assert ed.loc[("POSF-10289", "viscosity_kin_cSt", 253.15), "value"] == 6.5
    assert ed.loc[("POSF-10289", "viscosity_kin_cSt", 233.15), "value"] == 14.1
    assert ed.loc[("POSF-10325", "viscosity_kin_cSt", 253.15), "value"] == 4.5
    assert ed.loc[("POSF-10325", "viscosity_kin_cSt", 233.15), "value"] == 9.2

    # Boehm 2022 spot checks (neat fuels):
    bm = df[df["source_id"] == "boehm2022"].set_index(["sample_id", "property", "T_K"])
    v_a1_20 = bm.loc[("A-1 (POSF 10264)", "viscosity_kin_cSt", 253.15), "value"]
    val_a1_20 = v_a1_20.iloc[0] if hasattr(v_a1_20, "iloc") else v_a1_20
    assert abs(val_a1_20 - 3.4392) < 1e-3

    v_a1_40 = bm.loc[("A-1 (POSF 10264)", "viscosity_kin_cSt", 233.15), "value"]
    val_a1_40 = v_a1_40.iloc[0] if hasattr(v_a1_40, "iloc") else v_a1_40
    assert abs(val_a1_40 - 6.5151) < 1e-3

    v_a2_20 = bm.loc[("A-2 (POSF 10325)", "viscosity_kin_cSt", 253.15), "value"]
    val_a2_20 = v_a2_20.iloc[0] if hasattr(v_a2_20, "iloc") else v_a2_20
    assert abs(val_a2_20 - 4.4669) < 1e-3

    v_a2_40 = bm.loc[("A-2 (POSF 10325)", "viscosity_kin_cSt", 233.15), "value"]
    val_a2_40 = v_a2_40.iloc[0] if hasattr(v_a2_40, "iloc") else v_a2_40
    assert abs(val_a2_40 - 9.1359) < 1e-3


def test_temperature_and_canonical_units():
    df = L.load_cold_visc_v9()
    v = df[df["property"] == "viscosity_kin_cSt"]
    # All cold viscosity measurements must be stored in cSt
    assert (v["unit"] == "cSt").all()
    # Temperatures must be near 253.15 K (-20 C) or 233.15 K (-40 C)
    assert set(v["T_K"].round(2).unique()).issubset({253.15, 233.15})


def test_dynamic_density_rule():
    # If a value were in mPa.s without measured density, qualifier must be dynamic_no_density
    df = L.load_cold_visc_v9()
    dyn = df[df["property"] == "viscosity_dyn_mPas"]
    for _, r in dyn.iterrows():
        assert r["qualifier"] != "dynamic_no_density" or pd.isna(r["density_kg_m3"])
