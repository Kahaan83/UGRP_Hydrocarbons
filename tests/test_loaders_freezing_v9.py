"""Tests for v9 real-fuel freezing loader."""
from __future__ import annotations

import json
import pytest
from fuelblend.data import loaders_freezing_v9 as L


def test_row_counts():
    props, comps = L.load_freezing_v9()
    assert len(props) == 12
    assert (props["property"] == "freezing_point_C").sum() == 4
    assert (props["property"] == "flash_point_C").sum() == 4
    assert (props["property"] == "smoke_point_mm").sum() == 4
    assert len(comps) == 268
    assert comps["fuel_label"].nunique() == 4


def test_spot_values():
    props, comps = L.load_freezing_v9()
    p = props.set_index(["sample_id", "property"]).value
    
    # Table A-3, p14: POSF-10264 (A-1 JP-8)
    assert p[("POSF-10264", "freezing_point_C")] == -52.0
    assert p[("POSF-10264", "flash_point_C")] == 42.0
    
    # Table A-3, p14: POSF-10289 (A-3 JP-5)
    assert p[("POSF-10289", "freezing_point_C")] == -50.0
    assert p[("POSF-10289", "flash_point_C")] == 65.0
    
    # Table A-3, p14: POSF-10325 (A-2 Jet A)
    assert p[("POSF-10325", "freezing_point_C")] == -51.0
    
    # Table A-3, p14: POSF-4658
    assert p[("POSF-4658", "freezing_point_C")] == -48.0


def test_composition_closure_and_carbon_split():
    props, comps = L.load_freezing_v9()
    # Composition sums to 100 within 0.5 wt% per fuel
    for flabel, grp in comps.groupby("fuel_label"):
        assert abs(grp["mass_pct"].sum() - 100.0) < 0.5, f"{flabel} closure failed"
    # Carbon split is present (no_carbon_split is False)
    assert not comps["no_carbon_split"].any()
    assert (comps["n_carbon"] >= 7).all()


def test_duplicate_system_ids_resolved():
    props, comps = L.load_freezing_v9()
    p_sys = props.set_index("sample_id")["system_id"].to_dict()
    assert p_sys["POSF-10325"] == "fuel:ajftd:2"
    assert p_sys["POSF-10289"] == "fuel:ajftd:39"
    assert p_sys["POSF-4658"] == "fuel:ajftd:23"
    assert p_sys["POSF-10264"] == "fuel:afrl:posf10264"
