"""Tests for v10 real-fuel properties and composition loaders."""
from __future__ import annotations

from pathlib import Path
import pytest
from scripts.v9.build_unified_v9 import _label_table
from fuelblend.data.loaders_fuels_v10 import (
    load_fuel_properties_v10,
    load_fuel_composition_v10,
)

REPO = Path(__file__).resolve().parents[1]
FIXTURE_PROPS = REPO / "tests" / "data" / "fuels_v10_fixture.csv"
FIXTURE_COMPS = REPO / "tests" / "data" / "fuels_v10_composition_fixture.csv"


def test_minus_sign_variants():
    """−47, –47 and -47 all load as -47.0."""
    df = load_fuel_properties_v10(FIXTURE_PROPS)
    fz = df[df["property"] == "freezing_point_C"]
    assert len(fz) >= 3
    # All 3 fuels in fixture (POSF-10264, POSF-10289, POSF-4658) had minus variants
    vals = fz["value"].tolist()
    for v in vals[:3]:
        assert v == -47.0


def test_open_composition_flagged():
    """A fixture fuel summing to 97 % gets composition_open True."""
    df = load_fuel_composition_v10(FIXTURE_COMPS)
    open_fuel = df[df["sample_id"] == "POSF-OPEN"]
    assert len(open_fuel) > 0
    assert (open_fuel["composition_open"] == True).all()

    closed_fuel = df[df["sample_id"] == "POSF-CLOSED"]
    assert len(closed_fuel) > 0
    assert (closed_fuel["composition_open"] == False).all()


def test_class_only_flagged():
    """A hydrocarbon-type table with no carbon numbers gets class_only True."""
    df = load_fuel_composition_v10(FIXTURE_COMPS)
    class_fuel = df[df["sample_id"] == "POSF-CLASS"]
    assert len(class_fuel) > 0
    assert (class_fuel["class_only"] == True).all()

    closed_fuel = df[df["sample_id"] == "POSF-CLOSED"]
    assert (closed_fuel["class_only"] == False).all()


def test_posf_maps_to_ajftd():
    """POSF-10325 maps through canonical_id to the AJFTD id (read from label table, not hard-coded)."""
    labels = _label_table()
    expected_id = labels["10325"]

    df = load_fuel_properties_v10(FIXTURE_PROPS)
    posf_rows = df[df["sample_id"] == "POSF-10325"]
    assert len(posf_rows) > 0
    assert (posf_rows["system_id"] == expected_id).all()


def test_units():
    """cSt kept, mm2/s maps to cSt, degF converted."""
    df = load_fuel_properties_v10(FIXTURE_PROPS)

    # mm2/s -> cSt
    mm2_row = df[(df["sample_id"] == "POSF-10325") & (df["T_C"] == -20.0)].iloc[0]
    assert mm2_row["unit"] == "cSt"
    assert mm2_row["value"] == 4.5

    # cSt kept
    cst_row = df[(df["sample_id"] == "POSF-10325") & (df["T_C"] == -40.0)].iloc[0]
    assert cst_row["unit"] == "cSt"
    assert cst_row["value"] == 8.5

    # degF converted: 104 degF = (104 - 32) * 5/9 = 40.0 C
    degf_row = df[df["sample_id"] == "POSF-11498"].iloc[0]
    assert degf_row["unit"] == "C"
    assert pytest.approx(degf_row["value"], 0.01) == 40.0
