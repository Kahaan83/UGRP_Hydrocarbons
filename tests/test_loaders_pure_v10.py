"""Tests for v10 pure physical constants loader."""
from __future__ import annotations

from pathlib import Path
import pytest
from fuelblend.data.loaders_pure_v10 import (
    load_vapour_pressure_v10,
    load_lfl_pure_v10,
    load_pure_constants_v10,
)

REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "tests" / "data" / "pure_v10_fixture.csv"
DRAFT_VP = REPO / "data" / "raw" / "pure_v10" / "vapour_pressure_draft.csv"


def test_vp_units():
    """mmHg, bar and Pa all map to kPa."""
    df = load_vapour_pressure_v10(FIXTURE)
    # 760 mmHg -> 101.325 kPa
    row_mmhg = df[df["sample_id"] == "test1"].iloc[0]
    assert row_mmhg["unit"] == "kPa"
    assert pytest.approx(row_mmhg["value"], 0.01) == 101.325

    # 1.0 bar -> 100.0 kPa
    row_bar = df[df["sample_id"] == "test2"].iloc[0]
    assert row_bar["unit"] == "kPa"
    assert pytest.approx(row_bar["value"], 0.01) == 100.0

    # 100000.0 Pa -> 100.0 kPa
    row_pa = df[df["sample_id"] == "test3"].iloc[0]
    assert row_pa["unit"] == "kPa"
    assert pytest.approx(row_pa["value"], 0.01) == 100.0

    # 101.325 kPa -> 101.325 kPa
    row_kpa = df[df["sample_id"] == "test4"].iloc[0]
    assert row_kpa["unit"] == "kPa"
    assert pytest.approx(row_kpa["value"], 0.01) == 101.325


def test_lfl_measured_only():
    """Rows whose method contains estimate, Jones or homologue are rejected."""
    df = load_lfl_pure_v10(FIXTURE)
    sample_ids = df["sample_id"].tolist()
    # test5 is measured ASTM E681 -> KEPT
    assert "test5" in sample_ids
    # test6 is Jones estimate -> REJECTED
    assert "test6" not in sample_ids
    # test7 is homologue extrapolation -> REJECTED
    assert "test7" not in sample_ids
    # test8 is estimated from boiling point -> REJECTED
    assert "test8" not in sample_ids


def test_scanned_rows_forced_unsure():
    """Rows with extraction=image and confidence=clear load with confidence forced to unsure."""
    df = load_pure_constants_v10(FIXTURE)
    row_image = df[df["sample_id"] == "test9"].iloc[0]
    assert row_image["extraction"] == "image"
    assert row_image["confidence"] == "unsure"

    row_text = df[df["sample_id"] == "test10"].iloc[0]
    assert row_text["extraction"] == "text_layer"
    assert row_text["confidence"] == "clear"


def test_perhydrophenanthrene_draft_loads():
    """Vapour pressure draft has 25 points for perhydrophenanthrene."""
    df = load_vapour_pressure_v10(DRAFT_VP)
    assert len(df) == 25
    assert (df["compound"] == "perhydrophenanthrene").all()
    assert (df["unit"] == "kPa").all()
    assert pytest.approx(df["value"].min(), 0.01) == 8.292
    assert pytest.approx(df["value"].max(), 0.01) == 99.533
