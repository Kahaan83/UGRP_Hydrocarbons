"""Tests for Task D4: Hierarchical Stacked Model Evaluation."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "D4" / "stacked_eval.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "D4.md"


def test_d4_deliverables_exist():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_d4_stacked_eval_contents():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 20, f"Expected 20 rows (5 props x 4 models), got {len(df)}"

    expected_props = {
        "viscosity_kin_cSt_ln",
        "density_15C_kg_m3",
        "freezing_point_C",
        "flash_point_C",
        "heat_of_combustion_net_MJ_kg",
    }
    assert set(df["property"].unique()) == expected_props

    expected_models = {
        "Raw Physics",
        "Recalibrated Physics (LOO)",
        "Stacked Model (Physics + Recal + Shrunk Res)",
        "Pure Learned (Direct / Neural Extrap)",
    }
    assert set(df["model"].unique()) == expected_models


def test_d4_viscosity_recalibration_dominance():
    df = pd.read_csv(CSV_PATH)
    kv = df[df["property"] == "viscosity_kin_cSt_ln"].set_index("model")

    raw_mae = kv.loc["Raw Physics", "mae"]
    recal_mae = kv.loc["Recalibrated Physics (LOO)", "mae"]
    stacked_mae = kv.loc["Stacked Model (Physics + Recal + Shrunk Res)", "mae"]
    learned_mae = kv.loc["Pure Learned (Direct / Neural Extrap)", "mae"]

    # Recalibrated physics must outperform raw Refutas by >60%
    assert recal_mae < 0.15, f"Expected recal_mae < 0.15, got {recal_mae}"
    assert raw_mae > 0.30, f"Expected raw_mae > 0.30, got {raw_mae}"
    # Stacked model must preserve physics dominance (mae essentially equal to recalibrated physics)
    assert abs(stacked_mae - recal_mae) < 0.01
    # Pure learned (DeepSets extrapolation) must be substantially worse than recalibrated physics on real fuels
    assert learned_mae > recal_mae * 2.0


def test_d4_freezing_point_recalibration():
    df = pd.read_csv(CSV_PATH)
    fp = df[df["property"] == "freezing_point_C"].set_index("model")

    raw_mae = fp.loc["Raw Physics", "mae"]
    recal_mae = fp.loc["Recalibrated Physics (LOO)", "mae"]
    stacked_mae = fp.loc["Stacked Model (Physics + Recal + Shrunk Res)", "mae"]

    # Raw WAT baseline is ~14.6 C offset from freezing point
    assert raw_mae > 13.0
    # Recalibrated WAT - 14.58 C brings error down to ~4 C
    assert recal_mae < 5.0
    # Stacked model achieves <= recalibrated physics
    assert stacked_mae <= recal_mae + 0.1
