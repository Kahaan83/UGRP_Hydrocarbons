"""Unit tests for Phase 2 Octane Deliverables (T2.16 - T2.19)."""
from pathlib import Path
import pandas as pd

from fuelblend.data.validate import validate_measurement_df
from fuelblend.baselines.octane import blend_octane_linear_volume, blend_octane_interaction, calculate_octane_sensitivity

ROOT = Path(__file__).resolve().parents[1]
OCTANE_MEAS_CSV = ROOT / "data" / "processed" / "measurements_octane.csv"
CV_CSV = ROOT / "results" / "t2" / "octane_baseline_grouped_cv.csv"
SENS_CSV = ROOT / "results" / "t2" / "octane_sensitivity_report.csv"
NOTES_MD = ROOT / "docs" / "octane_data_notes.md"


def test_octane_measurements_labels_and_validation():
    assert OCTANE_MEAS_CSV.exists(), f"Missing {OCTANE_MEAS_CSV}"
    df = pd.read_csv(OCTANE_MEAS_CSV)
    assert len(df) >= 1000

    # T2.16 acceptance criterion: every row labeled as pure, blending_value, or measured_mixture
    assert "octane_data_type" in df.columns
    allowed_labels = {"pure", "blending_value", "measured_mixture"}
    assert set(df["octane_data_type"].unique()) == allowed_labels

    # Validate against strict measurement validator
    cols_to_val = [c for c in df.columns if c != "octane_data_type"]
    val_res = validate_measurement_df(df[cols_to_val], strict=True)
    assert val_res["valid"] is True


def test_octane_grouped_cv_results():
    assert CV_CSV.exists(), f"Missing {CV_CSV}"
    df = pd.read_csv(CV_CSV)
    assert len(df) == 2

    lin = df[df["model"].str.startswith("Linear")].iloc[0]
    bon = df[df["model"].str.startswith("Interaction")].iloc[0]

    # Interaction model improves upon linear baseline
    assert bon["mae"] < lin["mae"]
    assert bon["rmse"] < lin["rmse"]


def test_octane_sensitivity_report():
    assert SENS_CSV.exists(), f"Missing {SENS_CSV}"
    df = pd.read_csv(SENS_CSV)
    assert len(df) == 1
    assert df["mean_octane_sensitivity"].iloc[0] > 0.0


def test_octane_data_notes_exists():
    assert NOTES_MD.exists(), f"Missing {NOTES_MD}"
    text = NOTES_MD.read_text(encoding="utf-8")
    assert "CFR" in text
    assert "reproducibility" in text.lower()
