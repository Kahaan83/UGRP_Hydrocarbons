from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "C4" / "cold_flow_v2.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C4.md"


def test_cold_flow_v2_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_cold_flow_v2_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 4
    required_cols = [
        "property", "fuel_family", "model", "cv_scheme", "n_fuels",
        "mae_C", "mae_ci95", "rmse_C", "rmse_ci95", "bias_C", "pearson_r"
    ]
    for col in required_cols:
        assert col in df.columns

    df_models = df.set_index("model")

    # Verify that raw SLE reproduces v1 headline numbers
    raw_row = df_models.loc["v1 Raw SLE WAT Baseline (Schroeder-van Laar)"]
    assert abs(raw_row["mae_C"] - 14.58) < 0.10
    assert abs(raw_row["rmse_C"] - 15.79) < 0.10
    assert abs(raw_row["pearson_r"] - 0.443) < 0.01

    # Verify that LOO offset cuts MAE by > 65% down to ~4.05 C
    offset_row = df_models.loc["v2 Transparent Offset SLE Baseline (LOO CV)"]
    assert offset_row["mae_C"] < 5.0
    assert abs(offset_row["bias_C"]) < 0.10

    # Verify diesel cloud point is evaluated separately
    diesel_row = df_models.loc["v2 Heavy-Paraffin SLE Wax Appearance (Diesel Cloud Point)"]
    assert diesel_row["fuel_family"] == "diesel"
    assert diesel_row["property"] == "cloud_point_C"
    assert diesel_row["mae_C"] < 2.0
