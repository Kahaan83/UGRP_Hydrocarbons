from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "C1" / "hoc_eval.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C1.md"
GCV_V2_PATH = ROOT / "data" / "processed" / "v2" / "gcv_v2_MJ_kg.csv"


def test_hoc_v2_outputs_exist():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"
    assert GCV_V2_PATH.exists(), f"Missing {GCV_V2_PATH}"


def test_hoc_v2_evaluation_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 2
    required_cols = [
        "baseline_version", "tier_coverage", "n_fuels", "fuel_family",
        "mae_MJ_kg", "mae_ci95", "rmse_MJ_kg", "rmse_ci95", "mape_pct", "bias_MJ_kg", "error_reduction_pct"
    ]
    for col in required_cols:
        assert col in df.columns

    v1_row = df[df["baseline_version"].str.contains("v1")].iloc[0]
    v2_row = df[df["baseline_version"].str.contains("v2")].iloc[0]

    # Verify that v2 reduces MAE by > 60%
    assert v2_row["mae_MJ_kg"] < 0.25
    assert v1_row["mae_MJ_kg"] > 0.65
    assert v2_row["error_reduction_pct"] > 60.0

    # Verify that Dulong tier C was completely eliminated from gcv_v2
    gcv_v2 = pd.read_csv(GCV_V2_PATH)
    assert len(gcv_v2) == 1038
    assert (gcv_v2["tier"] == "C").sum() == 0
    assert (gcv_v2["tier"] == "B").sum() == 924
    assert (gcv_v2["tier"] == "A").sum() == 114
