from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "C7" / "d86_eval.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C7.md"


def test_d86_eval_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_d86_eval_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 6
    required_cols = [
        "cut", "description", "model", "cv_scheme", "n_fuels",
        "mae_C", "mae_ci95", "rmse_C", "rmse_ci95", "bias_C"
    ]
    for col in required_cols:
        assert col in df.columns

    cuts = set(df["cut"])
    assert "T10_C" in cuts
    assert "T50_C" in cuts
    assert "T90_C" in cuts

    # Verify that LOO recalibrated model achieves MAE < 5.0 C on all cuts
    recal = df[df["cv_scheme"] == "leave_one_fuel_out_cv"].set_index("cut")
    assert recal.loc["T10_C", "mae_C"] < 5.0
    assert recal.loc["T50_C", "mae_C"] < 5.0
    assert recal.loc["T90_C", "mae_C"] < 5.0

    # Verify that T50 zero-parameter physics has MAE ~ 4.05 C and near-zero bias
    t50_phys = df[(df["cut"] == "T50_C") & (df["cv_scheme"] == "zero_parameters")].iloc[0]
    assert abs(t50_phys["mae_C"] - 4.05) < 0.20
    assert abs(t50_phys["bias_C"]) < 0.50
