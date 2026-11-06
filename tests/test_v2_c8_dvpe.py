from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "C8" / "dvpe_eval.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C8.md"


def test_dvpe_eval_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_dvpe_eval_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 2
    required_cols = [
        "property", "model", "cv_scheme", "n_fuels",
        "mae_kpa", "mae_ci95", "rmse_kpa", "rmse_ci95", "bias_kpa"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify that real fuels were evaluated
    assert df["n_fuels"].iloc[0] == 13

    # Verify that LOO recalibration reduces bias
    raw_bias = df[df["cv_scheme"] == "zero_parameters"]["bias_kpa"].iloc[0]
    loo_bias = df[df["cv_scheme"] == "leave_one_fuel_out_cv"]["bias_kpa"].iloc[0]
    assert abs(loo_bias) < abs(raw_bias)
    assert abs(loo_bias) < 2.0
