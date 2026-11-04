from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "C2" / "flash_loso_v2.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C2.md"


def test_flash_loso_v2_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_flash_loso_v2_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 6
    required_cols = [
        "model", "has_learned_residual", "cv_scheme", "n_obs", "n_systems",
        "mae_C", "rmse_C", "bias_C", "over_prediction_rate_gt_5C", "worst_case_over_prediction_C"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify that learned residual reduces MAE on Raoult
    df_models = df.set_index("model")
    mae_raoult_base = df_models.loc["Raoult (Ideal, gamma=1)", "mae_C"]
    mae_raoult_res = df_models.loc["Raoult (Ideal, gamma=1) + HGB Residual", "mae_C"]
    assert mae_raoult_res < mae_raoult_base
    assert abs(mae_raoult_res - 3.9748) < 0.05

    # Verify that UNIFAC + residual caps worst-case over-prediction lower than Raoult base
    max_unifac_res = df_models.loc["UNIFAC-type (Group Contribution) + HGB Residual", "worst_case_over_prediction_C"]
    max_raoult_base = df_models.loc["Raoult (Ideal, gamma=1)", "worst_case_over_prediction_C"]
    assert max_unifac_res < max_raoult_base
    assert max_unifac_res < 30.0
