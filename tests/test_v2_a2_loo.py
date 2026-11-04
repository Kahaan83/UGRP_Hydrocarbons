"""Tests for Task A2: Leave-Fuel-Out Recalibration Check."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
LOO_CSV = ROOT / "results" / "v2" / "A2" / "recalibration_loo.csv"
PROGRESS_MD = ROOT / "docs" / "progress" / "v2" / "A2.md"


def test_a2_files_exist():
    assert LOO_CSV.exists(), "results/v2/A2/recalibration_loo.csv must exist"
    assert PROGRESS_MD.exists(), "docs/progress/v2/A2.md must exist"


def test_a2_models_evaluated():
    df = pd.read_csv(LOO_CSV)
    models = set(df["model"])
    assert "Physics (Raw Refutas)" in models
    assert "Recalibrated Physics (2-param 1000/T, LOO)" in models
    assert "Recalibrated Physics (3-param per-T bias, LOO)" in models
    assert "Recalibrated Physics (2-param, fit on lab only)" in models
    assert "Learned Residual: DeepSets (Attention)" in models
    assert "Learned Residual: Pairwise" in models

    # Check LOO recalibration metrics
    loo_row = df[df["model"] == "Recalibrated Physics (2-param 1000/T, LOO)"].iloc[0]
    assert 0.13 < loo_row["overall_rmse_ln"] < 0.15, "LOO 2-parameter RMSE must be ~0.1411"
    assert loo_row["pct_error_reduction"] > 65.0, "LOO error reduction must exceed 65%"

    # Check lab-only fit fails to reduce error significantly
    lab_row = df[df["model"] == "Recalibrated Physics (2-param, fit on lab only)"].iloc[0]
    assert lab_row["overall_rmse_ln"] > 0.38, "Lab-only fit cannot achieve low RMSE on real cold fuels"
