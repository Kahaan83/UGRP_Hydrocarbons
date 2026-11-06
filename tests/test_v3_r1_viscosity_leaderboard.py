from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r1_leaderboard_exists():
    assert (ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R1.md").exists()

def test_r1_models_present():
    df = pd.read_csv(ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv")
    models = set(df["model"].unique())
    expected = [
        "Raw Refutas (Zero-Param Physics)",
        "Recalibrated Refutas (2-Param LOO)",
        "Walther (ASTM D341 Raw Physics)",
        "Walther + Recalibration (LOO)",
        "Bilinear Pairwise (LOO)",
        "Stacked Physics + Residual",
        "DeepSets Neural Residual",
        "Pairwise Neural Residual"
    ]
    for m in expected:
        assert m in models
    assert "log_rmse" in df.columns
    assert "log_rmse_ci95" in df.columns
    assert "training_data" in df.columns
    assert "fold_discipline" in df.columns
