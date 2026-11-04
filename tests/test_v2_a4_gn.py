from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "A4" / "gn_baselines.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "A4_gn.md"


def test_gn_baselines_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_gn_baselines_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 5
    required_cols = [
        "model", "n_obs", "n_systems", "fold_discipline",
        "log_rmse", "log_mae", "mape_pct", "n_parameters"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify key hierarchy:
    # Per-system G < Refutas < Class-Pair G < Global G < GN (G=0)
    models = df.set_index("model")
    rmse_g0 = models.loc["Grunberg-Nissan (Ideal Log-Linear, G=0)", "log_rmse"]
    rmse_g_global = models.loc["Grunberg-Nissan (Global G, 5-Fold Grouped CV)", "log_rmse"]
    rmse_g_class = models.loc["Grunberg-Nissan (Class-Pair G, 5-Fold Grouped CV)", "log_rmse"]
    rmse_refutas = models.loc["Refutas Blending Index (Zero Parameters)", "log_rmse"]
    rmse_sys = models.loc["Grunberg-Nissan (Per-System G, In-Sample Bound)", "log_rmse"]

    assert rmse_g0 > rmse_g_global > rmse_g_class > rmse_refutas > rmse_sys
    assert abs(rmse_refutas - 0.0699) < 0.001
    assert abs(rmse_g0 - 0.1014) < 0.001
