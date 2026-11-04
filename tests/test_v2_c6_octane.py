from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "C6" / "octane_real_fuel.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C6.md"


def test_octane_real_fuel_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_octane_real_fuel_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 6
    required_cols = [
        "target_property", "target_name", "model", "cv_scheme",
        "n_fuels", "mae", "mae_ci95", "rmse", "rmse_ci95", "bias"
    ]
    for col in required_cols:
        assert col in df.columns

    # Check that both RON and MON are reported
    props = set(df["target_property"])
    assert "ron" in props
    assert "mon" in props

    # Check that recalibrated model reduces MAE
    ron_df = df[df["target_property"] == "ron"].set_index("model")
    mae_lin = ron_df.loc["Linear-by-Volume Octane Blending", "mae"]
    mae_recal = ron_df.loc["Recalibrated BON Baseline (LOO CV)", "mae"]
    assert mae_recal < mae_lin

    # Bias should be close to zero under LOO recalibration
    bias_recal = ron_df.loc["Recalibrated BON Baseline (LOO CV)", "bias"]
    assert abs(bias_recal) < 1.0
