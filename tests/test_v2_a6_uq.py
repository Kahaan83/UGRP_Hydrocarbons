from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "A6" / "uncertainty_growth_refit.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "A6_uq.md"


def test_uncertainty_growth_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_uncertainty_growth_refit_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 4
    required_cols = [
        "property", "predictor", "dataset", "v1_claimed_g",
        "v1_claimed_coverage", "v1_actual_coverage_real_fuels",
        "refit_g_in_sample_90", "refit_g_loo_mean", "held_out_loo_coverage", "status"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify that claimed coverage for DeepSets on real fuels failed severely
    ds_visc = df[(df["property"] == "viscosity_kin_cSt") & (df["predictor"] == "DeepSets Neural Residual")].iloc[0]
    assert ds_visc["v1_actual_coverage_real_fuels"] < 0.10  # Actual coverage ~6.4%
    assert ds_visc["refit_g_loo_mean"] > 40.0              # Needs ~50x inflation

    ds_dens = df[(df["property"] == "density_kg_m3") & (df["predictor"] == "DeepSets Neural Residual")].iloc[0]
    assert ds_dens["v1_actual_coverage_real_fuels"] < 0.20  # Actual coverage ~15.8%
    assert ds_dens["refit_g_loo_mean"] > 10.0              # Needs ~14.5x inflation

    # Verify that LOO refit achieves target ~90% coverage
    for _, row in df.iterrows():
        assert 0.85 <= row["held_out_loo_coverage"] <= 0.95
