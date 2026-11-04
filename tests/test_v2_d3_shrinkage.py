"""Unit tests for Task D3: Shrinkage Sweep on Real Fuels."""
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS_CSV = ROOT / "results" / "v2" / "D3" / "shrinkage_sweep.csv"
DOC_MD = ROOT / "docs" / "progress" / "v2" / "D3.md"


def test_deliverables_exist():
    assert RESULTS_CSV.exists(), f"Missing deliverable: {RESULTS_CSV}"
    assert DOC_MD.exists(), f"Missing progress document: {DOC_MD}"
    df = pd.read_csv(RESULTS_CSV)
    assert len(df) >= 50, f"Expected >= 50 rows in sweep, got {len(df)}"


def test_properties_and_alpha_range():
    df = pd.read_csv(RESULTS_CSV)
    targets = set(df["target"].unique())
    assert "kinematic_viscosity" in targets
    assert "density_15C_kg_m3" in targets
    assert "flash_point_C" in targets

    for tgt in targets:
        sub = df[df["target"] == tgt]
        alphas = set(np.round(sub["alpha_shrinkage"], 2))
        assert 0.0 in alphas, f"alpha=0.0 missing for {tgt}"
        assert 1.0 in alphas, f"alpha=1.0 missing for {tgt}"
        assert 0.5 in alphas, f"alpha=0.5 missing for {tgt}"


def test_loo_error_evaluation_not_just_training_residual():
    df = pd.read_csv(RESULTS_CSV)
    assert "loo_rmse" in df.columns
    assert "loo_mae" in df.columns
    assert "t4_mean_abs_residual" in df.columns

    # In T4.9, residual scaled trivially with alpha
    # In Task D3, LOO-RMSE reflects true out-of-sample error
    for tgt in ["kinematic_viscosity", "flash_point_C"]:
        sub = df[(df["target"] == tgt) & (df["condition"].str.contains("overall|closed"))]
        rmse_a0 = sub[sub["alpha_shrinkage"] == 0.0]["loo_rmse"].values[0]
        rmse_best = sub["loo_rmse"].min()
        assert rmse_best < rmse_a0 * 0.5, f"Expected >50% improvement at optimal alpha for {tgt}"


def test_viscosity_temperature_dependent_shrinkage():
    df = pd.read_csv(RESULTS_CSV)
    # Check -40 C vs +20 C optimal alpha
    sub_m40 = df[(df["target"] == "kinematic_viscosity") & (df["condition"] == "T = -40 C")]
    sub_p20 = df[(df["target"] == "kinematic_viscosity") & (df["condition"] == "T = 20 C")]

    opt_m40 = sub_m40.loc[sub_m40["loo_rmse"].idxmin()]["alpha_shrinkage"]
    opt_p20 = sub_p20.loc[sub_p20["loo_rmse"].idxmin()]["alpha_shrinkage"]

    # At -40 C, larger alpha is needed to fix underprediction than at +20 C
    assert opt_m40 > opt_p20, f"Expected opt_m40 ({opt_m40}) > opt_p20 ({opt_p20})"


def test_confidence_intervals_valid():
    df = pd.read_csv(RESULTS_CSV)
    for _, row in df.iterrows():
        ci_str = str(row["loo_rmse_ci95"])
        assert ci_str.startswith("[") and ci_str.endswith("]")
        parts = [float(x.strip()) for x in ci_str[1:-1].split(",")]
        assert len(parts) == 2
        lo, hi = parts[0], parts[1]
        assert lo <= hi, f"Invalid CI: {lo} > {hi}"
        assert lo - 0.05 <= row["loo_rmse"] <= hi + 0.05
