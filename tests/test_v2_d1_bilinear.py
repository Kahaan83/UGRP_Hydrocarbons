"""Unit tests for Task D1: Real-Fuel Hindcast of the Bilinear Pairwise Model."""
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS_CSV = ROOT / "results" / "v2" / "D1" / "bilinear_hindcast.csv"
DOC_MD = ROOT / "docs" / "progress" / "v2" / "D1.md"


def test_deliverables_exist():
    assert RESULTS_CSV.exists(), f"Missing deliverable: {RESULTS_CSV}"
    assert DOC_MD.exists(), f"Missing progress document: {DOC_MD}"
    df = pd.read_csv(RESULTS_CSV)
    assert len(df) >= 20, f"Expected >= 20 rows, got {len(df)}"


def test_required_targets_and_conditions_present():
    df = pd.read_csv(RESULTS_CSV)
    targets = set(df["target"].unique())
    assert "kinematic_viscosity" in targets
    assert "density_15C_kg_m3" in targets
    assert "flash_point_C" in targets
    assert "nhoc_net_MJ_kg" in targets

    visc_df = df[df["target"] == "kinematic_viscosity"]
    conditions = set(visc_df["condition"].unique())
    assert "overall (all T)" in conditions
    assert "T = -40 C" in conditions
    assert "T = -20 C" in conditions
    assert "T = 20 C" in conditions


def test_bilinear_pairwise_superiority_on_viscosity():
    df = pd.read_csv(RESULTS_CSV)
    sub = df[(df["target"] == "kinematic_viscosity") & (df["condition"] == "overall (all T)")]
    
    raw_refutas = sub[sub["model"] == "Raw Refutas Physics"]["rmse_log"].values[0]
    bilinear_loo = sub[sub["model"] == "Bilinear Pairwise (LOO)"]["rmse_log"].values[0]
    deepsets = sub[sub["model"] == "DeepSets Residual"]["rmse_log"].values[0]
    pairwise = sub[sub["model"] == "Pairwise Residual"]["rmse_log"].values[0]

    # Bilinear Pairwise must beat Raw Refutas by > 60% error reduction
    err_reduction = (1.0 - bilinear_loo / raw_refutas) * 100.0
    assert err_reduction > 60.0, f"Expected >60% error reduction, got {err_reduction:.1f}%"

    # Bilinear Pairwise LOO must beat lab-only DeepSets and Pairwise
    assert bilinear_loo < deepsets, f"Bilinear {bilinear_loo} should outperform DeepSets {deepsets}"
    assert bilinear_loo < pairwise, f"Bilinear {bilinear_loo} should outperform Pairwise {pairwise}"


def test_bilinear_flash_point_offset_elimination():
    df = pd.read_csv(RESULTS_CSV)
    sub = df[df["target"] == "flash_point_C"]
    
    raw_flash = sub[sub["model"] == "Raw Physics (Le Chatelier + Raoult)"]["mae_log"].values[0]
    bilinear_flash = sub[sub["model"] == "Bilinear Pairwise (LOO)"]["mae_log"].values[0]

    assert raw_flash > 10.0, "Raw physics Le Chatelier should have significant offset on D86 pseudo-compositions"
    assert bilinear_flash < 3.0, f"Bilinear LOO should reduce flash MAE to < 3.0 C, got {bilinear_flash}"


def test_confidence_intervals_valid():
    df = pd.read_csv(RESULTS_CSV)
    for _, row in df.iterrows():
        ci_str = str(row["rmse_log_ci95"])
        assert ci_str.startswith("[") and ci_str.endswith("]")
        parts = [float(x.strip()) for x in ci_str[1:-1].split(",")]
        assert len(parts) == 2
        lo, hi = parts[0], parts[1]
        assert lo <= hi, f"Invalid CI: {lo} > {hi}"
        # Point estimate must lie within or very near bootstrap interval
        assert lo - 0.05 <= row["rmse_log"] <= hi + 0.05
