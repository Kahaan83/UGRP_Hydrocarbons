"""
tests/test_v2_a5_distillation.py

Tests for Task A5: Distillation validation dataset and baseline evaluation.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

BASE = Path(__file__).resolve().parents[1]
REAL_FUELS_CSV = BASE / "data" / "processed" / "v2" / "distillation_real_fuels.csv"
BASELINE_CSV = BASE / "results" / "v2" / "A5" / "distillation_baseline.csv"


def test_distillation_real_fuels_deliverable_exists():
    assert REAL_FUELS_CSV.exists(), f"Missing deliverable: {REAL_FUELS_CSV}"
    df = pd.read_csv(REAL_FUELS_CSV)
    assert len(df) == 303, f"Expected 303 curves, found {len(df)}"
    
    required_cols = [
        "fuel_id",
        "source",
        "fuel_name",
        "fuel_class",
        "method",
        "T0_C",
        "T10_C",
        "T20_C",
        "T50_C",
        "T90_C",
        "FBP_C",
    ]
    for c in required_cols:
        assert c in df.columns, f"Missing required column: {c}"
        
    # Check that all 4 sources are present
    sources = set(df["source"].unique())
    assert sources == {"crc647", "av2718", "fortin2012", "ajftd"}


def test_distillation_curve_monotonicity():
    df = pd.read_csv(REAL_FUELS_CSV)
    cuts = ["T0_C", "T10_C", "T20_C", "T50_C", "T90_C", "FBP_C"]
    
    # Check that where cuts are present, temperature is strictly or weakly monotonically increasing
    monotonic_count = 0
    for _, r in df.iterrows():
        vals = [r[c] for c in cuts if pd.notna(r[c])]
        if len(vals) >= 3:
            diffs = np.diff(vals)
            # Allow tiny numerical tolerance of 0.5 C for experimental read noise
            assert (diffs >= -0.5).all(), f"Distillation curve not monotonic for {r['fuel_id']}: {vals}"
            monotonic_count += 1
    assert monotonic_count >= 200, f"Expected at least 200 verifiable curves, got {monotonic_count}"


def test_distillation_baseline_deliverable_exists():
    assert BASELINE_CSV.exists(), f"Missing deliverable: {BASELINE_CSV}"
    df = pd.read_csv(BASELINE_CSV)
    assert len(df) > 0, "Baseline evaluation CSV is empty"
    
    # Must have both D2887 vs D86 and RP1 surrogate evaluations
    eval_types = set(df["evaluation_type"].unique())
    assert "D2887_ideal_vs_D86_actual" in eval_types
    assert "RP1_surrogate_ideal_vs_ADC_actual" in eval_types


def test_distillation_thermodynamic_gap():
    df = pd.read_csv(BASELINE_CSV)
    d86_eval = df[df["evaluation_type"] == "D2887_ideal_vs_D86_actual"]
    
    # T50 should have relatively small discrepancy between D86 and D2887 (< 6.0 C MAE)
    t50 = d86_eval[d86_eval["recovery_point"] == "T50"]
    assert np.mean(np.abs(t50["error_C"])) < 6.0, f"Expected T50 MAE < 6.0 C, got {np.mean(np.abs(t50['error_C']))}"
    
    # T0 (IBP) and FBP should show massive systematic bias (> 30.0 C MAE)
    t0 = d86_eval[d86_eval["recovery_point"] == "T0"]
    fbp = d86_eval[d86_eval["recovery_point"] == "FBP"]
    assert np.mean(np.abs(t0["error_C"])) > 35.0, f"Expected IBP MAE > 35.0 C, got {np.mean(np.abs(t0['error_C']))}"
    assert np.mean(np.abs(fbp["error_C"])) > 35.0, f"Expected FBP MAE > 35.0 C, got {np.mean(np.abs(fbp['error_C']))}"
    
    # Experimental repeatability should be tight (< 2.0 C MAE across all cuts)
    rep_diffs = d86_eval["d86_repeatability_diff_C"].dropna()
    assert np.mean(np.abs(rep_diffs)) < 2.0, f"Expected repeatability MAE < 2.0 C, got {np.mean(np.abs(rep_diffs))}"
