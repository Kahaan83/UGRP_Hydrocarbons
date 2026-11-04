"""
tests/test_v2_a4_flash.py

Tests for Task A4: Flash point baseline (Le Chatelier) full evaluation.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from fuelblend.baselines.flash import (
    predict_flash_le_chatelier,
    psat_antoine,
    lfl_temperature_correction,
    check_antoine_validity,
)

BASE = Path(__file__).resolve().parents[1]
OUT_CSV = BASE / "results" / "v2" / "A4" / "flash_baseline_full.csv"


def test_flash_baseline_full_deliverable_exists():
    assert OUT_CSV.exists(), f"Missing deliverable: {OUT_CSV}"
    df = pd.read_csv(OUT_CSV)
    assert len(df) == 1006, f"Expected 1,006 mixture records, found {len(df)}"
    
    required_cols = [
        "row_id",
        "system_id",
        "n_comp",
        "is_hc",
        "coverage_status",
        "pred_flash_const_lfl_C",
        "err_const_lfl_C",
        "pred_flash_temp_lfl_C",
        "err_temp_lfl_C",
        "antoine_in_range_const",
        "has_c12_plus",
    ]
    for c in required_cols:
        assert c in df.columns, f"Missing required column: {c}"


def test_flash_baseline_coverage_counts():
    df = pd.read_csv(OUT_CSV)
    ok_df = df[df["coverage_status"] == "ok"]
    assert len(ok_df) == 126, f"Expected 126 predictable hydrocarbon records, found {len(ok_df)}"
    
    # All 126 predictable records must have is_hc == True
    assert ok_df["is_hc"].all(), "Predictable records must all be hydrocarbons"
    
    # 880 non-hydrocarbon mixtures must be flagged as missing components
    non_ok = df[df["coverage_status"] != "ok"]
    assert len(non_ok) == 1006 - 126


def test_flash_baseline_solver_accuracy_unit():
    # Sanity check: pure decane (A=6.199, B=1580.177, C=-70.46, LFL=0.75 vol%)
    p_decane, stat = predict_flash_le_chatelier([1.0], [(6.199, 1580.177, -70.46)], [0.75])
    assert stat == "ok"
    assert abs(p_decane - 47.4) < 0.5, f"Expected ~47.4 C for decane, got {p_decane}"


def test_flash_baseline_temperature_lfl_direction():
    # As T increases above 25 C, LFL(T) decreases, so flash point should be lower with temp correction
    p_const, _ = predict_flash_le_chatelier([1.0], [(6.199, 1580.177, -70.46)], [0.75], lfl_temp=False)
    p_temp, _ = predict_flash_le_chatelier([1.0], [(6.199, 1580.177, -70.46)], [0.75], lfl_temp=True)
    assert p_temp < p_const, "Temperature-dependent LFL should yield lower flash point for T > 25 C"


def test_sub_c12_high_accuracy():
    df = pd.read_csv(OUT_CSV)
    ok_df = df[df["coverage_status"] == "ok"]
    non_c12 = ok_df[~ok_df["has_c12_plus"]]
    assert len(non_c12) == 48, f"Expected 48 non-C12+ mixtures, got {len(non_c12)}"
    
    mae_non_c12 = float(np.mean(np.abs(non_c12["err_const_lfl_C"])))
    assert mae_non_c12 < 2.0, f"Expected MAE < 2.0 C on <C12 mixtures, got {mae_non_c12}"


def test_dangerous_overpredictions_flagged():
    df = pd.read_csv(OUT_CSV)
    ok_df = df[df["coverage_status"] == "ok"]
    
    # Check that dangerous over-predictions (pred > measured) are recorded
    dangerous = ok_df[ok_df["err_const_lfl_C"] > 0]
    assert len(dangerous) > 0, "Expected some dangerous over-predictions"
    max_dangerous_err = float(ok_df["err_const_lfl_C"].max())
    assert max_dangerous_err > 35.0, f"Expected max dangerous error > 35 C, got {max_dangerous_err}"
