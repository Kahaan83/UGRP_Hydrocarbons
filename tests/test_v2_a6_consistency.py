"""
tests/test_v2_a6_consistency.py

Tests for Task A6: Multi-property correlation audit and automated consistency checker.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from fuelblend.eval.consistency import PhysicalConsistencyChecker, ConsistencyResult

BASE = Path(__file__).resolve().parents[1]
CORR_CSV = BASE / "results" / "v2" / "A6" / "property_correlations.csv"


def test_property_correlations_deliverable_exists():
    assert CORR_CSV.exists(), f"Missing deliverable: {CORR_CSV}"
    df = pd.read_csv(CORR_CSV)
    assert len(df) == 121, f"Expected 121 correlation pairs (11x11), found {len(df)}"
    
    required_cols = [
        "property_1",
        "property_2",
        "n_observations",
        "pearson_r",
        "pearson_p_value",
        "spearman_rho",
        "spearman_p_value",
    ]
    for c in required_cols:
        assert c in df.columns, f"Missing required column: {c}"


def test_real_fuel_correlation_magnitudes():
    df = pd.read_csv(CORR_CSV)
    
    # T50 vs T90 distillation correlation should be positive and substantial (> 0.60)
    sub_50_90 = df[(df["property_1"] == "distillation_T50_C") & (df["property_2"] == "distillation_T90_C")]
    assert len(sub_50_90) == 1
    r_50_90 = float(sub_50_90["pearson_r"].iloc[0])
    assert r_50_90 > 0.60, f"Expected strong correlation between T50 and T90, got {r_50_90}"
    
    # Flash point vs T10 distillation should be positive (> 0.30) across 2,000+ fuels
    sub_fl_t10 = df[(df["property_1"] == "flash_point_C") & (df["property_2"] == "distillation_T10_C")]
    assert len(sub_fl_t10) == 1
    r_fl_t10 = float(sub_fl_t10["pearson_r"].iloc[0])
    assert r_fl_t10 > 0.30, f"Expected positive flash vs T10 correlation, got {r_fl_t10}"


def test_consistency_checker_valid_fuel():
    checker = PhysicalConsistencyChecker()
    # Typical Jet A properties
    valid_jet_a = {
        "flash_point_C": 48.0,
        "distillation_T10_C": 175.0,
        "distillation_T0_C": 155.0,
        "freezing_point_C": -47.0,
        "viscosity_minus20C_cSt": 4.5,
        "density_15C_kg_m3": 804.0,
        "aromatics_vol_pct": 18.0,
        "dcn": 45.0,
    }
    res = checker.check_candidate_blend(valid_jet_a)
    assert res.is_consistent is True, f"Valid Jet A flagged with violations: {res.violations}"
    assert len(res.violations) == 0


def test_consistency_checker_catches_flash_higher_than_t10():
    checker = PhysicalConsistencyChecker()
    # Flashing at 180 C with T10 = 170 C (thermodynamically impossible)
    invalid_fuel = {
        "flash_point_C": 180.0,
        "distillation_T10_C": 170.0,
        "freezing_point_C": -50.0,
    }
    res = checker.check_candidate_blend(invalid_fuel)
    assert res.is_consistent is False
    assert any("flash point" in v.lower() and "t10" in v.lower() for v in res.violations)


def test_consistency_checker_catches_frozen_viscosity_contradiction():
    checker = PhysicalConsistencyChecker()
    # Fuel freezes at -10 C, so claiming a liquid viscosity at -20 C is a physical contradiction
    invalid_fuel = {
        "freezing_point_C": -10.0,
        "viscosity_minus20C_cSt": 5.0,
    }
    res = checker.check_candidate_blend(invalid_fuel)
    assert res.is_consistent is False
    assert any("solid/wax" in v.lower() or "freezes at" in v.lower() for v in res.violations)


def test_consistency_checker_catches_high_aromatic_high_cetane():
    checker = PhysicalConsistencyChecker()
    # 40% aromatics with DCN = 65 (combustion kinetics violation)
    invalid_fuel = {
        "aromatics_vol_pct": 40.0,
        "dcn": 65.0,
    }
    res = checker.check_candidate_blend(invalid_fuel)
    assert res.is_consistent is False
    assert any("combustion kinetics" in v.lower() for v in res.violations)
