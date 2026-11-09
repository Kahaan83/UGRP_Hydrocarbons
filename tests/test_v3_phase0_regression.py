# tests/test_v3_phase0_regression.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_claims_register_v3_entries():
    content = (ROOT / "docs" / "claims_register.md").read_text()
    for clm in ["CLM-12", "CLM-13", "CLM-14", "CLM-15", "CLM-16", "CLM-17", "CLM-18"]:
        assert clm in content, f"Missing {clm} in docs/claims_register.md"

def test_r2_flash_anchoring_regression():
    r2 = pd.read_csv(ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv").set_index("evaluation_mode")
    row = r2.loc["Strict LOSO (Pure-Anchored Baseline)"]
    assert abs(float(row["mae_C"]) - 2.28) < 0.2
    assert float(row["over_prediction_rate_gt_5C"]) <= 0.05
    assert float(row["max_over_prediction_C"]) <= 5.5

def test_r3_real_fuels_flash_crc647_subset():
    """CRC 647 aviation fuels only (measured D93 flash + reconstructed compositions, kept in v4)."""
    r3 = pd.read_csv(ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv")
    crc = r3[r3["fuel_id"].astype(str).str.startswith("CRC647")]
    assert len(crc) == 57
    pm = crc[crc["method"].str.contains("D93|Pensky", na=False)]
    assert len(pm) >= 40
    assert pm["error_C"].mean() < -5.0  # conservative negative bias

def test_r1_viscosity_leaderboard_regression():
    r1 = pd.read_csv(ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv")
    recal = r1[(r1["model"] == "Recalibrated Refutas (2-Param LOO)") & (r1["condition"] == "overall (all T)")].iloc[0]
    bilinear = r1[(r1["model"] == "Bilinear Pairwise (LOO)") & (r1["condition"] == "overall (all T)")].iloc[0]
    assert abs(float(recal["log_rmse"]) - 0.1411) < 0.005
    assert abs(float(bilinear["log_rmse"]) - 0.1242) < 0.005

def test_r4_surrogates_regression():
    r4 = pd.read_csv(ROOT / "results" / "v3" / "R4" / "surrogate_check.csv")
    surrs = r4[r4["is_published_surrogate"] == True]
    assert len(surrs) == 11
    assert (surrs["all_constraints_passed"] == False).all()

def test_r5_freeze_regimes_regression():
    r5 = pd.read_csv(ROOT / "results" / "v3" / "R5" / "freeze_regime_test.csv").set_index("regime")
    assert float(r5.loc["Simple Blends (2-5 comps)", "mae_wat_minus_14C"]) > 8.0
    assert float(r5.loc["Simple Blends (2-5 comps)", "mae_coutinho_eutectic"]) < float(r5.loc["Simple Blends (2-5 comps)", "mae_wat_minus_14C"])

def test_r6_thermoml_audit_regression():
    r6 = pd.read_csv(ROOT / "results" / "v3" / "R6" / "filter_drop_breakdown.csv").set_index("rule")
    assert int(r6.loc["non_positive_value", "rows_matched"]) == 60

def test_r7_hoc_clarity_regression():
    r7 = pd.read_csv(ROOT / "results" / "v3" / "R7" / "hoc_propagation_error.csv").set_index("evaluation_setting")
    assert float(r7.loc["Group Additivity (Pure Structure Basis)", "propagated_input_uncertainty_MJ_kg"]) == 0.0
    assert float(r7.loc["ASTM D3338 (Predicted Density & Aromatics at Design Time)", "propagated_input_uncertainty_MJ_kg"]) > 0.05
