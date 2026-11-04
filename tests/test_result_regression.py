"""Task A10: Result Regression Test Suite.

Pins all verified Phase A headline results within tight numerical tolerances.
CI fails if any headline number drifts without an explicit claims register update.
"""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_a1_inverse_design_audit_regression():
    path = ROOT / "results" / "v2" / "A1" / "design_audit.csv"
    assert path.exists()
    df = pd.read_csv(path)

    # 1. Jet claimed dodecane + tetralin must fail freezing point (wat_C) and viscosity
    jet_claimed = df[df["design_id"] == "jet_claimed_text"].set_index("property")
    assert jet_claimed.loc["wat_C", "passed"] is False or str(jet_claimed.loc["wat_C", "passed"]).lower() == "false"
    assert abs(float(jet_claimed.loc["wat_C", "computed_value"]) - (-15.2)) < 0.2
    assert jet_claimed.loc["kin_visc_cSt@-20C", "passed"] is False or str(jet_claimed.loc["kin_visc_cSt@-20C", "passed"]).lower() == "false"
    assert "SOLID" in str(jet_claimed.loc["kin_visc_cSt@-20C", "computed_value"])

    # 2. Jet CSV design (hexadecane + biphenyl) freezes at +20.1 C
    jet_csv = df[df["design_id"] == "jet_csv"].set_index("property")
    assert abs(float(jet_csv.loc["wat_C", "computed_value"]) - 20.1) < 0.2


def test_a2_recalibration_loo_regression():
    path = ROOT / "results" / "v2" / "A2" / "recalibration_loo.csv"
    assert path.exists()
    df = pd.read_csv(path).set_index("model")

    raw_rmse = df.loc["Physics (Raw Refutas)", "overall_rmse_ln"]
    loo_rmse = df.loc["Recalibrated Physics (2-param 1000/T, LOO)", "overall_rmse_ln"]
    lab_rmse = df.loc["Recalibrated Physics (2-param, fit on lab only)", "overall_rmse_ln"]

    assert abs(raw_rmse - 0.4421) < 0.005
    assert abs(loo_rmse - 0.1411) < 0.005
    assert abs(lab_rmse - 0.4090) < 0.005


def test_a3_cetane_basis_regression():
    path = ROOT / "results" / "v2" / "A3" / "cetane_basis_rerun.csv"
    assert path.exists()
    df = pd.read_csv(path)

    vol_mae = df[df["blending_basis"].str.startswith("volume")]["mae"].iloc[0]
    mass_mae = df[df["blending_basis"].str.startswith("mass")]["mae"].iloc[0]
    mole_mae = df[df["blending_basis"].str.startswith("mole")]["mae"].iloc[0]

    assert abs(vol_mae - 3.488) < 0.01
    assert abs(mass_mae - 3.731) < 0.01
    assert abs(mole_mae - 5.707) < 0.01
    assert vol_mae < mass_mae < mole_mae


def test_a4_gn_baselines_regression():
    path = ROOT / "results" / "v2" / "A4" / "gn_baselines.csv"
    assert path.exists()
    df = pd.read_csv(path).set_index("model")

    gn_zero = df.loc["Grunberg-Nissan (Ideal Log-Linear, G=0)", "log_rmse"]
    refutas = df.loc["Refutas Blending Index (Zero Parameters)", "log_rmse"]
    gn_class = df.loc["Grunberg-Nissan (Class-Pair G, 5-Fold Grouped CV)", "log_rmse"]
    gn_sys = df.loc["Grunberg-Nissan (Per-System G, In-Sample Bound)", "log_rmse"]

    assert abs(gn_zero - 0.1014) < 0.002
    assert abs(refutas - 0.0699) < 0.002
    assert abs(gn_class - 0.0890) < 0.002
    assert abs(gn_sys - 0.0412) < 0.002


def test_a5_hoc_recheck_regression():
    path = ROOT / "results" / "v2" / "A5" / "hoc_recheck.csv"
    assert path.exists()
    df = pd.read_csv(path).set_index("evaluation")

    phys_mae = df.loc["zero_parameters", "mae_MJ_kg"]
    loo_mae = df.loc["leave_one_fuel_out_cv", "mae_MJ_kg"]

    assert abs(phys_mae - 0.6889) < 0.002
    assert abs(loo_mae - 0.0860) < 0.005


def test_a6_uq_growth_refit_regression():
    path = ROOT / "results" / "v2" / "A6" / "uncertainty_growth_refit.csv"
    assert path.exists()
    df = pd.read_csv(path)

    ds_visc = df[(df["property"] == "viscosity_kin_cSt") & (df["predictor"] == "DeepSets Neural Residual")].iloc[0]
    assert ds_visc["v1_actual_coverage_real_fuels"] < 0.10
    assert ds_visc["refit_g_loo_mean"] > 40.0


def test_a7_flash_margin_semantics_regression():
    path = ROOT / "results" / "v2" / "A7" / "flash_margin_semantics.csv"
    assert path.exists()
    df = pd.read_csv(path)

    row_34 = df[df["safety_margin_C"] == 3.4].iloc[0]
    assert 0.88 <= row_34["model_obs_coverage"] <= 0.92
    assert row_34["model_system_coverage"] < 0.85
    assert row_34["uncovered_systems"] == 5


def test_a8_flash_data_usage_regression():
    path = ROOT / "results" / "v2" / "A8" / "flash_usage_audit.csv"
    assert path.exists()
    df = pd.read_csv(path)

    kinds = df[df["category"] == "sample_kind"].set_index("sample_kind")
    assert kinds.loc["pure", "n_observations"] == 2726
    assert kinds.loc["real_fuel", "n_observations"] == 2487
    assert kinds.loc["binary", "n_observations"] + kinds.loc["multi", "n_observations"] == 126


def test_a9_deterministic_hindcast_regression():
    path = ROOT / "results" / "v2" / "A9" / "hindcast_summary.csv"
    assert path.exists()
    df = pd.read_csv(path)

    row_m20 = df[(df["kind"] == "kv") & (df["T_C"] == -20.0)].iloc[0]
    m20_rmse = float(row_m20["extrap_deepsets_rmse_ln"])
    assert abs(m20_rmse - 0.2720) < 0.005
