"""Unit tests for T2.5 - T2.8 Flash Point Expansion and Activity Coefficient Baselines."""
from pathlib import Path
import pandas as pd

from fuelblend.data.validate import validate_measurement_df

ROOT = Path(__file__).resolve().parents[1]
EXP_FLASH_CSV = ROOT / "data" / "processed" / "measurements_flash_expanded.csv"
ACT_CSV = ROOT / "results" / "t2" / "flash_loso_activity_coeff.csv"
METHOD_CSV = ROOT / "results" / "t2" / "flash_loso_by_method.csv"


def test_flash_expanded_schema():
    assert EXP_FLASH_CSV.exists(), f"Missing {EXP_FLASH_CSV}"
    df = pd.read_csv(EXP_FLASH_CSV)
    assert len(df) >= 5000

    # Strict validator must pass
    val = validate_measurement_df(df, strict=True)
    assert val["valid"] is True

    # Method must be filled for every row
    assert df["method"].isna().sum() == 0
    assert (df["method"].str.strip() != "").all()


def test_activity_coefficient_baseline():
    assert ACT_CSV.exists(), f"Missing {ACT_CSV}"
    df = pd.read_csv(ACT_CSV)
    assert len(df) == 2

    raoult = df[df["baseline_model"].str.contains("Raoult")].iloc[0]
    act = df[df["baseline_model"].str.contains("Activity Coeff")].iloc[0]

    # Activity coefficient model reduces RMSE and dangerous over-prediction rate
    assert act["rmse_C"] < raoult["rmse_C"]
    assert act["overpred_gt_5C_pct"] < raoult["overpred_gt_5C_pct"]


def test_per_method_evaluation():
    assert METHOD_CSV.exists(), f"Missing {METHOD_CSV}"
    df = pd.read_csv(METHOD_CSV)
    assert len(df) >= 1
    assert "test_method" in df.columns
    assert "act_coeff_rmse_C" in df.columns
