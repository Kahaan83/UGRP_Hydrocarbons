"""Unit test for Consistency-Check Distillation (T2.20, T2.21)."""
from pathlib import Path
import pandas as pd

from fuelblend.data.validate import validate_measurement_df

ROOT = Path(__file__).resolve().parents[1]
DIST_CSV = ROOT / "data" / "processed" / "measurements_distillation.csv"
CHECK_CSV = ROOT / "results" / "t2" / "distillation_consistency_check.csv"


def test_distillation_measurements_validation():
    assert DIST_CSV.exists(), f"Missing {DIST_CSV}"
    df = pd.read_csv(DIST_CSV)
    assert len(df) == 342  # 57 fuels * 6 cuts
    val = validate_measurement_df(df, strict=True)
    assert val["valid"] is True


def test_distillation_consistency_check():
    assert CHECK_CSV.exists(), f"Missing {CHECK_CSV}"
    df = pd.read_csv(CHECK_CSV)
    assert len(df) == 57
    assert (df["distillation_slope_t90_t10"] > 0).all()
