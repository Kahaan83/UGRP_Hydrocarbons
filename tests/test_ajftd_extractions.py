"""Unit test for AJFTD smoke point and distillation extractions (T22 and T28)."""
from pathlib import Path
import pandas as pd

def test_smoke_point():
    p = Path("data/processed/smoke_point_ajftd.csv")
    assert p.exists(), "smoke_point_ajftd.csv missing"
    df = pd.read_csv(p)
    assert len(df) == 100, f"Expected 100 smoke points, got {len(df)}"
    assert set(["fuel_id", "fuel_label", "value", "unit", "method", "tier"]).issubset(df.columns)
    assert (df["unit"] == "mm").all()
    assert (df["method"] == "D1322").all()
    assert (df["tier"] == "A").all()
    # Smoke points must be within physical jet fuel range (15 - 35 mm)
    assert (df["value"] >= 15.0).all() and (df["value"] <= 35.0).all()
    print("test_smoke_point passed.")

def test_distillation():
    p = Path("data/processed/distillation_ajftd.csv")
    assert p.exists(), "distillation_ajftd.csv missing"
    df = pd.read_csv(p)
    assert len(df) == 99, f"Expected 99 distillation curves, got {len(df)}"
    expected_cols = ["fuel_id", "T0_C", "T10_C", "T20_C", "T30_C", "T40_C", "T50_C", "T60_C", "T70_C", "T80_C", "T90_C", "T100_C", "method", "tier"]
    for col in expected_cols:
        assert col in df.columns, f"Missing column {col}"
    assert (df["method"] == "D86").all()
    assert (df["tier"] == "A").all()
    # T100 < 320 C for all fuels
    assert (df["T100_C"] < 320.0).all()
    # Monotonicity holds for 98/99 fuels (1 known data anomaly flagged in fuel 433)
    assert df["is_monotonic"].sum() == 98
    print("test_distillation passed.")

if __name__ == "__main__":
    test_smoke_point()
    test_distillation()
    print("All AJFTD extraction tests passed successfully.")
