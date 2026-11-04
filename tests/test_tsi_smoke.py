"""
test_tsi_smoke.py - Test suite for TSI and smoke point reference table (T19, T20, T42).
"""
import os
import json
import pandas as pd

def test_tsi_table_exists():
    path = "data/processed/tsi_table.csv"
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert len(df) >= 30, f"Expected >= 30 rows, got {len(df)}"

def test_tsi_columns():
    df = pd.read_csv("data/processed/tsi_table.csv")
    expected_cols = [
        "inchikey", "smiles", "compound_name", "tsi", "ysi",
        "smoke_point_mm", "source", "tier", "derived", "censored", "note"
    ]
    for col in expected_cols:
        assert col in df.columns, f"Missing column {col}"

def test_tsi_anchors():
    df = pd.read_csv("data/processed/tsi_table.csv")
    
    # Check Methylcyclohexane anchor (TSI = 5.0)
    mch = df[df["compound_name"] == "methylcyclohexane"]
    assert len(mch) >= 1
    assert abs(mch.iloc[0]["tsi"] - 5.0) < 1e-2

    # Check 1-Methylnaphthalene anchor (TSI = 100.0)
    mn = df[df["compound_name"] == "1-methylnaphthalene"]
    assert len(mn) >= 1
    assert abs(mn.iloc[0]["tsi"] - 100.0) < 1e-2

def test_tsi_ranges():
    df = pd.read_csv("data/processed/tsi_table.csv")
    
    # Toluene TSI ~ 36-39
    tol = df[df["compound_name"] == "toluene"]
    assert len(tol) >= 1
    assert 36.0 <= tol.iloc[0]["tsi"] <= 39.5

    # Benzene TSI ~ 25-29
    ben = df[df["compound_name"] == "benzene"]
    assert len(ben) >= 1
    assert 25.0 <= ben.iloc[0]["tsi"] <= 29.0

    # n-alkanes (hexane ~ 2.5, heptane ~ 2.8)
    hex_row = df[df["compound_name"] == "n-hexane"]
    assert len(hex_row) >= 1
    assert abs(hex_row.iloc[0]["tsi"] - 2.5) < 0.2

    hep_row = df[df["compound_name"] == "n-heptane"]
    assert len(hep_row) >= 1
    assert abs(hep_row.iloc[0]["tsi"] - 2.8) < 0.2

    # Check smoke points are positive
    valid_sp = df["smoke_point_mm"].dropna()
    assert (valid_sp > 0).all()

    # Check tier values
    assert set(df["tier"].unique()).issubset({"A", "B", "C"})

def test_tsi_constants_json():
    path = "data/processed/tsi_constants.json"
    assert os.path.exists(path), f"Missing {path}"
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    assert "boehm2022" in data
    assert "mensch2010" in data
    assert abs(data["boehm2022"]["a"] - (-3.40)) < 1e-3
    assert abs(data["boehm2022"]["b"] - 3.727) < 1e-3
    assert abs(data["mensch2010"]["a"] - (-4.75)) < 1e-3
    assert abs(data["mensch2010"]["b"] - 4.051) < 1e-3

if __name__ == "__main__":
    test_tsi_table_exists()
    test_tsi_columns()
    test_tsi_anchors()
    test_tsi_ranges()
    test_tsi_constants_json()
    print("All tests in test_tsi_smoke.py passed successfully!")
