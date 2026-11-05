# tests/test_v3_r3_flash_hindcast.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r3_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R3.md").exists()

def test_r3_real_fuels_flash_metrics():
    df = pd.read_csv(ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv")
    assert len(df) >= 10
    assert "method" in df.columns
    assert "error_C" in df.columns
    assert "fuel_family" in df.columns
    methods = df["method"].unique()
    assert any("D93" in m or "Pensky" in m for m in methods)
