# tests/test_v3_r5_freeze_regimes.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r5_freeze_regimes_exist():
    assert (ROOT / "results" / "v3" / "R5" / "freeze_regime_test.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R5.md").exists()

def test_r5_regime_breakdown():
    df = pd.read_csv(ROOT / "results" / "v3" / "R5" / "freeze_regime_test.csv").set_index("regime")
    assert "Simple Blends (2-5 comps)" in df.index
    assert "Multicomponent Real Fuels (>20 comps)" in df.index
    assert "mae_wat_minus_14C" in df.columns
    assert "mae_coutinho_eutectic" in df.columns
