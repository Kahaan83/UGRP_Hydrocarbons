# tests/test_v3_phase3_designs.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_purchasable_pool_integrity():
    df = pd.read_csv(ROOT / "data" / "processed" / "v3" / "purchasable_compounds.csv")
    assert "cas_number" in df.columns
    assert "supplier" in df.columns
    assert "purity_pct" in df.columns
    assert "is_single_isomer" in df.columns
    single_isomers = df[df["is_single_isomer"] == True]
    assert len(single_isomers) >= 12
