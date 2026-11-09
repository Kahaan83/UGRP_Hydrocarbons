"""QUARANTINED (v4 Task 3). Claim CLM-26 (R3 real-fuel flash hindcast, mixed 70-fuel set).

The 70-fuel hindcast mixes 57 CRC 647 fuels with 8 SURR_* and 5 POSF fuels whose compositions are
placeholder/target-derived in the ledger (dropped from v4). The CRC-only version of this check stays
active in tests/test_v3_phase0_regression.py.
"""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]

def test_r3_real_fuels_flash_regression():
    r3 = pd.read_csv(ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv")
    assert len(r3) >= 50
    # Pensky-Martens closed-cup negative bias check
    pm = r3[r3["method"].str.contains("D93|Pensky", na=False)]
    assert len(pm) >= 40
    assert pm["error_C"].mean() < -5.0  # conservative negative bias
