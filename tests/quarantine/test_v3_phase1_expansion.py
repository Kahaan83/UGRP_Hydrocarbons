"""QUARANTINED (v4 Task 3). Claim CLM-29 (S1/S2/S4 Phase-1 corpus counts and coverage matrix).

The gasoline/diesel counts and the S4 coverage matrix count fuels of the unfiltered v3 inventory,
most of which have target-derived or placeholder compositions (results/v4/P1/provenance_ledger.csv;
55 gasolines, 11 diesels and 18 jets are dropped from data/real_fuels/v4). The split-leakage and
oxygenate-registry tests are unaffected and stay active in tests/test_v3_phase1_expansion.py.
"""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]

def test_phase1_deliverables_exist():
    assert (ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv").exists()
    assert (ROOT / "data" / "real_fuels" / "v3" / "compositions_v3.csv").exists()
    assert (ROOT / "results" / "v3" / "S4" / "coverage_matrix.csv").exists()
    assert (ROOT / "splits" / "v3" / "real_fuels_v3_splits.json").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "PHASE_1_SUMMARY.md").exists()

def test_s1_gasoline_count():
    inv = pd.read_csv(ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv")
    gasolines = inv[inv["fuel_family"] == "gasoline"]
    assert len(gasolines) >= 30, f"Expected >= 30 gasolines, got {len(gasolines)}"
    assert "ethanol_vol_pct" in inv.columns or "oxygenate_pct" in inv.columns

def test_s2_diesel_count():
    inv = pd.read_csv(ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv")
    diesels = inv[inv["fuel_family"] == "diesel"]
    assert len(diesels) >= 15, f"Expected >= 15 diesels, got {len(diesels)}"
    # Verify cetane labeling exists
    assert "cetane_test_type" in inv.columns or any("cetane" in c.lower() for c in inv.columns)

def test_s4_coverage_matrix():
    cov = pd.read_csv(ROOT / "results" / "v3" / "S4" / "coverage_matrix.csv")
    assert "fuels_with_composition_and_measurement" in cov.columns
    assert "status" in cov.columns
    # Check that flags are present if count < 10
    below_10 = cov[cov["fuels_with_composition_and_measurement"] < 10]
    if len(below_10) > 0:
        assert (below_10["status"] == "FLAG_BELOW_10").all()

