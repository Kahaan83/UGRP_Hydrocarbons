"""QUARANTINED (v4 Task 3). Claims CLM-23 (U6 v3 inverse designs) and CLM-24 (U2 "30/30 robustness" check).

The v3 designs were produced by scripts/v3/run_v3_inverse_designs.py with its own flash/freeze/DVPE
models and margins that disagree with configs/uncertainty_v3.yaml (CLM-28); the tests below pin those
designs. The purchasable-pool integrity test does not depend on them and stays active in
tests/test_v3_phase3_designs.py.
"""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]

def test_phase3_deliverables_exist():
    assert (ROOT / "data" / "processed" / "v3" / "purchasable_compounds.csv").exists()
    assert (ROOT / "configs" / "specs_astm.yaml").exists()
    assert (ROOT / "configs" / "specs_bis.yaml").exists()
    assert (ROOT / "results" / "v3" / "U6" / "designs_v3.csv").exists()
    assert (ROOT / "results" / "v3" / "U6" / "design_robustness_audit.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "PHASE_3_SUMMARY.md").exists()

def test_v3_designs_pass_constraints():
    df = pd.read_csv(ROOT / "results" / "v3" / "U6" / "designs_v3.csv")
    assert len(df) >= 3
    for _, row in df.iterrows():
        assert row["all_constraints_passed"] == True or str(row["all_constraints_passed"]).lower() == "true"
