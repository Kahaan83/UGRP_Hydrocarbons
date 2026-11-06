# tests/test_v3_phase3_designs.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_phase3_deliverables_exist():
    assert (ROOT / "data" / "processed" / "v3" / "purchasable_compounds.csv").exists()
    assert (ROOT / "configs" / "specs_astm.yaml").exists()
    assert (ROOT / "configs" / "specs_bis.yaml").exists()
    assert (ROOT / "results" / "v3" / "U6" / "designs_v3.csv").exists()
    assert (ROOT / "results" / "v3" / "U6" / "design_robustness_audit.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "PHASE_3_SUMMARY.md").exists()

def test_purchasable_pool_integrity():
    df = pd.read_csv(ROOT / "data" / "processed" / "v3" / "purchasable_compounds.csv")
    assert "cas_number" in df.columns
    assert "supplier" in df.columns
    assert "purity_pct" in df.columns
    assert "is_single_isomer" in df.columns
    single_isomers = df[df["is_single_isomer"] == True]
    assert len(single_isomers) >= 12

def test_v3_designs_pass_constraints():
    df = pd.read_csv(ROOT / "results" / "v3" / "U6" / "designs_v3.csv")
    assert len(df) >= 3
    for _, row in df.iterrows():
        assert row["all_constraints_passed"] == True or str(row["all_constraints_passed"]).lower() == "true"
