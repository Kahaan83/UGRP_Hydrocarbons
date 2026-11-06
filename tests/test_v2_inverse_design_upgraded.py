from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DESIGNS_PATH = ROOT / "results" / "v2" / "improvement" / "upgraded_designs.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "upgraded_inverse_designs.md"


def test_upgraded_designs_exist():
    assert DESIGNS_PATH.exists()
    assert DOC_PATH.exists()


def test_all_upgraded_designs_pass_constraints():
    df = pd.read_csv(DESIGNS_PATH)
    # Check that status is PASS for all constraints across all 4 designs
    assert (df["status"] == "PASS").all()
    assert (df["margin_satisfied"] == True).all()

    # Check specific properties for Jet A 3-comp
    jet3 = df[df["design_id"] == "jet_minimal_3comp"].set_index("property")
    assert jet3.loc["flash_point_C", "computed_value"] >= 38.0
    assert jet3.loc["freezing_point_C", "computed_value"] <= -40.0
    assert jet3.loc["viscosity_minus20C_recal_cSt", "computed_value"] <= 8.0
