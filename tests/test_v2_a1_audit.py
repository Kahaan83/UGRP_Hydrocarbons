"""Tests for Task A1: Independent Audit of Inverse Designs."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
AUDIT_CSV = ROOT / "results" / "v2" / "A1" / "design_audit.csv"
PROGRESS_MD = ROOT / "docs" / "progress" / "v2" / "A1.md"


def test_a1_audit_files_exist():
    assert AUDIT_CSV.exists(), "results/v2/A1/design_audit.csv must exist"
    assert PROGRESS_MD.exists(), "docs/progress/v2/A1.md must exist"


def test_a1_audit_constraints_evaluated():
    df = pd.read_csv(AUDIT_CSV)
    required_designs = {"jet_csv", "jet_claimed_text", "diesel_csv", "gasoline_csv"}
    assert set(df["design_id"]).issuperset(required_designs)

    # Check Jet A claimed text
    jet_claimed = df[df["design_id"] == "jet_claimed_text"]
    props = dict(zip(jet_claimed["property"], jet_claimed["passed"]))

    # Must verify freezing point / WAT fails
    assert not props["wat_C"], "Dodecane + tetralin must FAIL freezing point <= -40 C"
    # Must verify viscosity at -20 C fails
    assert not props["kin_visc_cSt@-20C"], "Dodecane + tetralin must FAIL viscosity at -20 C due to solidification"
    # Density must pass
    assert props["density_kg_m3"], "Dodecane + tetralin passes density 775-840 kg/m3"

    # Gasoline CSV must fail density and freezing point
    gas_df = df[df["design_id"] == "gasoline_csv"]
    gas_props = dict(zip(gas_df["property"], gas_df["passed"]))
    assert not gas_props["density_kg_m3"], "Gasoline candidate pool must fail density limit"
    assert not gas_props["wat_C"], "Gasoline blend must fail freezing point"
