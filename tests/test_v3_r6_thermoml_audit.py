# tests/test_v3_r6_thermoml_audit.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r6_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R6" / "filter_drop_breakdown.csv").exists()
    assert (ROOT / "results" / "v3" / "R6" / "property_distributions_before_after.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R6.md").exists()

def test_r6_negative_excess_viscosity_audited():
    df = pd.read_csv(ROOT / "results" / "v3" / "R6" / "filter_drop_breakdown.csv").set_index("rule")
    assert "non_positive_value" in df.index
    content = (ROOT / "docs" / "progress" / "v3" / "R6.md").read_text()
    assert "excess viscosity" in content.lower()
    assert "absolute viscosity" in content.lower()
