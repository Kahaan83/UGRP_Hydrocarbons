# tests/test_v3_r4_surrogates.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r4_surrogate_check_exists():
    assert (ROOT / "results" / "v3" / "R4" / "surrogate_check.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R4.md").exists()

def test_r4_all_11_surrogates_audited():
    df = pd.read_csv(ROOT / "results" / "v3" / "R4" / "surrogate_check.csv")
    surrs = df[df["is_published_surrogate"] == True]
    assert len(surrs) == 11
    # Check that failed constraints are identified
    assert "all_constraints_passed" in df.columns
    assert "failed_constraints" in df.columns
    # All 11 published surrogates must fail at least one ASTM constraint
    assert (surrs["all_constraints_passed"] == False).all()
    # Designed blends must pass
    designs = df[df["is_published_surrogate"] == False]
    assert len(designs) == 3
    assert (designs["all_constraints_passed"] == True).all()

