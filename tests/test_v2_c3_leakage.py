from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C3.md"
CHECKS_CSV = ROOT / "results" / "t0" / "leakage_checks.csv"


def test_c3_doc_exists():
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"
    assert CHECKS_CSV.exists(), f"Missing {CHECKS_CSV}"


def test_c3_leakage_checks_pass():
    df = pd.read_csv(CHECKS_CSV).set_index("check_id")

    # Core leakage checks must pass:
    # C1: 0 CRC rows in training
    assert df.loc["C1", "verdict"] == "PASS"
    assert int(df.loc["C1", "affected_rows"]) == 0

    # C2: 0 near-duplicates
    assert df.loc["C2", "verdict"] == "PASS"
    assert int(df.loc["C2", "affected_rows"]) == 0

    # C3b: 0 mixture measurements in Antoine refit
    assert df.loc["C3b", "verdict"] == "PASS"
    assert int(df.loc["C3b", "affected_rows"]) == 0

    # C5: T1 test suite passes
    assert df.loc["C5", "verdict"] == "PASS"
