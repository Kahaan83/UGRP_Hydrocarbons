from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "A12" / "redistribution_audit.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "A12.md"


def test_redistribution_audit_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_redistribution_audit_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) >= 24
    required_cols = [
        "source_id", "source_name", "license_terms", "redistribution_flag",
        "tracked_files_in_git", "tracked_file_count", "risk_assessment", "recommended_action"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify that AJFTD is flagged for action required
    ajftd_row = df[df["source_id"] == "ajftd"].iloc[0]
    assert "ACTION_REQUIRED" in ajftd_row["risk_assessment"]
    assert ajftd_row["tracked_file_count"] >= 1
    assert "ajftd_data" in ajftd_row["tracked_files_in_git"]

    # Verify that open sources are compliant
    thermoml = df[df["source_id"] == "thermoml"].iloc[0]
    assert "COMPLIANT" in thermoml["risk_assessment"]
    assert thermoml["tracked_file_count"] == 0
