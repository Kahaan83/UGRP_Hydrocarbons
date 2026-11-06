from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "real_fuels" / "candidates_v2.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "B1.md"


def test_candidates_inventory_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_candidates_inventory_counts():
    df = pd.read_csv(CSV_PATH)
    required_cols = [
        "candidate_id", "fuel_name", "fuel_family", "source",
        "composition_method", "resolved_components_est", "properties_reported",
        "access_type", "license", "ingestion_status"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify at least 15 to 20 fuels per family
    counts = df["fuel_family"].value_counts()
    assert counts.get("jet", 0) >= 20
    assert counts.get("gasoline", 0) >= 15
    assert counts.get("diesel", 0) >= 15
    assert len(df) >= 50
