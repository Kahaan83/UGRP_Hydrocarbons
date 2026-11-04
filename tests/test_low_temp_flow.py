"""Test suite for low temperature flow properties of real fuels and blends."""

import os
import pandas as pd

PROCESSED_FILE = os.path.join(
    os.path.dirname(__file__), "..", "data", "processed", "low_temp_flow_properties_fuels.csv"
)

def test_file_exists_and_non_empty():
    assert os.path.exists(PROCESSED_FILE), f"Missing {PROCESSED_FILE}"
    df = pd.read_csv(PROCESSED_FILE)
    assert len(df) > 0, "File is empty"
    assert len(df) >= 150, f"Expected >= 150 records, got {len(df)}"

def test_columns_and_schema():
    df = pd.read_csv(PROCESSED_FILE)
    expected_cols = [
        "fuel_id", "fuel_name", "source", "property", "value",
        "unit", "qualifier", "method", "wax_content_pct",
        "n_paraffin_content_pct", "tier", "derived", "note"
    ]
    assert list(df.columns) == expected_cols, f"Columns mismatch: {list(df.columns)}"

def test_units_and_properties():
    df = pd.read_csv(PROCESSED_FILE)
    assert (df["unit"] == "C").all(), "All units must be Celsius ('C')"
    valid_props = {"pour_point_C", "freezing_point_C", "cloud_point_C", "cfpp_C"}
    assert set(df["property"].unique()).issubset(valid_props), f"Unexpected properties: {df['property'].unique()}"

def test_physically_sound_bounds():
    df = pd.read_csv(PROCESSED_FILE)
    # Jet and diesel fuels flow properties should be between -100 °C and +20 °C
    assert (df["value"] >= -100.0).all(), f"Values below -100 C found: {df[df['value'] < -100]}"
    assert (df["value"] <= 20.0).all(), f"Values above +20 C found: {df[df['value'] > 20]}"

def test_sources_represented():
    df = pd.read_csv(PROCESSED_FILE)
    sources = set(df["source"].unique())
    # Should include CRC, AJFTD, NASA, and Kirgina
    assert any("crc" in s.lower() for s in sources), "Missing CRC source"
    assert "ajftd" in sources, "Missing AJFTD source"
    assert any("19860018789" in s for s in sources), "Missing NASA NTRS source"
    assert any("kirgina" in s.lower() for s in sources), "Missing Kirgina source"

if __name__ == "__main__":
    test_file_exists_and_non_empty()
    test_columns_and_schema()
    test_units_and_properties()
    test_physically_sound_bounds()
    test_sources_represented()
    print("All low temp flow tests passed successfully!")
