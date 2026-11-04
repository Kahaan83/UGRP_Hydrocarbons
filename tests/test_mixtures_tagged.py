"""Test for T1.4 tagged mixture table."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAGGED_PATH = ROOT / "data" / "processed" / "mixtures_tagged.csv"


def test_tagged_mixtures_counts():
    assert TAGGED_PATH.exists(), f"Missing {TAGGED_PATH}"
    df = pd.read_csv(TAGGED_PATH)
    assert len(df) == 510, f"Expected 510 tagged systems, got {len(df)}"
    assert df["system_id"].nunique() == 510, "System IDs must be unique"

    # All systems must have valid fuel_family
    allowed_families = {"jet", "diesel", "gasoline", "model-mixture"}
    assert set(df["fuel_family"].unique()).issubset(allowed_families)
    assert df["fuel_family"].isna().sum() == 0

    # real_fuel must be boolean
    assert set(df["real_fuel"].unique()).issubset({True, False})
    assert df["real_fuel"].isna().sum() == 0
