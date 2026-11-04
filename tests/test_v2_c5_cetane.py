from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CURATED_CSV = ROOT / "data" / "processed" / "cetane_pure_curated.csv"
MISSING_MD = ROOT / "data" / "MISSING.md"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "C5.md"


def test_cetane_c5_outputs_exist():
    assert CURATED_CSV.exists(), f"Missing {CURATED_CSV}"
    assert MISSING_MD.exists(), f"Missing {MISSING_MD}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_cetane_c5_content():
    df = pd.read_csv(CURATED_CSV)
    assert len(df) >= 88

    # Verify primary reference standards exist
    iks = set(df["inchikey"])
    # 1-methylnaphthalene (DCN = 0.0)
    row_1mn = df[df["compound_name"].str.contains("1-methylnaphthalene", na=False)]
    assert len(row_1mn) > 0
    assert float(row_1mn["cetane_number"].iloc[0]) == 0.0

    # Isooctane
    row_iso = df[df["compound_name"].str.contains("isooctane", na=False)]
    assert len(row_iso) > 0
    assert float(row_iso["cetane_number"].iloc[0]) in (15.0, 17.5)

    # Verify missing gaps logged in MISSING.md
    missing_text = MISSING_MD.read_text(encoding="utf-8")
    assert "NREL Cetane Compendium" in missing_text
    assert "cetane_pure_curated.csv" in missing_text
