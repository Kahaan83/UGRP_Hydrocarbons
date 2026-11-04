from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "A7" / "flash_margin_semantics.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "A7.md"


def test_flash_margin_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_flash_margin_semantics_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) >= 7
    required_cols = [
        "safety_margin_C", "analysis_level", "model_obs_coverage",
        "model_system_coverage", "uncovered_systems", "worst_case_violation_C"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify key finding: 3.4 C margin has ~89.7% observation coverage but only 83.9% system coverage
    row_34 = df[df["safety_margin_C"] == 3.4].iloc[0]
    assert 0.88 <= row_34["model_obs_coverage"] <= 0.92
    assert row_34["model_system_coverage"] < 0.85
    assert row_34["uncovered_systems"] == 5

    # Verify worst-case over-prediction exceeds 11 C at margin 0.5 C
    row_05 = df[df["safety_margin_C"] == 0.5].iloc[0]
    assert row_05["worst_case_violation_C"] > 10.0
