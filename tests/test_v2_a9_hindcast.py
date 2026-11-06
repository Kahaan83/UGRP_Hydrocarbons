"""Tests for Task A9: Deterministic Hindcast."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
HIND_V2 = ROOT / "results" / "v2" / "A9" / "hindcast_predictions.csv"
PROGRESS_MD = ROOT / "docs" / "progress" / "v2" / "A9.md"
REPRO_SH = ROOT / "scripts" / "reproduce_hindcast.sh"


def test_a9_files_exist():
    assert HIND_V2.exists(), "results/v2/A9/hindcast_predictions.csv must exist"
    assert PROGRESS_MD.exists(), "docs/progress/v2/A9.md must exist"
    assert REPRO_SH.exists(), "scripts/reproduce_hindcast.sh must exist"


def test_a9_hindcast_predictions_match_baseline():
    baseline = pd.read_csv(ROOT / "results" / "hindcast" / "hindcast_predictions.csv")
    reproduced = pd.read_csv(HIND_V2)

    assert len(baseline) == len(reproduced) == 228
    # Assert primary model predictions match within numerical precision (< 1e-3)
    for col in ["physics", "extrap_deepsets", "extrap_pairwise"]:
        diff = (baseline[col] - reproduced[col]).abs().max()
        assert diff < 1e-3, f"Column {col} diff {diff} exceeds numerical precision tolerance"
