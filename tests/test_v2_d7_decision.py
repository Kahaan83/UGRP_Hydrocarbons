"""Tests for Task D7: Master Model Strategy Decision Note."""
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
DECISION_MD = ROOT / "docs" / "model_decision.md"


def test_d7_decision_doc_exists():
    assert DECISION_MD.exists(), f"Missing {DECISION_MD}"


def test_d7_required_table_citations():
    text = DECISION_MD.read_text(encoding="utf-8")

    # Ground rule: Every conclusion must point to a specific deliverable table
    required_tables = [
        "results/v2/A2/recalibration_loo.csv",
        "results/v2/A4/gn_baselines.csv",
        "results/v2/A6/uncertainty_growth_refit.csv",
        "results/v2/C1/hoc_eval.csv",
        "results/v2/C2/flash_loso_v2.csv",
        "results/v2/C4/cold_flow_v2.csv",
        "results/v2/D4/stacked_eval.csv",
        "results/v2/D5/error_vs_n.csv",
        "results/v2/D6/truncation.csv",
    ]

    for table in required_tables:
        assert table in text, f"Decision note missing citation to {table}"


def test_d7_primary_model_decisions():
    text = DECISION_MD.read_text(encoding="utf-8")

    # Primary model designations must be unequivocal
    assert "2-Parameter Recalibrated Refutas" in text
    assert "Ideal Volume Blending" in text
    assert "Group Additivity GCV Baseline" in text
    assert "Transparent Offset SLE WAT" in text

    # Architectural scope rule: neural residuals scoped to lab mixtures with shrinkage
    assert "ambient" in text.lower()
    assert "shrinkage" in text.lower()
    assert "recalibrated physics" in text.lower()
