# tests/test_v3_r2_flash_anchoring.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r2_audit_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R2.md").exists()

def test_r2_strict_loso_metrics():
    df = pd.read_csv(ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv").set_index("evaluation_mode")
    assert "Strict LOSO (Pure-Anchored Baseline)" in df.index
    row = df.loc["Strict LOSO (Pure-Anchored Baseline)"]
    assert float(row["mae_C"]) < 3.0
    assert float(row["over_prediction_rate_gt_5C"]) <= 0.05
    assert float(row["max_over_prediction_C"]) < 10.0
