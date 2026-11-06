"""Unit tests for Task D2: Class-Level Interaction Learning."""
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS_CSV = ROOT / "results" / "v2" / "D2" / "class_interactions.csv"
DOC_MD = ROOT / "docs" / "progress" / "v2" / "D2.md"


def test_deliverables_exist():
    assert RESULTS_CSV.exists(), f"Missing deliverable: {RESULTS_CSV}"
    assert DOC_MD.exists(), f"Missing progress document: {DOC_MD}"
    df = pd.read_csv(RESULTS_CSV)
    assert len(df) >= 20, f"Expected >= 20 rows, got {len(df)}"


def test_required_models_present():
    df = pd.read_csv(RESULTS_CSV)
    models = set(df["model"].unique())
    assert "Raw Refutas Physics" in models
    assert any("Lab Class Interaction" in m for m in models)
    assert any("Real-Fuel Class Interaction" in m for m in models)
    assert any("Recalibrated Physics (Task A2" in m for m in models)
    assert any("Recalibrated Physics + Class Interaction" in m for m in models)


def test_recalibrated_physics_dominates_lab_class_interactions():
    df = pd.read_csv(RESULTS_CSV)
    sub = df[df["condition"] == "overall (all T)"]
    
    raw = sub[sub["model"] == "Raw Refutas Physics"]["rmse_log"].values[0]
    lab_trans = sub[sub["model"].str.contains("Lab Class Interaction")]["rmse_log"].values[0]
    recal = sub[sub["model"].str.contains(r"Recalibrated Physics \(Task A2")]["rmse_log"].values[0]
    comb = sub[sub["model"].str.contains(r"Recalibrated Physics \+ Class Interaction")]["rmse_log"].values[0]

    # Lab class interactions should not improve over raw refutas on real fuels (zero-shot failure)
    assert lab_trans >= raw * 0.95, "Lab class interactions should fail to resolve cold-flow slope zero-shot"

    # Recalibrated physics should strongly outperform lab transfer
    assert recal < 0.15, f"Recalibrated physics should achieve < 0.15 log-RMSE, got {recal}"
    assert comb <= recal + 1e-4, f"Combined model should match or refine recalibrated physics: {comb} vs {recal}"


def test_confidence_intervals_valid():
    df = pd.read_csv(RESULTS_CSV)
    for _, row in df.iterrows():
        ci_str = str(row["rmse_log_ci95"])
        assert ci_str.startswith("[") and ci_str.endswith("]")
        parts = [float(x.strip()) for x in ci_str[1:-1].split(",")]
        assert len(parts) == 2
        lo, hi = parts[0], parts[1]
        assert lo <= hi, f"Invalid CI: {lo} > {hi}"
        assert lo - 0.05 <= row["rmse_log"] <= hi + 0.05
