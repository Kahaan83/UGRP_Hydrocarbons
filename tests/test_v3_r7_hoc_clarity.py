from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_r7_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R7" / "hoc_propagation_error.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R7.md").exists()


def test_r7_design_time_metrics():
    df = pd.read_csv(ROOT / "results" / "v3" / "R7" / "hoc_propagation_error.csv").set_index("evaluation_setting")
    assert "Group Additivity (Pure Structure Basis)" in df.index
    assert "ASTM D3338 (Measured Density & Aromatics)" in df.index
    assert "ASTM D3338 (Predicted Density & Aromatics at Design Time)" in df.index
    # Group additivity has no propagated input error
    ga_prop = float(df.loc["Group Additivity (Pure Structure Basis)", "propagated_input_uncertainty_MJ_kg"])
    assert ga_prop == 0.0
    # D3338 with predicted inputs must show non-zero propagated uncertainty
    d3338_prop = float(df.loc["ASTM D3338 (Predicted Density & Aromatics at Design Time)", "propagated_input_uncertainty_MJ_kg"])
    assert d3338_prop > 0.05
