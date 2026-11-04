"""Tests for Task D5: Error vs Number of Components Across Regimes."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "D5" / "error_vs_n.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "D5.md"


def test_d5_deliverables_exist():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_d5_buckets_and_regimes():
    df = pd.read_csv(CSV_PATH)
    expected_buckets = {
        "1 (pure)",
        "2 (binary)",
        "3 (ternary)",
        "4-10",
        "10-20",
        "20-40",
        "40+",
    }
    assert set(df["complexity_bucket"].unique()) == expected_buckets

    expected_regimes = {
        "synthetic_benchmark_blends",
        "measured_composition_fuels",
        "reconstructed_composition_fuels",
    }
    assert set(df["data_regime"].unique()) == expected_regimes


def test_d5_model_rank_inversion():
    """Verify the crucial finding: learned beats physics on binaries, but recalibrated physics beats learned on real fuels."""
    df = pd.read_csv(CSV_PATH)

    # 1. On lab binaries (N=2), learned residual outperforms raw physics
    b2 = df[df["complexity_bucket"] == "2 (binary)"].set_index("model")
    raw_b2 = b2.loc["Raw Physics (Refutas)", "log_mae"]
    ds_b2 = b2.loc["Learned Residual: DeepSets", "log_mae"]
    assert ds_b2 < raw_b2, f"Expected DeepSets {ds_b2} < Raw Physics {raw_b2} on binaries"

    # 2. On many-component real fuels (N=20-40), recalibrated physics wins overwhelmingly over DeepSets & Pairwise
    b_real = df[df["complexity_bucket"] == "20-40"].set_index("model")
    recal_real = b_real.loc["Recalibrated Physics (1000/T Slope)", "log_mae"]
    ds_real = b_real.loc["Learned Residual: DeepSets", "log_mae"]
    pw_real = b_real.loc["Learned Residual: Pairwise", "log_mae"]

    assert recal_real < 0.15, f"Expected recal_real < 0.15, got {recal_real}"
    assert ds_real > recal_real * 2.0, f"Expected DeepSets {ds_real} > 2x recalibrated {recal_real}"
    assert pw_real > recal_real * 2.0, f"Expected Pairwise {pw_real} > 2x recalibrated {recal_real}"
