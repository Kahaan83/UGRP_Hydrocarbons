"""Test for T5.1 System-Level Cluster Bootstrap Summary."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATH = ROOT / "results" / "t5" / "system_cluster_bootstrap_summary.csv"


def test_cluster_bootstrap_summary():
    assert SUMMARY_PATH.exists(), f"Missing {SUMMARY_PATH}"
    df = pd.read_csv(SUMMARY_PATH)
    assert len(df) >= 7

    # All rows must declare cluster bootstrap method
    for method in df["ci_method"]:
        assert "cluster bootstrap" in str(method).lower()

    # Flash LOSO must show statistically significant improvement
    flash_row = df[df["target"] == "flash_point_C"].iloc[0]
    assert flash_row["ci_95_hi"] < 0.0, "Flash LOSO CI upper bound must be negative (HGB better than baseline)"
