"""Test for T5.2 Recalibrated-Physics Baseline Ablation."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_CSV = ROOT / "results" / "t5" / "recalibrated_physics_ablation.csv"


def test_recalibrated_physics_results():
    assert OUT_CSV.exists(), f"Missing {OUT_CSV}"
    df = pd.read_csv(OUT_CSV)
    assert len(df) >= 4

    raw_phys = df[df["type"] == "physics"].iloc[0]
    recal = df[df["type"] == "recalibrated_physics"]

    # Recalibrated physics must outperform raw physics significantly
    for _, r in recal.iterrows():
        assert r["overall_rmse_ln"] < raw_phys["overall_rmse_ln"]
        assert r["pct_error_reduction"] > 50.0  # Recovers >50% error reduction
