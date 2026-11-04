"""Unit test for T2.15 DCN Blending Baselines."""
from pathlib import Path
import pandas as pd

from fuelblend.baselines.cetane import blend_dcn_linear

ROOT = Path(__file__).resolve().parents[1]
RESULTS_PATH = ROOT / "results" / "t2" / "cetane_baseline_grouped_cv.csv"


def test_blend_dcn_linear():
    pure = [20.0, 60.0]
    fracs = [0.5, 0.5]
    pred = blend_dcn_linear(fracs, pure)
    assert pred == 40.0


def test_cetane_grouped_cv_results():
    assert RESULTS_PATH.exists(), f"Missing {RESULTS_PATH}"
    df = pd.read_csv(RESULTS_PATH)
    assert len(df) == 3
    assert set(df["blending_basis"]) == {"mole", "mass", "volume"}

    mass_row = df[df["blending_basis"] == "mass"].iloc[0]
    mole_row = df[df["blending_basis"] == "mole"].iloc[0]

    # Mass basis outperforms mole basis
    assert mass_row["mae"] < mole_row["mae"]
    assert mass_row["mae"] < 4.0
