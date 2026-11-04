"""Unit test for T2.4 Heat of Combustion Real-Fuel Evaluation."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EVAL_CSV = ROOT / "results" / "t2" / "heat_of_combustion_eval.csv"


def test_heat_of_combustion_eval_results():
    assert EVAL_CSV.exists(), f"Missing {EVAL_CSV}"
    df = pd.read_csv(EVAL_CSV)
    assert len(df) == 2

    base_row = df[df["model"].str.startswith("Physics Baseline")].iloc[0]
    assert base_row["fuel_family"] == "jet"
    assert base_row["n_fuels"] == 57
    assert base_row["mape_pct"] < 2.0  # Physics baseline MAPE is under 2%
    assert base_row["mae_MJ_kg"] < 1.0  # MAE is under 1 MJ/kg
