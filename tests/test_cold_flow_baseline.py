"""Unit test for Cold Flow Baselines & Scoping Decision (T2.9 - T2.12)."""
from pathlib import Path
import pandas as pd
import numpy as np

from fuelblend.baselines.cold_flow import (
    schroeder_van_laar_T,
    wax_appearance_temperature_C,
)

ROOT = Path(__file__).resolve().parents[1]
RESULTS_PATH = ROOT / "results" / "t2" / "cold_flow_baseline_eval.csv"
SCOPE_PATH = ROOT / "docs" / "scope_decisions.md"


def test_schroeder_van_laar_calculation():
    # Octane (Tm = 216.4 K, dHfus ~ 21.4 kJ/mol) at x = 1.0 (pure)
    t_pure = schroeder_van_laar_T(1.0, 216.4, 21400.0)
    assert np.isclose(t_pure, 216.4)

    # Diluted solute (x = 0.1) must precipitate at lower temperature
    t_dilute = schroeder_van_laar_T(0.1, 216.4, 21400.0)
    assert t_dilute < t_pure
    assert t_dilute > 150.0


def test_wax_appearance_temperature():
    # Binary of light solvent + heavy wax:
    # Solute A: Tm = 200 K, dH = 15 kJ/mol, x = 0.95
    # Solute B (Wax): Tm = 300 K, dH = 35 kJ/mol, x = 0.05
    x = [0.95, 0.05]
    tm = [200.0, 300.0]
    hf = [15000.0, 35000.0]

    wat_C = wax_appearance_temperature_C(x, tm, hf)
    assert wat_C > -100.0
    assert wat_C < 26.85  # Below pure wax Tm (300 K = 26.85 C)


def test_cold_flow_eval_table_exists():
    assert RESULTS_PATH.exists(), f"Missing {RESULTS_PATH}"
    df = pd.read_csv(RESULTS_PATH)
    assert len(df) == 1
    assert df["n_fuels"].iloc[0] == 57
    assert df["pearson_r"].iloc[0] > 0.3  # Positive correlation with freezing point


def test_scope_decision_pour_point():
    assert SCOPE_PATH.exists(), f"Missing {SCOPE_PATH}"
    content = SCOPE_PATH.read_text(encoding="utf-8")
    assert "Decision SD-01" in content
    assert "Pour point" in content
    assert "scoped out" in content.lower()
