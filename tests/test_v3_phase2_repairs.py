# tests/test_v3_phase2_repairs.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_phase2_deliverables_exist():
    assert (ROOT / "results" / "v3" / "G2" / "octane_lump_model_eval.csv").exists()
    assert (ROOT / "results" / "v3" / "G6" / "dvpe_activity_eval.csv").exists()
    assert (ROOT / "configs" / "uncertainty_v3.yaml").exists()
    assert (ROOT / "results" / "v3" / "G8" / "diesel_cetane_cold_flow_eval.csv").exists()
    assert (ROOT / "results" / "v3" / "G10" / "thermoml_integration_eval.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "PHASE_2_SUMMARY.md").exists()

def test_octane_lump_model_improvement():
    df = pd.read_csv(ROOT / "results" / "v3" / "G2" / "octane_lump_model_eval.csv").set_index("model")
    mae_linear = float(df.loc["Linear-by-Mole Baseline", "ron_mae"])
    mae_ghosh = float(df.loc["Ghosh-Style Lump Model (LOO)", "ron_mae"])
    assert mae_ghosh < mae_linear
    assert mae_ghosh <= 3.5  # Substantial improvement toward literature 1-2 ON

def test_dvpe_activity_improvement():
    df = pd.read_csv(ROOT / "results" / "v3" / "G6" / "dvpe_activity_eval.csv").set_index("model")
    mae_ideal = float(df.loc["Ideal Raoult Baseline", "dvpe_mae_kpa"])
    mae_act = float(df.loc["Activity-Corrected DVPE Model (LOO)", "dvpe_mae_kpa"])
    assert mae_act < mae_ideal
