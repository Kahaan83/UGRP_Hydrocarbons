from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
ABLATION_CSV = ROOT / "results" / "v2" / "improvement" / "paper_ablation_table.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "THERMODYNAMIC_IMPROVEMENTS_SUMMARY.md"


def test_ablation_files_exist():
    assert ABLATION_CSV.exists()
    assert DOC_PATH.exists()


def test_ablation_monotonic_error_reduction():
    df = pd.read_csv(ABLATION_CSV).set_index("version")
    assert len(df) == 3

    # Flash point progression: v1 > v2 > v2_thermo
    assert df.loc["v1_raw_physics", "flash_mae_C"] > df.loc["v2_recalibrated", "flash_mae_C"] > df.loc["v2_thermodynamic_anchored", "flash_mae_C"]
    assert df.loc["v1_raw_physics", "flash_max_overprediction_C"] > df.loc["v2_recalibrated", "flash_max_overprediction_C"] > df.loc["v2_thermodynamic_anchored", "flash_max_overprediction_C"]

    # Freezing point progression: v1 > v2 > v2_thermo
    assert df.loc["v1_raw_physics", "freeze_wat_mae_C"] > df.loc["v2_recalibrated", "freeze_wat_mae_C"] > df.loc["v2_thermodynamic_anchored", "freeze_wat_mae_C"]

    # Cryogenic viscosity progression: v1 > v2 > v2_thermo
    assert df.loc["v1_raw_physics", "visc_m40C_log_rmse"] > df.loc["v2_recalibrated", "visc_m40C_log_rmse"] > df.loc["v2_thermodynamic_anchored", "visc_m40C_log_rmse"]

    # Gasoline RON progression: v1 > v2 > v2_thermo
    assert df.loc["v1_raw_physics", "gasoline_ron_mae"] > df.loc["v2_recalibrated", "gasoline_ron_mae"] > df.loc["v2_thermodynamic_anchored", "gasoline_ron_mae"]
