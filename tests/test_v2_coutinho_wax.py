from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
EVAL_PATH = ROOT / "results" / "v2" / "improvement" / "freezing_point_coutinho_eval.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "freezing_point_coutinho_upgrade.md"


def test_coutinho_files_exist():
    assert EVAL_PATH.exists()
    assert DOC_PATH.exists()


def test_coutinho_eutectic_depression():
    from fuelblend.thermo.coutinho_wax import predict_freezing_point_coutinho, get_pure_thermal_properties
    thermal_map = get_pure_thermal_properties()
    
    # Binary n-C10 + n-C16 mixture (10 wt% n-C16, 90 wt% n-C10)
    # n-C10 inchikey: DIOQZVSQGTUSAI-UHFFFAOYSA-N, n-C16 inchikey: DCAYPVUWAIABOU-UHFFFAOYSA-N
    iks = ["DIOQZVSQGTUSAI-UHFFFAOYSA-N", "DCAYPVUWAIABOU-UHFFFAOYSA-N"]
    tf, delta_e, meta = predict_freezing_point_coutinho([0.90, 0.10], iks, thermal_map)
    
    # Verify eutectic depression occurs (delta_e > 5 C)
    assert delta_e > 5.0
    assert tf < meta["ideal_wat_C"]


def test_coutinho_eval_error_reduction():
    df = pd.read_csv(EVAL_PATH).set_index("model")
    mae_ideal = df.loc["Ideal Single-Solute WAT (Schroeder-van Laar)", "mae_C"]
    mae_cout = df.loc["Coutinho Multi-Component Non-Ideal Wax Model", "mae_C"]
    
    # Coutinho model must reduce MAE by > 65% over raw ideal WAT
    assert (mae_ideal - mae_cout) / mae_ideal > 0.65
    assert mae_cout < 4.5
