from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "A5" / "hoc_recheck.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "A5_hoc.md"


def test_hoc_recheck_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_hoc_recheck_content():
    df = pd.read_csv(CSV_PATH)
    assert len(df) == 3
    required_cols = [
        "model", "evaluation", "n_fuels", "features",
        "mae_MJ_kg", "mae_ci95", "rmse_MJ_kg", "rmse_ci95", "mape_pct", "bias_MJ_kg"
    ]
    for col in required_cols:
        assert col in df.columns

    evals = df.set_index("evaluation")
    mae_phys = evals.loc["zero_parameters", "mae_MJ_kg"]
    mae_in = evals.loc["in_sample_reproduced", "mae_MJ_kg"]
    mae_loo = evals.loc["leave_one_fuel_out_cv", "mae_MJ_kg"]

    # Verify physical ordering: in-sample < LOO < physics
    assert mae_in < mae_loo < mae_phys
    assert abs(mae_phys - 0.6889) < 0.001
    assert abs(mae_in - 0.0783) < 0.001
    assert abs(mae_loo - 0.0860) < 0.002
    # Verify significant error reduction (> 80%) under LOO CV
    assert (mae_phys - mae_loo) / mae_phys > 0.80
