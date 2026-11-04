"""Tests for Task A3: Cetane Basis Bug Check & Physics Verification."""
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from fuelblend.baselines.cetane import blend_dcn_linear, mass_to_volume_fractions, predict_mixture_dcn

ROOT = Path(__file__).resolve().parents[1]
RERUN_CSV = ROOT / "results" / "v2" / "A3" / "cetane_basis_rerun.csv"
PROGRESS_MD = ROOT / "docs" / "progress" / "v2" / "A3.md"


def test_a3_files_exist():
    assert RERUN_CSV.exists(), "results/v2/A3/cetane_basis_rerun.csv must exist"
    assert PROGRESS_MD.exists(), "docs/progress/v2/A3.md must exist"


def test_mass_to_volume_conversion():
    # If densities are equal, volume fractions must equal mass fractions
    w = [0.5, 0.5]
    rho_equal = [800.0, 800.0]
    v_eq = mass_to_volume_fractions(w, rho_equal)
    np.testing.assert_allclose(v_eq, [0.5, 0.5])

    # If densities differ (e.g. heptane 684 vs tetralin 970), volume fractions must differ
    rho_diff = [684.0, 970.0]
    v_diff = mass_to_volume_fractions(w, rho_diff)
    # Lighter heptane occupies more volume (> 50%)
    assert v_diff[0] > 0.55
    assert v_diff[1] < 0.45
    assert np.isclose(v_diff.sum(), 1.0)


def test_a3_mass_and_volume_differ_and_volume_is_superior():
    df = pd.read_csv(RERUN_CSV).set_index("blending_basis")

    mass_mae = df.loc["mass", "mae"]
    vol_mae = df.loc["volume (true liquid density)", "mae"]
    mass_rmse = df.loc["mass", "rmse"]
    vol_rmse = df.loc["volume (true liquid density)", "rmse"]

    # Acceptance criterion: Mass and volume results differ unless a test proves they must be identical
    assert abs(mass_mae - vol_mae) > 0.1, "Mass and volume MAE must differ with real densities"
    assert abs(mass_rmse - vol_rmse) > 0.1, "Mass and volume RMSE must differ with real densities"

    # Volume basis improves MAE and RMSE
    assert vol_mae < mass_mae, f"Volume basis MAE ({vol_mae}) should outperform mass basis ({mass_mae})"
    assert vol_rmse < mass_rmse, f"Volume basis RMSE ({vol_rmse}) should outperform mass basis ({mass_rmse})"
