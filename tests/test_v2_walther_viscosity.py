from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
EVAL_PATH = ROOT / "results" / "v2" / "improvement" / "viscosity_walther_eval.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "viscosity_walther_upgrade.md"




def test_walther_calculation():
    from fuelblend.thermo.walther_viscosity import fit_walther_parameters, walther_viscosity, extrapolate_cryogenic_walther
    # Typical Jet A: nu(+20 C) = 1.6 cSt, nu(-20 C) = 4.2 cSt
    nu_m40 = extrapolate_cryogenic_walther(-40.0, 20.0, 1.60, -20.0, 4.20)
    # At -40 C, viscosity should be in realistic range ~10-15 cSt
    assert 8.0 <= nu_m40 <= 20.0


def test_walther_eval_content():
    df = pd.read_csv(EVAL_PATH).set_index("model")
    assert "Physics Baseline (Raw Refutas)" in df.index
    assert "ASTM D341 Walther Viscosity Curvature Model" in df.index
