from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CLEAN_CSV = ROOT / "data" / "processed" / "v2" / "thermoml_clean_hc.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "thermoml_outlier_filtering.md"


def test_clean_thermoml_exists():
    assert CLEAN_CSV.exists()
    assert DOC_PATH.exists()


def test_clean_thermoml_quality():
    df = pd.read_csv(CLEAN_CSV)
    assert len(df) > 40000

    # 1. No non-positive values
    assert (df["value"] > 0).all()

    # 2. No extreme temperatures
    assert (df["T_K"] >= 170.0).all()
    assert (df["T_K"] <= 550.0).all()

    # 3. No high pressures (> 500 kPa)
    p_valid = df["P_kPa"].dropna()
    assert (p_valid <= 500.0).all()

    # 4. Clean densities are liquid-phase
    dens = df[df["prop"] == "Mass density, kg/m3"]
    assert (dens["value"] >= 550.0).all()
    assert (dens["value"] <= 1100.0).all()

    # 5. Clean viscosities are non-negative and finite
    visc = df[df["prop"].str.startswith("Viscosity")]
    assert (visc["value"] >= 1e-5).all()
    assert (visc["value"] <= 2.0).all()
