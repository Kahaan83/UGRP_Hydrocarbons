"""Tests for fuelblend.data.loaders_pure_v9.

Verifies:
- All 4 pure draft files load correctly.
- LFL rows strictly require lfl_source in {'measured', 'tabulated'}; estimates like Jones rule are rejected.
- Canonical units: density in kg/m3, viscosity in mPa*s, vapor pressure in kPa, LFL in vol %.
"""
from pathlib import Path
import pandas as pd
import pytest

from fuelblend.data.loaders_pure_v9 import (
    load_vapour_pressure_pure_v9,
    load_lfl_pure_v9,
    load_viscosity_pure_v9,
    load_density_pure_v9,
)

ROOT = Path(__file__).resolve().parents[1]
PURE_DIR = ROOT / "data" / "raw" / "pure_v9"


def test_load_vapour_pressure_pure_v9():
    df = load_vapour_pressure_pure_v9()
    assert len(df) >= 20
    assert (df["T_K"] > 250).all() and (df["T_K"] < 600).all()
    assert (df["P_kPa"] > 0).all() and (df["P_kPa"] < 200).all()
    assert (df["confidence"] == "clear").all()


def test_load_lfl_pure_v9_refuses_estimates(tmp_path):
    # Test valid draft
    df = load_lfl_pure_v9()
    assert len(df) == 4
    assert (df["lfl_vol_pct"] > 0.4).all() and (df["lfl_vol_pct"] < 2.0).all()
    assert set(df["lfl_source"].str.lower()).issubset({"measured", "tabulated"})

    # Test that Jones-rule estimate row is refused
    bad_df = pd.DataFrame([{
        "name": "test_compound",
        "smiles": "CC",
        "lfl_vol_pct": 1.23,
        "T_C": 25.0,
        "lfl_source": "estimated: 0.55*Cst (Jones rule)",
        "ref": "test",
        "page": "p1",
        "table": "t1",
        "confidence": "clear"
    }])
    p = tmp_path / "bad_lfl.csv"
    bad_df.to_csv(p, index=False)
    loaded = load_lfl_pure_v9(p)
    assert len(loaded) == 0, "Jones-rule estimate was admitted but should be refused!"


def test_load_viscosity_pure_v9():
    df = load_viscosity_pure_v9()
    assert len(df) >= 10
    assert (df["value_mPas"] > 0.5).all() and (df["value_mPas"] < 20.0).all()
    assert "1-methylnaphthalene" in df["name"].values
    assert "perhydrofluorene" in df["name"].values


def test_load_density_pure_v9():
    df = load_density_pure_v9()
    assert len(df) >= 10
    assert (df["value_kg_m3"] > 800).all() and (df["value_kg_m3"] < 1100).all()
    assert "1-methylnaphthalene" in df["name"].values
    assert "2,6-dimethylnaphthalene" in df["name"].values
