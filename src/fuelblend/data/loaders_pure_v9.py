"""Loader: Pure physical constants for candidate compounds (v9 Target 3).

Inputs:
  data/raw/pure_v9/vapour_pressure_draft.csv
  data/raw/pure_v9/lfl_draft.csv
  data/raw/pure_v9/viscosity_draft.csv
  data/raw/pure_v9/density_draft.csv

Rules:
* All numbers are parsed to float, strictly validated (finite, positive).
* LFL rows: must have lfl_source in {'measured', 'tabulated'}; Jones-rule estimates (0.55*Cst) are refused.
* Canonical units:
  - density: kg/m3
  - viscosity: mPa*s (dynamic) or cSt (kinematic)
  - vapour pressure: kPa
  - lfl: vol %
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RAW_DIR = ROOT / "data" / "raw" / "pure_v9"

VP_PATH = RAW_DIR / "vapour_pressure_draft.csv"
LFL_PATH = RAW_DIR / "lfl_draft.csv"
VISC_PATH = RAW_DIR / "viscosity_draft.csv"
DENS_PATH = RAW_DIR / "density_draft.csv"

ALLOWED_LFL_SOURCES = {"measured", "tabulated"}


def load_vapour_pressure_pure_v9(path: Path = VP_PATH) -> pd.DataFrame:
    """Load pure-component vapour pressure draft data."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    required = {"name", "smiles", "T_K", "P_kPa", "source", "ref", "page", "table", "confidence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in vapour pressure draft: {missing}")
    
    df["T_K"] = pd.to_numeric(df["T_K"], errors="raise")
    df["P_kPa"] = pd.to_numeric(df["P_kPa"], errors="raise")
    
    # Validation
    valid = (df["T_K"] > 0) & (df["P_kPa"] > 0) & np.isfinite(df["T_K"]) & np.isfinite(df["P_kPa"])
    return df[valid].copy()


def load_lfl_pure_v9(path: Path = LFL_PATH) -> pd.DataFrame:
    """Load pure-component LFL draft data, strictly refusing estimates."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    required = {"name", "smiles", "lfl_vol_pct", "T_C", "lfl_source", "ref", "page", "table", "confidence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in LFL draft: {missing}")

    df["lfl_vol_pct"] = pd.to_numeric(df["lfl_vol_pct"], errors="raise")
    df["T_C"] = pd.to_numeric(df["T_C"], errors="raise")

    # Filter out Jones-rule or any estimate source
    admitted = df["lfl_source"].str.lower().isin(ALLOWED_LFL_SOURCES)
    valid = (df["lfl_vol_pct"] > 0) & (df["lfl_vol_pct"] < 20.0) & admitted
    return df[valid].copy()


def load_viscosity_pure_v9(path: Path = VISC_PATH) -> pd.DataFrame:
    """Load pure-component viscosity draft data."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    required = {"name", "smiles", "T_K", "value_mPas", "visc_type", "ref", "page", "table", "confidence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in viscosity draft: {missing}")

    df["T_K"] = pd.to_numeric(df["T_K"], errors="raise")
    df["value_mPas"] = pd.to_numeric(df["value_mPas"], errors="raise")

    valid = (df["T_K"] > 0) & (df["value_mPas"] > 0) & np.isfinite(df["T_K"]) & np.isfinite(df["value_mPas"])
    return df[valid].copy()


def load_density_pure_v9(path: Path = DENS_PATH) -> pd.DataFrame:
    """Load pure-component density draft data."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    required = {"name", "smiles", "T_K", "value_kg_m3", "ref", "page", "table", "confidence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in density draft: {missing}")

    df["T_K"] = pd.to_numeric(df["T_K"], errors="raise")
    df["value_kg_m3"] = pd.to_numeric(df["value_kg_m3"], errors="raise")

    valid = (df["T_K"] > 0) & (df["value_kg_m3"] > 500) & (df["value_kg_m3"] < 1500)
    return df[valid].copy()
