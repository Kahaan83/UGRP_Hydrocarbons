"""Loaders for v10 pure physical constants data.

Targeted properties:
- Vapor pressure (converts mmHg, bar, Pa to kPa)
- LFL (measured only; rejects Jones, homologue, or any estimates)
- Viscosity and Density
- Scanned rows (extraction=image and confidence=clear forced to unsure)
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW_DIR = REPO / "data" / "raw" / "pure_v10"

REJECTED_LFL_KEYWORDS = ("estimate", "jones", "homologue")


def _clean_numeric(val: Any) -> float:
    """Normalize minus signs and parse float."""
    if pd.isna(val):
        return np.nan
    s = str(val).strip()
    s = s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
    return float(s)


def _convert_pressure(val: float, unit: str) -> tuple[float, str]:
    """Convert pressure to kPa."""
    u = unit.strip().lower()
    if u in ("mmhg", "torr"):
        return val * (101.325 / 760.0), "kPa"
    elif u == "bar":
        return val * 100.0, "kPa"
    elif u in ("pa", "pascal", "n/m2", "n/m^2"):
        return val / 1000.0, "kPa"
    elif u == "kpa":
        return val, "kPa"
    return val, unit


def load_pure_constants_v10(path: Path | None = None) -> pd.DataFrame:
    """Load pure physical constants data table."""
    if path is not None:
        p = Path(path)
        if not p.exists():
            return pd.DataFrame()
        df = pd.read_csv(p)
    else:
        files = sorted(RAW_DIR.glob("*_draft.csv"))
        dfs = []
        for f in files:
            if f.exists() and f.stat().st_size > 0:
                try:
                    d = pd.read_csv(f)
                    if len(d) > 0:
                        dfs.append(d)
                except Exception:
                    pass
        if not dfs:
            return pd.DataFrame()
        df = pd.concat(dfs, ignore_index=True)

    if len(df) == 0:
        return pd.DataFrame()

    rows = []
    for _, r in df.iterrows():
        prop = str(r.get("property", "")).strip().lower()
        method = str(r.get("method", "")).strip().lower()

        # Rule: LFL must be measured; reject estimates, Jones rule, or homologue extrapolations
        if "lfl" in prop or prop == "lower_flammability_limit":
            if any(kw in method for kw in REJECTED_LFL_KEYWORDS):
                continue

        val = _clean_numeric(r["value"])
        unit = str(r.get("unit", "")).strip()

        # Pressure unit conversion
        if "vapor_pressure" in prop or prop in ("vp", "vapour_pressure"):
            val, unit = _convert_pressure(val, unit)

        # Scanned rows forced to unsure
        extraction = str(r.get("extraction", "text_layer")).strip()
        confidence = str(r.get("confidence", "clear")).strip()
        if extraction == "image" and confidence == "clear":
            confidence = "unsure"

        # Temperatures
        t_c_raw = r.get("T_C")
        t_c = None
        t_k = None
        if pd.notna(t_c_raw) and str(t_c_raw).strip() != "":
            t_c = _clean_numeric(t_c_raw)
            t_k = round(t_c + 273.15, 3)
        elif "T_K" in r and pd.notna(r["T_K"]) and str(r["T_K"]).strip() != "":
            t_k = _clean_numeric(r["T_K"])
            t_c = round(t_k - 273.15, 3)

        compound = str(r.get("compound", r.get("name", ""))).strip()
        sample_id = str(r.get("sample_id", compound)).strip()

        rows.append({
            "compound": compound,
            "name": compound,
            "inchikey": str(r.get("inchikey", "")).strip(),
            "source_id": str(r.get("source_id", "v10")).strip(),
            "table_or_page": str(r.get("table_or_page", "")).strip(),
            "sample_id": sample_id,
            "kind": str(r.get("kind", "pure")).strip(),
            "property": str(r.get("property", "")).strip(),
            "value": val,
            "unit": unit,
            "T_C": t_c,
            "T_K": t_k,
            "P_kPa": val if "vapor_pressure" in prop or prop in ("vp", "vapour_pressure") else np.nan,
            "method": str(r.get("method", "")).strip(),
            "qualifier": str(r.get("qualifier", "")).strip(),
            "note": str(r.get("note", "")).strip(),
            "page": r.get("page"),
            "table": str(r.get("table", "")).strip(),
            "confidence": confidence,
            "raw_text": str(r.get("raw_text", "")),
            "extraction": extraction,
        })

    return pd.DataFrame(rows)


def load_vapour_pressure_v10(path: Path | None = None) -> pd.DataFrame:
    """Load pure-component vapour pressure data, converting units to kPa."""
    default_p = RAW_DIR / "vapour_pressure_draft.csv"
    target_path = Path(path) if path is not None else default_p
    df = load_pure_constants_v10(target_path)
    if len(df) == 0:
        return df
    vp_mask = df["property"].str.lower().str.contains("vapor_pressure|vapour_pressure|vp")
    return df[vp_mask].copy()


def load_lfl_pure_v10(path: Path | None = None) -> pd.DataFrame:
    """Load pure-component LFL data, refusing estimates."""
    default_p = RAW_DIR / "lfl_draft.csv"
    target_path = Path(path) if path is not None else default_p
    df = load_pure_constants_v10(target_path)
    if len(df) == 0:
        return df
    lfl_mask = df["property"].str.lower().str.contains("lfl|lower_flammability_limit")
    return df[lfl_mask].copy()


def load_viscosity_pure_v10(path: Path | None = None) -> pd.DataFrame:
    """Load pure-component viscosity data."""
    default_p = RAW_DIR / "viscosity_draft.csv"
    target_path = Path(path) if path is not None else default_p
    df = load_pure_constants_v10(target_path)
    if len(df) == 0:
        return df
    visc_mask = df["property"].str.lower().str.contains("visc")
    return df[visc_mask].copy()


def load_density_pure_v10(path: Path | None = None) -> pd.DataFrame:
    """Load pure-component density data."""
    default_p = RAW_DIR / "density_draft.csv"
    target_path = Path(path) if path is not None else default_p
    df = load_pure_constants_v10(target_path)
    if len(df) == 0:
        return df
    dens_mask = df["property"].str.lower().str.contains("dens")
    return df[dens_mask].copy()
