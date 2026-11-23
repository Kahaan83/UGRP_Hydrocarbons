"""Loaders for v10 real-fuel properties and composition data.

Supports:
- AFRL-TR-2020-0017 (afrl2020)
- AIAA 2017-0146 (edwards2017)
- Minus-sign variants (Unicode minus U+2212, en-dash U+2013, em-dash U+2014, ASCII hyphen)
- Unit normalization (cSt, mm2/s -> cSt, degF -> C)
- Composition closure checking (flags composition_open if |sum - 100| >= 0.5)
- Hydrocarbon type tables without carbon numbers (class_only)
- Mapping sample_id through canonical_id() to AJFTD fuel ids.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW_DIR = REPO / "data" / "raw" / "fuels_v10"


def _clean_numeric(val: Any) -> float:
    """Normalize minus signs and parse float."""
    if pd.isna(val):
        return np.nan
    s = str(val).strip()
    s = s.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
    return float(s)


def _get_canonical_id(source_id: str, sample_id: str) -> str:
    """Resolve sample_id via canonical_id in scripts.v9.build_unified_v9."""
    try:
        import sys
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from scripts.v9.build_unified_v9 import canonical_id
        return canonical_id(source_id, sample_id)
    except Exception:
        # Fallback if label table is unavailable
        d = re.search(r"(\d{4,5})", str(sample_id))
        if d and "POSF" in str(sample_id).upper():
            return f"fuel:posf:{d.group(1)}"
        return f"fuel:v10:{source_id}:{sample_id}"


def load_fuel_properties_v10(path: Path | None = None) -> pd.DataFrame:
    """Load v10 real-fuel properties in measurements.csv schema.

    If path is None, scans data/raw/fuels_v10/ for verified or draft properties CSVs.
    """
    if path is not None:
        files = [Path(path)]
    else:
        files = sorted(RAW_DIR.glob("*_properties*.csv"))

    if not files:
        return pd.DataFrame()

    dfs = []
    for f in files:
        if f.exists() and f.stat().st_size > 0:
            try:
                df = pd.read_csv(f)
                if len(df) > 0:
                    dfs.append(df)
            except Exception:
                pass

    if not dfs:
        return pd.DataFrame()

    raw = pd.concat(dfs, ignore_index=True)
    if len(raw) == 0:
        return pd.DataFrame()

    rows = []
    for _, r in raw.iterrows():
        sid = str(r["sample_id"]).strip()
        source_id = str(r.get("source_id", "v10")).strip()
        canonical_sys = _get_canonical_id(source_id, sid)
        prop = str(r["property"]).strip()
        val = _clean_numeric(r["value"])

        # Handle units
        unit = str(r.get("unit", "")).strip()
        if unit in ("degF", "°F", "F"):
            val = (val - 32.0) * 5.0 / 9.0
            unit = "C"
        elif unit in ("mm2/s", "mm^2/s"):
            unit = "cSt"

        # Temperature
        t_c_raw = r.get("T_C")
        t_c = None
        t_k = None
        if pd.notna(t_c_raw) and str(t_c_raw).strip() != "":
            t_c = _clean_numeric(t_c_raw)
            t_k = round(t_c + 273.15, 2)

        qual = r.get("qualifier")
        qual = "" if pd.isna(qual) else str(qual)

        rows.append({
            "property": prop,
            "system_id": canonical_sys,
            "sample_id": sid,
            "x_key": "1.0",
            "T_key": f"{t_k:.1f}" if t_k is not None else "amb",
            "P_key": "101",
            "T_C": t_c,
            "T_K": t_k,
            "P_kPa": 101.325,
            "kind": r.get("kind", "real_fuel"),
            "n_comp": 1,
            "is_hc": True,
            "comp_inchikeys": json.dumps([]),
            "x_mole": json.dumps([]),
            "w_mass": json.dumps([1.0]),
            "fuel_id": sid,
            "ref": f"{source_id}_{r.get('table_or_page', '')}",
            "composition_level": "detailed_gcxgc",
            "comp_coverage": 1.0,
            "comp_basis": "mass",
            "value": val,
            "unit": unit,
            "source_id": source_id,
            "qualifier": qual,
            "confidence": r.get("confidence", "clear"),
            "table_or_page": r.get("table_or_page", ""),
            "raw_text": r.get("raw_text", ""),
            "extraction": r.get("extraction", "text_layer"),
        })

    return pd.DataFrame(rows)


def load_fuel_composition_v10(path: Path | None = None) -> pd.DataFrame:
    """Load v10 real-fuel compositions in AJFTD GCxGC long layout.

    If path is None, scans data/raw/fuels_v10/ for verified or draft composition CSVs.
    """
    if path is not None:
        files = [Path(path)]
    else:
        files = sorted(RAW_DIR.glob("*_composition*.csv"))

    if not files:
        return pd.DataFrame()

    dfs = []
    for f in files:
        if f.exists() and f.stat().st_size > 0:
            try:
                df = pd.read_csv(f)
                if len(df) > 0:
                    dfs.append(df)
            except Exception:
                pass

    if not dfs:
        return pd.DataFrame()

    raw = pd.concat(dfs, ignore_index=True)
    if len(raw) == 0:
        return pd.DataFrame()

    parsed_rows = []
    for _, r in raw.iterrows():
        sid = str(r.get("sample_id", r.get("fuel_id", ""))).strip()
        source_id = str(r.get("source_id", "v10")).strip()
        canonical_sys = _get_canonical_id(source_id, sid)

        nc_raw = r.get("n_carbon")
        is_no_carbon = pd.isna(nc_raw) or str(nc_raw).strip() in ("", "0", "nan")
        n_carbon = 0 if is_no_carbon else int(_clean_numeric(nc_raw))

        val = _clean_numeric(r["value"])

        parsed_rows.append({
            "fuel_id": canonical_sys,
            "sample_id": sid,
            "fuel_label": str(r.get("fuel_label", sid)).strip(),
            "hc_class": str(r.get("class", r.get("hc_class", ""))).strip(),
            "n_carbon": n_carbon,
            "cn_open_ended": False,
            "metric_raw": "mass_pct",
            "provenance": source_id,
            "mass_pct": val,
            "source_id": source_id,
            "no_carbon_split": is_no_carbon,
            "confidence": r.get("confidence", "clear"),
            "raw_text": r.get("raw_text", ""),
            "extraction": r.get("extraction", "text_layer"),
        })

    df_out = pd.DataFrame(parsed_rows)

    # Compute per-fuel flags: composition_open and class_only
    # composition_open: True if |sum(mass_pct) - 100| >= 0.5
    fuel_totals = df_out.groupby("sample_id")["mass_pct"].sum()
    open_map = {sid: bool(abs(tot - 100.0) >= 0.5) for sid, tot in fuel_totals.items()}
    df_out["composition_open"] = df_out["sample_id"].map(open_map)

    # class_only: True if all rows for this fuel have no carbon split
    class_only_map = df_out.groupby("sample_id")["no_carbon_split"].all().to_dict()
    df_out["class_only"] = df_out["sample_id"].map(class_only_map)

    return df_out
