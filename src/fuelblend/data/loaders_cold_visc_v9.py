"""Loader for v9 real-fuel cold viscosity data (Edwards 2017 and Boehm 2022).

Reads data/raw/cold_visc_v9/:
- edwards2017_cold_visc_draft.csv (or edwards2017_cold_visc.csv)
- boehm2022_cold_visc_draft.csv (or boehm2022_cold_visc.csv)
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW_DIR = REPO / "data" / "raw" / "cold_visc_v9"

CANONICAL_SYSTEM_IDS = {
    "POSF-10325": "fuel:ajftd:2",
    "A-2 (POSF 10325)": "fuel:ajftd:2",
    "POSF-10289": "fuel:ajftd:39",
    "A-3 (POSF 10289)": "fuel:ajftd:39",
    "A-3 (POSF10289)": "fuel:ajftd:39",
    "POSF-4658": "fuel:ajftd:23",
    "POSF-10264": "fuel:afrl:posf10264",
    "A-1 (POSF 10264)": "fuel:afrl:posf10264",
    "C-1 (POSF 13718)": "fuel:ajftd:24570",
    "C-7 (POSF 13727)": "fuel:ajftd:24567",
    "C-9 (POSF 12968)": "fuel:ajftd:143",
    "HRJ Camelina (POSF 7720)": "fuel:fuellib:posf7720",
    "C-3 (POSF 12959)": "fuel:afrl:posf12959",
    "C-8 (POSF 13491)": "fuel:afrl:posf13491",
    "POSF 12945": "fuel:afrl:posf12945",
    "POSF12945": "fuel:afrl:posf12945",
}


def _find_file(base_name: str) -> Path:
    final = RAW_DIR / f"{base_name}.csv"
    if final.exists():
        return final
    draft = RAW_DIR / f"{base_name}_draft.csv"
    if draft.exists():
        return draft
    raise FileNotFoundError(f"Neither {final} nor {draft} found in {RAW_DIR}")


def load_cold_visc_v9() -> pd.DataFrame:
    """Load v9 real-fuel cold viscosity measurements."""
    edwards_file = _find_file("edwards2017_cold_visc")
    boehm_file = _find_file("boehm2022_cold_visc")

    df_ed = pd.read_csv(edwards_file)
    df_bm = pd.read_csv(boehm_file)
    raw = pd.concat([df_ed, df_bm], ignore_index=True)

    rows = []
    for _, r in raw.iterrows():
        sid = str(r["sample_id"])
        canonical_sys = CANONICAL_SYSTEM_IDS.get(sid, f"fuel:{r['source_id']}:{sid}")
        prop = r["property"]
        val = float(r["value"])
        unit = str(r.get("unit", ""))

        t_c = float(r["T_C"])
        t_k = round(t_c + 273.15, 2)

        qual = str(r.get("qualifier", ""))
        if pd.isna(qual):
            qual = ""

        # Dynamic viscosity conversion check
        if prop == "viscosity_dyn_mPas" and "density" not in str(r.get("note", "")):
            qual = "dynamic_no_density"

        rows.append({
            "property": prop,
            "system_id": canonical_sys,
            "sample_id": sid,
            "x_key": "1.0",
            "T_key": f"{t_k:.1f}",
            "P_key": "101.3",
            "T_K": t_k,
            "P_kPa": 101.325,
            "kind": r.get("kind", "real_fuel"),
            "n_comp": 1 if r.get("kind") == "real_fuel" else 2,
            "is_hc": True,
            "comp_inchikeys": json.dumps([]),
            "x_mole": json.dumps([]),
            "w_mass": r.get("fractions", "[1.0]"),
            "fuel_id": sid,
            "ref": f"{r['source_id']}_{r['doi']}",
            "composition_level": "fuel_id_only",
            "comp_coverage": 1.0,
            "comp_basis": r.get("fraction_basis", "mass"),
            "value": val,
            "unit": unit,
            "vmin": val,
            "vmax": val,
            "q25": val,
            "q75": val,
            "n_obs": 1,
            "n_refs": 1,
            "n_sources": 1,
            "sources": r["source_id"],
            "n_outlier": 0,
            "conflict": False,
            "source_id": r["source_id"],
            "qualifier": qual,
            "confidence": r.get("confidence", "clear"),
            "table_or_page": r.get("table_or_page", ""),
        })

    return pd.DataFrame(rows)
