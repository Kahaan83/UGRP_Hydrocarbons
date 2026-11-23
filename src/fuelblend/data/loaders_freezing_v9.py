"""Loader for v9 real-fuel freezing and composition data (Edwards 2017 AFRL reference fuels).

Reads data/raw/freezing_v9/:
- edwards2017_properties_draft.csv (or edwards2017_properties.csv)
- edwards2017_composition_draft.csv (or edwards2017_composition.csv)

Returns (properties_df, compositions_df) where:
- properties_df has the measurements_dedup schema with canonical system_ids
- compositions_df has the ajftd_gcxgc_composition.csv schema with no_carbon_split flag
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW_DIR = REPO / "data" / "raw" / "freezing_v9"

# Identity map for Edwards 2017 reference fuels
CANONICAL_SYSTEM_IDS = {
    "POSF-10325": "fuel:ajftd:2",
    "POSF-10289": "fuel:ajftd:39",
    "POSF-4658": "fuel:ajftd:23",
    "POSF-10264": "fuel:afrl:posf10264",
}


def _find_file(base_name: str) -> Path:
    final = RAW_DIR / f"{base_name}.csv"
    if final.exists():
        return final
    draft = RAW_DIR / f"{base_name}_draft.csv"
    if draft.exists():
        return draft
    raise FileNotFoundError(f"Neither {final} nor {draft} found in {RAW_DIR}")


def load_freezing_v9() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load v9 real-fuel freezing properties and compositions."""
    prop_path = _find_file("edwards2017_properties")
    comp_path = _find_file("edwards2017_composition")

    raw_props = pd.read_csv(prop_path)
    raw_comps = pd.read_csv(comp_path)

    # 1. Build properties in measurements_dedup schema
    prop_rows = []
    for _, r in raw_props.iterrows():
        sid = r["sample_id"]
        canonical_sys = CANONICAL_SYSTEM_IDS.get(sid, f"fuel:edwards2017:{sid}")
        prop = r["property"]
        val = float(r["value"])

        # Temperature handling
        t_c = r.get("T_C")
        t_k = float(t_c) + 273.15 if pd.notna(t_c) and t_c != "" else None

        prop_rows.append({
            "property": prop,
            "system_id": canonical_sys,
            "sample_id": sid,
            "x_key": "1.0",
            "T_key": f"{t_k:.1f}" if t_k is not None else "amb",
            "P_key": "101.3",
            "T_K": t_k,
            "P_kPa": 101.325,
            "kind": "real_fuel",
            "n_comp": 1,
            "is_hc": True,
            "comp_inchikeys": json.dumps([]),
            "x_mole": json.dumps([]),
            "w_mass": json.dumps([1.0]),
            "fuel_id": sid,
            "ref": f"{r['source_id']}_{r['doi']}",
            "composition_level": "detailed_gcxgc",
            "comp_coverage": 1.0,
            "comp_basis": "mass",
            "value": val,
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
            "qualifier": r.get("qualifier", ""),
            "confidence": r.get("confidence", "clear"),
            "table_or_page": r.get("table_or_page", ""),
        })
    props_df = pd.DataFrame(prop_rows)

    # 2. Build compositions in ajftd_gcxgc_composition schema
    comp_rows = []
    for _, r in raw_comps.iterrows():
        sid = r["fuel_id"]
        canonical_sys = CANONICAL_SYSTEM_IDS.get(sid, f"fuel:edwards2017:{sid}")
        nc = r["n_carbon"]
        no_split = pd.isna(nc) or nc == 0 or nc == ""

        comp_rows.append({
            "fuel_id": canonical_sys,
            "sample_id": sid,
            "fuel_label": r["fuel_label"],
            "hc_class": r["class"],
            "n_carbon": int(nc) if not no_split else 0,
            "cn_open_ended": False,
            "metric_raw": "mass_pct",
            "provenance": r["source_id"],
            "mass_pct": float(r["value"]),
            "source_id": r["source_id"],
            "no_carbon_split": no_split,
            "confidence": r.get("confidence", "clear"),
        })
    comps_df = pd.DataFrame(comp_rows)

    return props_df, comps_df
