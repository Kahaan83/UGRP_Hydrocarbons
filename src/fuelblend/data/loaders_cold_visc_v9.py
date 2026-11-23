"""Loader for v9 real-fuel cold viscosity (Boehm et al. 2022, sheet `This Work`, verified against the workbook).

Reads data/raw/cold_visc_v9/boehm2022_cold_visc.csv.  Only NEAT fuels (one component, kind real_fuel) are
returned: blends of a real fuel with pure compounds have no InChIKey composition and are dropped here.
The Edwards 2017 draft is NOT read (source not retrievable, see results/v9/verification.md).

system_id here is the provisional ``fuel:boehm2022:<sample>``; scripts/v9/build_unified_v9.py rewrites it to the
canonical fuel id through the identity map.  kinematic viscosity stays cSt, density is kg/m3 (g/cc x 1000 in the
raw file), T_K = T_C + 273.15.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "data" / "raw" / "cold_visc_v9" / "boehm2022_cold_visc.csv"


def load_cold_visc_v9(path: Path = RAW) -> pd.DataFrame:
    raw = pd.read_csv(path)
    rows = []
    for _, r in raw.iterrows():
        if r["kind"] != "real_fuel" or len(json.loads(r["components"])) != 1:
            continue
        t_k = round(float(r["T_C"]) + 273.15, 2)
        qual = r.get("qualifier")
        qual = "" if pd.isna(qual) else str(qual)
        val = float(r["value"])
        rows.append({
            "property": r["property"], "system_id": f"fuel:boehm2022:{r['sample_id']}",
            "sample_id": r["sample_id"], "x_key": "1.0", "T_key": f"{t_k:.1f}", "P_key": "101",
            "T_K": t_k, "P_kPa": 101.325, "kind": "real_fuel", "n_comp": 1, "is_hc": True,
            "comp_inchikeys": json.dumps([]), "x_mole": json.dumps([]), "w_mass": json.dumps([1.0]),
            "fuel_id": r["sample_id"], "ref": f"{r['source_id']}_{r['doi']}",
            "composition_level": "fuel_id_only", "comp_coverage": 1.0, "comp_basis": "mass",
            "value": val, "unit": r["unit"], "source_id": r["source_id"], "qualifier": qual,
            "confidence": r["confidence"], "table_or_page": r["table_or_page"],
        })
    return pd.DataFrame(rows)
