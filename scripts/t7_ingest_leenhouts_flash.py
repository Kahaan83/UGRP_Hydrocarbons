#!/usr/bin/env python
"""T7: ingest the Leenhouts flash-point MIXTURE rows.

    python scripts/t7_ingest_leenhouts_flash.py

Parses flashpoint_mix_exp.csv into canonical schema:
- comp_inchikeys: JSON list of sorted InChIKeys
- x_mole: JSON list of normalized mole fractions aligned to comp_inchikeys
- w_mass: JSON list of mass fractions computed via mole_to_mass
- is_hc: boolean indicating if all components are hydrocarbons (C, H only)
- system_id: md5 hash of '|'.join(comp_inchikeys)[:12]
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
from common import is_hydrocarbon, mole_to_mass  # noqa: E402

SRC = BASE / "data/interim/leenhouts_benchmark/datasets/datasets/flashpoint/processed_data/flashpoint_mix_exp.csv"
OUT = BASE / "data/processed/t7_leenhouts_flash_mix.csv"
REPORT = BASE / "reports/T7_LEENHOUTS_FLASH.md"


def main() -> None:
    print(f"Reading {SRC}...")
    df = pd.read_csv(SRC)
    print(f"Raw mix rows: {len(df)}")

    cols = ["component_inchi_0", "component_inchi_1", "component_inchi_2", "component_inchi_3"]
    fracs = ["fraction_0", "fraction_1", "fraction_2", "fraction_3"]

    out = []
    for r in df.itertuples():
        pairs = []
        for c, f in zip(cols, fracs):
            v = getattr(r, c)
            fv = getattr(r, f)
            if pd.notna(v) and pd.notna(fv) and float(fv) > 0:
                inchi_str = str(v).strip()
                ik = Chem.InchiToInchiKey(inchi_str)
                mol = Chem.MolFromInchi(inchi_str)
                mw = Descriptors.MolWt(mol) if mol else 100.0
                is_c_h = is_hydrocarbon(mol) if mol else False
                pairs.append({"ik": ik, "x": float(fv), "mw": mw, "is_hc": is_c_h})
        if not pairs:
            continue
        pairs.sort(key=lambda t: t["ik"])
        keys = [p["ik"] for p in pairs]
        x_raw = [p["x"] for p in pairs]
        x_tot = sum(x_raw)
        x = [v / x_tot for v in x_raw]
        mws = [p["mw"] for p in pairs]
        w = mole_to_mass(x, mws)
        all_hc = all(p["is_hc"] for p in pairs)
        sys_id = hashlib.md5("|".join(keys).encode("ascii")).hexdigest()[:12]
        out.append({
            "property": "flash_point_C",
            "fuel_id": f"leenhouts_mix_{hashlib.md5(str(float(r.value)).encode()).hexdigest()[:6]}",
            "system_id": sys_id,
            "kind": "binary" if len(keys) == 2 else "multi",
            "comp_inchikeys": json.dumps(keys),
            "x_mole": json.dumps([round(float(v), 6) for v in x]),
            "w_mass": json.dumps([round(float(v), 6) for v in w]),
            "value": float(r.value),
            "T_K": 298.15,
            "P_kPa": 101.325,
            "n_comp": len(keys),
            "is_hc": all_hc,
            "source": "leenhouts",
            "ref": "Leenhouts et al. benchmark (Zenodo 19914760); mixture flash (T7 ingestion)",
            "raw_value": float(r.value),
            "raw_unit": "Celcius",
            "quality_weight": 1.0,
            "tier": "L3",
            "method": "Le Chatelier baseline; mixture (2-4 components, molar fractions)"
        })

    od = pd.DataFrame(out)
    od.to_csv(OUT, index=False)
    print(f"Saved {OUT}. Total rows: {len(od)}")
    print(f"is_hc distribution:\n{od['is_hc'].value_counts().to_dict()}")
    print(f"Sample comp_inchikeys: {od['comp_inchikeys'].iloc[0]}")
    print(f"Sample system_id: {od['system_id'].iloc[0]}")


if __name__ == "__main__":
    main()
