#!/usr/bin/env python
"""
scripts/v2/assemble_distillation_real_fuels.py

Assembles and standardizes real-fuel distillation data from:
1. CRC Report 647: 57 aviation fuels (ASTM D86 test 1, test 2, and ASTM D2887)
2. CRC Report AV-27-18: 19 aviation fuels (ASTM D2887)
3. Fortin et al. 2012 (Energy & Fuels): 14 fuels (11 RP-1 blends, 1 RP-1 4572, 2 RP-2 fuels, ADC method)
4. AJFTD: 99 alternative jet fuels (ASTM D86)

Standardizes recovery points to:
T0_C (IBP), T5_C, T10_C, T20_C, T30_C, T40_C, T50_C, T60_C, T70_C, T80_C, T90_C, T95_C, FBP_C (T100).

Outputs: data/processed/v2/distillation_real_fuels.csv
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
OUT_DIR = BASE / "data" / "processed" / "v2"
OUT_CSV = OUT_DIR / "distillation_real_fuels.csv"

CUTS = ["T0", "T5", "T10", "T20", "T30", "T40", "T50", "T60", "T70", "T80", "T90", "T95", "T100"]


def load_crc647() -> list[dict]:
    p = BASE / "data" / "raw" / "pour_point" / "crc647_fuel_properties.csv"
    if not p.exists():
        print(f"Warning: {p} not found")
        return []
    df = pd.read_csv(p)
    d = df[df["property"].str.startswith("distillation_")].copy()
    
    # Identify unique combinations of fuel_code and method
    rows = []
    groups = d.groupby(["fuel_code", "fuel_name", "method"])
    for (f_code, f_name, method), grp in groups:
        rec = {
            "fuel_id": f"crc647_{f_code}",
            "source": "crc647",
            "fuel_name": f_name,
            "fuel_class": "Aviation Kerosene" if "Jet" in f_name or "JP" in f_name else "Aviation Fuel",
            "method": method,
        }
        for c in CUTS:
            val = grp[grp["property"] == f"distillation_{c}_C"]["value_std"].values
            col_name = "FBP_C" if c == "T100" else f"{c}_C"
            rec[col_name] = float(val[0]) if len(val) > 0 and pd.notna(val[0]) else np.nan
        rows.append(rec)
    return rows


def load_av2718() -> list[dict]:
    p = BASE / "data" / "processed" / "crc_fuels.csv"
    if not p.exists():
        print(f"Warning: {p} not found")
        return []
    df = pd.read_csv(p)
    av = df[(df["source"] == "av2718") & df["property"].str.startswith("distillation_")].copy()
    
    rows = []
    groups = av.groupby(["fuel_id", "fuel_name", "method"])
    for (f_id, f_name, method), grp in groups:
        rec = {
            "fuel_id": f"av2718_{f_id}",
            "source": "av2718",
            "fuel_name": f_name,
            "fuel_class": "Jet-A / Synthetic Kerosene",
            "method": method,
        }
        for c in CUTS:
            val = grp[grp["property"] == f"distillation_{c}_C"]["value"].values
            col_name = "FBP_C" if c == "T100" else f"{c}_C"
            rec[col_name] = float(val[0]) if len(val) > 0 and pd.notna(val[0]) else np.nan
        rows.append(rec)
    return rows


def load_fortin() -> list[dict]:
    p = BASE / "data" / "processed" / "fortin2012_rp1_blends.csv"
    if not p.exists():
        print(f"Warning: {p} not found")
        return []
    df = pd.read_csv(p)
    
    rows = []
    groups = df.groupby(["fuel_id", "fuel_class"])
    for (f_id, f_class), grp in groups:
        rec = {
            "fuel_id": f"fortin2012_{f_id.replace(' ', '_')}",
            "source": "fortin2012",
            "fuel_name": f"Fortin et al. 2012 - {f_id}",
            "fuel_class": f_class,
            "method": "Advanced Distillation Curve (ADC)",
        }
        # T0 from T0.025 or initial_boiling_sustained
        t0 = grp[grp["property"] == "distillation_T0.025_C"]["value_std"].values
        if len(t0) == 0:
            t0 = grp[grp["property"] == "initial_boiling_sustained_C"]["value_std"].values
        rec["T0_C"] = float(t0[0]) if len(t0) > 0 and pd.notna(t0[0]) else np.nan
        
        for c in ["T5", "T10", "T20", "T30", "T40", "T50", "T60", "T70", "T80", "T90"]:
            val = grp[grp["property"] == f"distillation_{c}_C"]["value_std"].values
            rec[f"{c}_C"] = float(val[0]) if len(val) > 0 and pd.notna(val[0]) else np.nan
            
        rec["T95_C"] = np.nan
        rec["FBP_C"] = np.nan
        rows.append(rec)
    return rows


def load_ajftd() -> list[dict]:
    p = BASE / "data" / "processed" / "distillation_ajftd.csv"
    if not p.exists():
        print(f"Warning: {p} not found")
        return []
    df = pd.read_csv(p)
    
    rows = []
    for r in df.itertuples():
        rec = {
            "fuel_id": f"ajftd_{r.fuel_id}",
            "source": "ajftd",
            "fuel_name": str(r.fuel_label),
            "fuel_class": "Alternative Jet Fuel",
            "method": str(r.method) if pd.notna(r.method) else "ASTM D86",
        }
        for c in ["T0", "T5", "T10", "T20", "T30", "T40", "T50", "T60", "T70", "T80", "T90", "T95"]:
            rec[f"{c}_C"] = getattr(r, f"{c}_C", np.nan)
        rec["FBP_C"] = getattr(r, "T100_C", np.nan)
        rows.append(rec)
    return rows


def main():
    print("Assembling real fuel distillation datasets...")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    crc_rows = load_crc647()
    print(f"  Loaded {len(crc_rows)} curve records from CRC 647")
    
    av_rows = load_av2718()
    print(f"  Loaded {len(av_rows)} curve records from CRC AV-27-18")
    
    fortin_rows = load_fortin()
    print(f"  Loaded {len(fortin_rows)} curve records from Fortin et al. 2012")
    
    ajftd_rows = load_ajftd()
    print(f"  Loaded {len(ajftd_rows)} curve records from AJFTD")
    
    all_rows = crc_rows + av_rows + fortin_rows + ajftd_rows
    out_df = pd.DataFrame(all_rows)
    
    # Sort and reorder columns
    cols = [
        "fuel_id",
        "source",
        "fuel_name",
        "fuel_class",
        "method",
        "T0_C",
        "T5_C",
        "T10_C",
        "T20_C",
        "T30_C",
        "T40_C",
        "T50_C",
        "T60_C",
        "T70_C",
        "T80_C",
        "T90_C",
        "T95_C",
        "FBP_C",
    ]
    out_df = out_df[cols]
    out_df.to_csv(OUT_CSV, index=False)
    print(f"\nSuccessfully assembled {len(out_df)} standardized distillation curves into {OUT_CSV}")
    print("\nBreakdown by source:")
    print(out_df["source"].value_counts())
    print("\nBreakdown by method:")
    print(out_df["method"].value_counts())


if __name__ == "__main__":
    main()
