#!/usr/bin/env python3
"""parse_boehm2022_cold_viscosity.py - Comprehensive multi-section parser for Boehm et al. (2022).

Extracts all 221 mixture compositions across Sections A (binary/pure), B (ternary 50/50 bases),
and C/D (real fuel surrogates) at 5 temperatures (-40 C, -20 C, -10 C, 0 C, +15 C).
Source: Boehm et al. Frontiers in Energy Research 2022, DOI: 10.3389/fenrg.2022.1074699.
Raw Excel: data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx

Outputs:
  - data/interim/external_staging/boehm2022_raw_extracted.csv

NOTE: Staged output ONLY. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import json
import pandas as pd
import numpy as np

RAW_EXCEL = Path("data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx")
STAGING_DIR = Path("data/interim/external_staging")

def clean_str(val):
    if pd.isna(val):
        return ""
    return str(val).strip()

def main():
    if not RAW_EXCEL.exists():
        print(f"ERROR: {RAW_EXCEL} not found.")
        return

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_excel(RAW_EXCEL, sheet_name="This Work", header=6)

    current_section = "A"
    current_components = []
    records = []

    for i in range(len(df)):
        r = df.iloc[i]
        c0 = clean_str(r.iloc[0])
        c1 = clean_str(r.iloc[1])
        c2 = clean_str(r.iloc[2])
        c3 = clean_str(r.iloc[3])
        c4 = clean_str(r.iloc[4])

        if c0 in ["A", "B", "C", "D"]:
            current_section = c0

        # Check if row is a component header
        try:
            float(c2)
            is_data = True
        except ValueError:
            is_data = False

        if not is_data:
            if c2 != "" and c2 != "nan":
                current_components = [c for c in [c2, c3, c4] if c != "" and c != "nan"]
            continue

        # Data row
        v2 = float(c2) if c2 != "" and c2 != "nan" else 0.0
        v3 = float(c3) if c3 != "" and c3 != "nan" else 0.0
        v4 = float(c4) if c4 != "" and c4 != "nan" else 0.0

        # Map component volume fractions
        comp_vol_pct = {}
        if current_section == "A":
            # Direct pure or binary
            if len(current_components) >= 1 and v2 > 0:
                comp_vol_pct[current_components[0]] = v2
            if len(current_components) >= 2 and v3 > 0:
                comp_vol_pct[current_components[1]] = v3
            if len(current_components) >= 3 and v4 > 0:
                comp_vol_pct[current_components[2]] = v4

        elif current_section == "B":
            # Ternary mixtures where one component is a 50/50 binary mixture
            # e.g., current_components[0] = 'ethylbenzene', current_components[1] = 'butylcyclohexane/heptylbenzene 50/50'
            for idx_c, c_name in enumerate(current_components):
                vol = [v2, v3, v4][idx_c]
                if vol <= 0.0:
                    continue
                if " 50/50" in c_name or "/" in c_name:
                    base_str = c_name.replace(" 50/50", "")
                    parts = [p.strip() for p in base_str.split("/")]
                    for p in parts:
                        comp_vol_pct[p] = comp_vol_pct.get(p, 0.0) + vol * 0.5
                else:
                    comp_vol_pct[c_name] = comp_vol_pct.get(c_name, 0.0) + vol

        else:
            # Section C & D: Surrogates and real fuels
            for idx_c, c_name in enumerate(current_components):
                vol = [v2, v3, v4][idx_c]
                if vol > 0:
                    comp_vol_pct[c_name] = vol

        # Temperatures and properties:
        # Col 5, 6: -40C (dens, visc)
        # Col 7, 8: -20C
        # Col 9, 10: -10C
        # Col 11, 12: 0C
        # Col 13, 14: 15C
        temps = [
            (-40.0, 233.15, r.iloc[5], r.iloc[6]),
            (-20.0, 253.15, r.iloc[7], r.iloc[8]),
            (-10.0, 263.15, r.iloc[9], r.iloc[10]),
            (0.0, 273.15, r.iloc[11], r.iloc[12]),
            (15.0, 288.15, r.iloc[13], r.iloc[14]),
        ]

        for t_c, t_k, dens_val, visc_val in temps:
            dens = float(dens_val) if pd.notna(dens_val) and str(dens_val).strip() != "" else np.nan
            visc = float(visc_val) if pd.notna(visc_val) and str(visc_val).strip() != "" else np.nan

            if pd.isna(dens) and pd.isna(visc):
                continue

            records.append({
                "excel_row": i,
                "section": current_section,
                "composition_json": json.dumps(comp_vol_pct),
                "n_comp": len(comp_vol_pct),
                "T_C": t_c,
                "T_K": t_k,
                "density_g_cc": dens,
                "viscosity_cSt": visc,
                "provenance": "Boehm et al. (2022) HEAT Lab"
            })

    out_df = pd.DataFrame(records)
    out_path = STAGING_DIR / "boehm2022_raw_extracted.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Extracted {len(out_df)} measurement condition points to {out_path}")
    print(f"Viscosity valid rows: {out_df['viscosity_cSt'].notna().sum()}")
    print(f"Density valid rows: {out_df['density_g_cc'].notna().sum()}")
    print(out_df.groupby(["section", "T_C"])["viscosity_cSt"].count().unstack(fill_value=0))

if __name__ == "__main__":
    main()
