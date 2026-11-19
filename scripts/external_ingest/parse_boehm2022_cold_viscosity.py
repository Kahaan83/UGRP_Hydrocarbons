#!/usr/bin/env python3
"""parse_boehm2022_cold_viscosity.py - Parse cold kinematic viscosity and density from Boehm et al. (2022).

Source: Boehm et al. "Viscosity and Density Measurements of Pure Hydrocarbons and Binary/Ternary Blends
at Sub-Ambient Temperatures" Frontiers in Energy Research 2022. DOI: 10.3389/fenrg.2022.1074699.
Raw Excel: data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx

Outputs:
  - data/interim/external_staging/boehm2022_raw_extracted.csv

NOTE: Staged output ONLY. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import pandas as pd

RAW_EXCEL = Path("data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx")
STAGING_DIR = Path("data/interim/external_staging")


def main():
    if not RAW_EXCEL.exists():
        print(f"ERROR: {RAW_EXCEL} not found.")
        return

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_excel(RAW_EXCEL, sheet_name="This Work", header=7)

    # Columns of interest
    # 2: 2,2,4,4,6,8,8-heptamethylnonane
    # 3: methylcyclohexane
    # 4: hexylbenzene
    c1_name = "2,2,4,4,6,8,8-heptamethylnonane"
    c2_name = "methylcyclohexane"
    c3_name = "hexylbenzene"

    temp_cols = [
        (-40, 233.15, "Density at -40C, g/cc", "Viscosity at -40C, cSt"),
        (-20, 253.15, "Density at -20C, g/cc", "Viscosity at -20C, cSt"),
        (-10, 263.15, "Density at -10C, g/cc", "Viscosity at -10C, cSt"),
        (0, 273.15, "Density at 0C, g/cc", "Viscosity at 0C, cSt"),
        (15, 288.15, "Density at 15C, g/cc", "Viscosity at 15C, cSt"),
    ]

    records = []
    for idx, row in df.iterrows():
        try:
            v1 = float(row[c1_name]) if pd.notna(row[c1_name]) else 0.0
            v2 = float(row[c2_name]) if pd.notna(row[c2_name]) else 0.0
            v3 = float(row[c3_name]) if pd.notna(row[c3_name]) else 0.0
        except ValueError:
            continue

        tot = v1 + v2 + v3
        if tot <= 0.0:
            continue

        for t_c, t_k, dens_col, visc_col in temp_cols:
            dens = float(row[dens_col]) if dens_col in row and pd.notna(row[dens_col]) else float("nan")
            visc = float(row[visc_col]) if visc_col in row and pd.notna(row[visc_col]) else float("nan")

            if pd.notna(dens) or pd.notna(visc):
                records.append({
                    "row_index": idx,
                    "c1_heptamethylnonane_pct": v1,
                    "c2_methylcyclohexane_pct": v2,
                    "c3_hexylbenzene_pct": v3,
                    "T_C": t_c,
                    "T_K": t_k,
                    "density_g_cc": dens,
                    "viscosity_cSt": visc,
                    "provenance": "Boehm et al. (2022) HEAT Lab",
                })

    df_out = pd.DataFrame(records)
    out_csv = STAGING_DIR / "boehm2022_raw_extracted.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"Staged {len(df_out)} low-T measurement rows -> {out_csv}")
    print(df_out["T_C"].value_counts().sort_index())


if __name__ == "__main__":
    main()
