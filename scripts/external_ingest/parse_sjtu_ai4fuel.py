#!/usr/bin/env python3
"""parse_sjtu_ai4fuel.py - Parse experimental fuel properties from the SJTU AI4fuel / LG-Transformer database.

Source: Zhang, Huang et al. "LG-Transformer: learned-graph transformer framework enabling
diverse physicochemical properties prediction toward fuel design" Nature Communications 2026.
DOI: 10.1038/s41467-026-73853-z. Figshare: 10.6084/m9.figshare.29571449.
Raw Excel: data/raw/git/AI4fuel/all_data.xlsx

Outputs:
  - data/interim/external_staging/ai4fuel/ (individual property CSVs)

NOTE: Staged output ONLY. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import pandas as pd

RAW_EXCEL = Path("data/raw/git/AI4fuel/all_data.xlsx")
OUT_DIR = Path("data/interim/external_staging/ai4fuel")

# Map sheet names to standard property names and value column identifiers
SHEETS_MAP = {
    "RON": ("ron", "RON_Measured"),
    "MON": ("mon", "MON_Measured"),
    "CN": ("cetane_number_D613", "CN_Measured"),
    "DCN": ("dcn", "DCN"),
    "LFL": ("lfl_vol_pct", "LFL (vol.%)_Measured"),
    "UFL": ("ufl_vol_pct", "UFL (vol.%)_Measured"),
    "FP": ("flash_point_C", "FP (oC)_CRC handbook"),
    "Viscosity": ("viscosity_dyn_mPas", "Viscosity,  mPa*s@20oC_Measured"),
    "Density": ("density_g_cm3", "Density (g/cm3)@15~30oC_CRC handbook"),
    "MP": ("melting_point_C", "Tm (oC)_CRC handbook"),
    "BP": ("boiling_point_C", "Tb (oC)_CRC handbook"),
    "VP": ("vapor_pressure_kPa", "VP@25 C (kPa)_CRC handbook"),
    "YSI": ("ysi", "YSI_Unified_Measured"),
    "LHV": ("lhv_MJ_kg", "LHV"),
    "EOV": ("enthalpy_vaporization_kJ_mol", "Enthalpy of Vaporization(25C)/kJ mol-1"),
}


def main():
    if not RAW_EXCEL.exists():
        print(f"ERROR: {RAW_EXCEL} not found.")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    xls = pd.ExcelFile(RAW_EXCEL)
    print(f"Loaded {RAW_EXCEL}. Processing {len(SHEETS_MAP)} property sheets...")

    total_rows = 0
    for sheet_name, (prop_key, val_col_hint) in SHEETS_MAP.items():
        if sheet_name not in xls.sheet_names:
            print(f"Warning: Sheet {sheet_name} not found in Excel file.")
            continue

        df = pd.read_excel(xls, sheet_name=sheet_name)
        # Find the measured value column
        target_col = None
        for col in df.columns:
            if col == val_col_hint or (val_col_hint.split()[0] in str(col) and "pred" not in str(col).lower()):
                target_col = col
                break

        if target_col is None:
            # fallback to 6th column (index 5)
            target_col = df.columns[5]

        # Extract standard columns: Class, Name, CAS, SMILES, Formula, Measured_Value
        clean_df = pd.DataFrame({
            "class": df["Class"] if "Class" in df else "",
            "name": df["Name"] if "Name" in df else "",
            "cas": df["CAS"] if "CAS" in df else "",
            "smiles": df["SMILES"] if "SMILES" in df else "",
            "formula": df["Formula"] if "Formula" in df else "",
            "value": pd.to_numeric(df[target_col], errors="coerce"),
            "property": prop_key,
            "raw_column": str(target_col),
            "provenance": "SJTU AI4fuel (Nature Commun. 2026)",
        }).dropna(subset=["value", "smiles"])

        out_csv = OUT_DIR / f"ai4fuel_{prop_key}.csv"
        clean_df.to_csv(out_csv, index=False)
        print(f"  [{prop_key:25s}] {len(clean_df):4d} valid measured rows -> {out_csv.name}")
        total_rows += len(clean_df)

    print(f"Successfully extracted {total_rows} total property measurements across {len(SHEETS_MAP)} sheets!")


if __name__ == "__main__":
    main()
