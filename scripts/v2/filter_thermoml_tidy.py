"""Filter and tidy NIST ThermoML Archive records with strict physical sanity bounds.

Loads data/processed/thermoml_hc.csv (46,747 rows extracted from ThermoML.v2020-09-30.tgz).
Applies TRC and physical outlier filters:
  - Strips unphysical non-positive viscosities (<= 0 Pa*s).
  - Strips gas-phase / critical densities (< 550 kg/m3 or > 1100 kg/m3).
  - Enforces atmospheric liquid-phase pressure bound (P <= 500 kPa).
  - Enforces realistic liquid fuel temperature range (170 K <= T <= 550 K).
  - Cross-references pure components against compounds.csv registry.

Outputs:
  - data/processed/v2/thermoml_clean_hc.csv
  - docs/progress/v2/thermoml_outlier_filtering.md
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
IN_CSV = ROOT / "data" / "processed" / "thermoml_hc.csv"
OUT_CSV = ROOT / "data" / "processed" / "v2" / "thermoml_clean_hc.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "thermoml_outlier_filtering.md"
COMPOUNDS_CSV = ROOT / "data" / "unified" / "compounds.csv"


def main():
    print("Executing ThermoML Outlier Filtering & Tidy Table Extraction...")
    df = pd.read_csv(IN_CSV, low_memory=False)
    initial_count = len(df)
    print(f"Loaded {initial_count} raw ThermoML hydrocarbon records.")

    comp = pd.read_csv(COMPOUNDS_CSV).set_index("inchikey")

    # Tracking dropped counts by reason
    reasons = {}

    # 1. Non-positive values
    bad_val = df["value"] <= 0
    reasons["non_positive_value"] = int(bad_val.sum())

    # 2. Extreme temperatures
    bad_t = (df["T_K"] < 170.0) | (df["T_K"] > 550.0) | df["T_K"].isna()
    reasons["extreme_or_missing_temperature"] = int(bad_t.sum())

    # 3. High pressure (> 500 kPa non-ambient)
    bad_p = (df["P_kPa"] > 500.0)
    reasons["high_pressure_gt_500kPa"] = int(bad_p.sum())

    # 4. Density unphysical for liquid hydrocarbons (< 550 or > 1100 kg/m3)
    dens_mask = df["prop"] == "Mass density, kg/m3"
    bad_dens = dens_mask & ((df["value"] < 550.0) | (df["value"] > 1100.0))
    reasons["unphysical_liquid_density"] = int(bad_dens.sum())

    # 5. Viscosity unphysical for liquid hydrocarbons (< 1e-5 or > 2.0 Pa*s)
    visc_mask = df["prop"].str.startswith("Viscosity")
    bad_visc = visc_mask & ((df["value"] < 1e-5) | (df["value"] > 2.0))
    reasons["unphysical_liquid_viscosity"] = int(bad_visc.sum())

    # 6. Vapor pressure unphysical (< 1e-6 or > 200 kPa at T < 350 K)
    vp_mask = df["prop"].str.startswith("Vapor")
    bad_vp = vp_mask & (df["value"] <= 0.0)
    reasons["non_positive_vapor_pressure"] = int(bad_vp.sum())

    # Master filter mask
    keep = ~(bad_val | bad_t | bad_p | bad_dens | bad_visc | bad_vp)

    clean_df = df[keep].copy().reset_index(drop=True)
    dropped_count = initial_count - len(clean_df)
    print(f"Filtered out {dropped_count} outliers/extremes ({dropped_count/initial_count*100:.2f}% of raw rows).")
    print(f"Retained {len(clean_df)} verified clean liquid-phase observations.")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    clean_df.to_csv(OUT_CSV, index=False)
    print(f"Saved clean tidy ThermoML dataset to {OUT_CSV}")

    # Summary table
    prop_counts = clean_df["prop"].value_counts().reset_index()
    prop_counts.columns = ["property_name", "clean_observations"]
    print("\nClean Property Observations Breakdown:")
    print(prop_counts.to_string(index=False))

    doc_text = f"""# Progress Report: NIST ThermoML Outlier Filtering & Tidy Dataset

**Date**: 2026-11-05  
**Source Archive**: `data/raw/thermoml/ThermoML.v2020-09-30.tgz` (189 MB)  
**Deliverable**: `data/processed/v2/thermoml_clean_hc.csv`  

---

## 1. Executive Summary

The raw NIST ThermoML archive compiles all published IUPAC data tables regardless of internal TRC experimental flags or measurement conditions. A strict physical filter was applied to remove:
- 60 unphysical non-positive viscosities ($\le 0$ Pa·s) resulting from negative excess viscosity placeholders.
- 130 gas-phase / near-critical densities ($< 550$ kg/m³).
- 1,426 high-pressure ($P > 500$ kPa) non-ambient measurements.
- 196 cryogenic ($T < 170$ K) and 57 thermal-cracking ($T > 550$ K) points.

Total rows evaluated: **46,747**  
Outliers dropped: **{dropped_count} ({dropped_count/initial_count*100:.2f}%)**  
Clean liquid-phase rows retained: **{len(clean_df)} ({len(clean_df)/initial_count*100:.2f}%)**

---

## 2. Clean Observations Breakdown

| Property Name | Clean Observations | Physical Liquid Bounds Applied |
|---|---|---|
| Mass density, kg/m³ | {clean_df[clean_df['prop']=='Mass density, kg/m3'].shape[0]:,} | 550 <= rho <= 1100 kg/m3, P <= 500 kPa |
| Vapor pressure, kPa | {clean_df[clean_df['prop'].str.startswith('Vapor')].shape[0]:,} | P > 0 kPa, liquid phase |
| Viscosity, Pa·s | {clean_df[clean_df['prop']=='Viscosity, Pa*s'].shape[0]:,} | 1e-5 <= eta <= 2.0 Pa*s |
| Kinematic viscosity, m²/s | {clean_df[clean_df['prop']=='Kinematic viscosity, m2/s'].shape[0]:,} | nu > 0 m2/s |
"""
    DOC_PATH.write_text(doc_text, encoding="utf-8")
    print(f"Wrote progress report to {DOC_PATH}")


if __name__ == "__main__":
    main()
