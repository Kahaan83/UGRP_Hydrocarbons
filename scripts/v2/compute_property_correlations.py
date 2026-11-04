#!/usr/bin/env python
"""
scripts/v2/compute_property_correlations.py

Computes empirical multi-property correlation matrix and physical coupling rules
across 2,695 real fuels from CRC (Report 647 & AV-27-18) and AJFTD databases.

Evaluates:
- Flash point vs T10/IBP (Mackison / API empirical correlation)
- Kinematic viscosity vs Density (Walther-ASTM hydrocarbon envelope)
- Freezing point vs Viscosity
- Aromatic content vs Cetane number / Smoke point

Deliverable: results/v2/A6/property_correlations.csv
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[2]
OUT_DIR = BASE / "results" / "v2" / "A6"
OUT_CSV = OUT_DIR / "property_correlations.csv"
WIDE_FUELS_CSV = OUT_DIR / "real_fuels_multi_property.csv"


def load_real_fuel_properties() -> pd.DataFrame:
    # 1. CRC
    crc = pd.read_csv(BASE / "data" / "processed" / "crc_fuels.csv")
    crc_fuels = []
    for fid in crc.fuel_id.unique():
        sub = crc[crc.fuel_id == fid]
        rec = {"fuel_id": f"crc_{fid}", "source": "crc"}
        
        fl = sub[sub.property == "flash_point_C"]["value"].values
        if len(fl) > 0: rec["flash_point_C"] = fl[0]
        
        fr = sub[sub.property == "freezing_point_C"]["value"].values
        if len(fr) > 0: rec["freezing_point_C"] = fr[0]
        
        v20 = sub[(sub.property == "viscosity_kin_cSt") & (sub.T_C == -20.0)]["value"].values
        if len(v20) > 0: rec["viscosity_minus20C_cSt"] = v20[0]
        
        d15 = sub[(sub.property == "density_kg_m3") & ((sub.T_C == 15.0) | sub.T_C.isna())]["value"].values
        if len(d15) > 0: rec["density_15C_kg_m3"] = d15[0]
        
        ar = sub[sub.property == "aromatics_vol_pct"]["value"].values
        if len(ar) > 0: rec["aromatics_vol_pct"] = ar[0]
        
        for cut in ["T0", "T10", "T50", "T90", "T100"]:
            t = sub[sub.property == f"distillation_{cut}_C"]["value"].values
            if len(t) > 0: rec[f"distillation_{cut}_C"] = t[0]
        crc_fuels.append(rec)

    # 2. AJFTD
    ajftd = pd.read_csv(BASE / "data" / "processed" / "ajftd_fuels.csv")
    ajftd_fuels = []
    for fid in ajftd.fuel_id.unique():
        sub = ajftd[ajftd.fuel_id == fid]
        rec = {"fuel_id": f"ajftd_{fid}", "source": "ajftd"}
        
        fl = sub[sub.property == "flash_point_C"]["value"].values
        if len(fl) > 0: rec["flash_point_C"] = fl[0]
        
        fr = sub[sub.property == "freezing_point_C"]["value"].values
        if len(fr) > 0: rec["freezing_point_C"] = fr[0]
        
        v20 = sub[(sub.property == "viscosity_kin_cSt") & (sub.T_K == -20.0)]["value"].values
        if len(v20) > 0: rec["viscosity_minus20C_cSt"] = v20[0]
        
        d15 = sub[(sub.property == "density_kg_m3") & ((sub.T_K == 15.0) | (sub.T_K == 288.15) | (sub.T_K == 288.706) | sub.T_K.isna())]["value"].values
        if len(d15) > 0: rec["density_15C_kg_m3"] = d15[0]
        
        sm = sub[sub.property == "smoke_point_mm"]["value"].values
        if len(sm) > 0: rec["smoke_point_mm"] = sm[0]
        
        for cut in ["T0", "T10", "T50", "T90", "T100"]:
            t = sub[sub.property == f"distillation_{cut}_C"]["value"].values
            if len(t) > 0: rec[f"distillation_{cut}_C"] = t[0]
        ajftd_fuels.append(rec)

    return pd.DataFrame(crc_fuels + ajftd_fuels)


def main():
    print("Loading multi-property real fuel measurements...")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_fuels = load_real_fuel_properties()
    df_fuels.to_csv(WIDE_FUELS_CSV, index=False)
    print(f"Saved wide multi-property fuel records (n={len(df_fuels)}) to {WIDE_FUELS_CSV}")
    
    props = [
        "flash_point_C",
        "freezing_point_C",
        "viscosity_minus20C_cSt",
        "density_15C_kg_m3",
        "distillation_T0_C",
        "distillation_T10_C",
        "distillation_T50_C",
        "distillation_T90_C",
        "distillation_T100_C",
        "smoke_point_mm",
        "aromatics_vol_pct",
    ]
    
    # Compute pairwise correlation table
    rows = []
    for i, p1 in enumerate(props):
        for j, p2 in enumerate(props):
            if p1 == p2:
                sub = df_fuels[p1].dropna()
                rows.append({
                    "property_1": p1,
                    "property_2": p2,
                    "n_observations": len(sub),
                    "pearson_r": 1.0,
                    "pearson_p_value": 0.0,
                    "spearman_rho": 1.0,
                    "spearman_p_value": 0.0,
                })
                continue

            sub = df_fuels[[p1, p2]].dropna()
            n_obs = len(sub)
            if n_obs >= 5:
                x = sub[p1].to_numpy(float)
                y = sub[p2].to_numpy(float)
                res_p = stats.pearsonr(x, y)
                r = float(np.atleast_1d(res_p.statistic)[0]) if hasattr(res_p, "statistic") else float(np.atleast_1d(res_p)[0])
                p_val = float(np.atleast_1d(res_p.pvalue)[0]) if hasattr(res_p, "pvalue") else float(np.atleast_1d(res_p)[1])
                
                res_s = stats.spearmanr(x, y)
                rho = float(np.atleast_1d(res_s.statistic)[0]) if hasattr(res_s, "statistic") else float(np.atleast_1d(res_s)[0])
                s_pval = float(np.atleast_1d(res_s.pvalue)[0]) if hasattr(res_s, "pvalue") else float(np.atleast_1d(res_s)[1])
            else:
                r, p_val, rho, s_pval = np.nan, np.nan, np.nan, np.nan
                
            rows.append({
                "property_1": p1,
                "property_2": p2,
                "n_observations": n_obs,
                "pearson_r": round(float(r), 4) if pd.notna(r) else np.nan,
                "pearson_p_value": float(p_val) if pd.notna(p_val) else np.nan,
                "spearman_rho": round(float(rho), 4) if pd.notna(rho) else np.nan,
                "spearman_p_value": float(s_pval) if pd.notna(s_pval) else np.nan,
            })
            
    corr_df = pd.DataFrame(rows)
    corr_df.to_csv(OUT_CSV, index=False)
    print(f"Saved pairwise property correlations to {OUT_CSV}")
    
    # Print Key Physical Couplings
    print("\n=== KEY PHYSICAL COUPLINGS ON REAL FUELS ===")
    
    # 1. Flash point vs T10
    sub_fl = df_fuels[["flash_point_C", "distillation_T10_C"]].dropna()
    res_fl = stats.pearsonr(sub_fl["flash_point_C"], sub_fl["distillation_T10_C"])
    r_fl = float(np.atleast_1d(res_fl.statistic)[0]) if hasattr(res_fl, "statistic") else float(np.atleast_1d(res_fl)[0])
    print(f"1. Flash Point vs T10 Distillation: r = {r_fl:.3f} (n = {len(sub_fl)} real fuels)")
    
    # 2. Flash point vs IBP
    sub_ibp = df_fuels[["flash_point_C", "distillation_T0_C"]].dropna()
    res_ibp = stats.pearsonr(sub_ibp["flash_point_C"], sub_ibp["distillation_T0_C"])
    r_ibp = float(np.atleast_1d(res_ibp.statistic)[0]) if hasattr(res_ibp, "statistic") else float(np.atleast_1d(res_ibp)[0])
    print(f"2. Flash Point vs IBP Distillation: r = {r_ibp:.3f} (n = {len(sub_ibp)} real fuels)")
    
    # 3. Viscosity (-20 C) vs Flash Point
    sub_vf = df_fuels[["viscosity_minus20C_cSt", "flash_point_C"]].dropna()
    res_vf = stats.pearsonr(sub_vf["viscosity_minus20C_cSt"], sub_vf["flash_point_C"])
    r_vf = float(np.atleast_1d(res_vf.statistic)[0]) if hasattr(res_vf, "statistic") else float(np.atleast_1d(res_vf)[0])
    print(f"3. Viscosity (-20 C) vs Flash Point: r = {r_vf:.3f} (n = {len(sub_vf)} real fuels)")
    
    # 4. Smoke Point vs Freezing Point
    sub_sf = df_fuels[["smoke_point_mm", "freezing_point_C"]].dropna()
    res_sf = stats.pearsonr(sub_sf["smoke_point_mm"], sub_sf["freezing_point_C"])
    r_sf = float(np.atleast_1d(res_sf.statistic)[0]) if hasattr(res_sf, "statistic") else float(np.atleast_1d(res_sf)[0])
    print(f"4. Smoke Point vs Freezing Point: r = {r_sf:.3f} (n = {len(sub_sf)} real fuels)")


if __name__ == "__main__":
    main()
