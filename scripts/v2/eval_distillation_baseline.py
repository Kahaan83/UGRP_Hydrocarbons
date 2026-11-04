#!/usr/bin/env python
"""
scripts/v2/eval_distillation_baseline.py

Evaluates baseline distillation prediction models against real fuel measurements:
1. Ideal Boiling Point Distribution (ASTM D2887 / SimDist) as baseline for ASTM D86 batch distillation
   evaluated across 57 real aviation fuels from CRC Report 647.
2. Repeatability of ASTM D86 (Test 1 vs Test 2) on the same 57 fuels as empirical experimental noise floor.
3. Ideal cumulative boiling point curve of the 14-component NIST surrogate against Fortin et al. 2012
   measured Advanced Distillation Curve (ADC) for RP-1 4572.

Outputs: results/v2/A5/distillation_baseline.csv
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
REAL_FUELS_CSV = BASE / "data" / "processed" / "v2" / "distillation_real_fuels.csv"
OUT_DIR = BASE / "results" / "v2" / "A5"
OUT_CSV = OUT_DIR / "distillation_baseline.csv"

# 14-component NIST surrogate model for RP-1 (NISTIR 6646, Table 2 & Table 8)
# Component: (name, mole_pct, boiling_point_C, liquid_density_g_cm3, mol_wt)
RP1_SURROGATE = [
    ("3-ethyl-4,4-dimethyl-2-pentene", 9.98, 138.0, 0.745, 126.24),
    ("2-methylnonane", 2.32, 166.8, 0.726, 142.28),
    ("3-methyldecane", 10.84, 188.0, 0.742, 156.31),
    ("Cyclodecene", 2.11, 195.0, 0.880, 138.25),
    ("4-methyl-4-undecene", 10.45, 208.0, 0.760, 168.32),
    ("1-dodecene", 2.64, 213.4, 0.758, 168.32),
    ("2-methyldecalin", 22.35, 215.0, 0.889, 152.28),
    ("dodecane", 1.93, 216.3, 0.749, 170.33),
    ("1-tridecene", 1.45, 232.0, 0.766, 182.35),
    ("2-methylnaphthalene", 5.10, 241.0, 1.006, 142.20),
    ("Cyclododecane", 4.27, 243.0, 0.864, 168.32),
    ("2,7,10-trimethyldodecane", 10.38, 245.0, 0.772, 212.41),
    ("Heptylcyclohexane", 14.22, 247.0, 0.812, 182.35),
    ("n-hexadecane", 1.95, 287.0, 0.773, 226.44),
]


def compute_surrogate_ideal_distillation() -> dict[str, float]:
    """Compute volume-fraction weighted cumulative distillation curve of RP-1 surrogate."""
    # Compute volume fractions
    vols = []
    for name, mole_pct, bp, rho, mw in RP1_SURROGATE:
        mass = mole_pct * mw
        vol = mass / rho
        vols.append((bp, vol))
    
    # Sort by boiling point
    vols.sort(key=lambda t: t[0])
    tot_vol = sum(v for _, v in vols)
    
    # Build cumulative curve
    cum_vol = 0.0
    curve = []
    for bp, vol in vols:
        cum_vol += (vol / tot_vol) * 100.0
        curve.append((cum_vol, bp))
        
    # Interpolate at standard cuts: T5, T10, T20, T50, T85, T90
    cuts = [5.0, 10.0, 20.0, 50.0, 85.0, 90.0]
    preds = {}
    for c in cuts:
        # Step interpolation (first component exceeding cut)
        for cv, bp in curve:
            if cv >= c:
                preds[f"T{int(c)}_C"] = bp
                break
    return preds


def main():
    print("Loading real fuel distillation dataset...")
    df = pd.read_csv(REAL_FUELS_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    rows = []
    
    # -------------------------------------------------------------------------
    # 1. CRC 647: D2887 (Ideal SimDist) vs D86 (Actual Batch Distillation)
    # -------------------------------------------------------------------------
    crc = df[df["source"] == "crc647"].copy()
    fuels = crc["fuel_id"].unique()
    
    cuts = ["T0", "T10", "T20", "T50", "T90", "FBP"]
    
    for f in fuels:
        f_sub = crc[crc["fuel_id"] == f]
        d86_1 = f_sub[f_sub["method"] == "ASTM D86 (test 1)"]
        d86_2 = f_sub[f_sub["method"] == "ASTM D86 (test 2)"]
        d2887 = f_sub[f_sub["method"] == "ASTM D2887"]
        
        if len(d86_1) == 0 or len(d2887) == 0:
            continue
            
        f_name = d86_1["fuel_name"].values[0]
        
        for cut in cuts:
            col = f"{cut}_C"
            val_d86 = d86_1[col].values[0] if col in d86_1.columns and pd.notna(d86_1[col].values[0]) else np.nan
            val_d2887 = d2887[col].values[0] if col in d2887.columns and pd.notna(d2887[col].values[0]) else np.nan
            val_d86_rep = d86_2[col].values[0] if len(d86_2) > 0 and col in d86_2.columns and pd.notna(d86_2[col].values[0]) else np.nan
            
            if pd.notna(val_d86) and pd.notna(val_d2887):
                rows.append({
                    "evaluation_type": "D2887_ideal_vs_D86_actual",
                    "fuel_id": f,
                    "fuel_name": f_name,
                    "recovery_point": cut,
                    "measured_D86_C": val_d86,
                    "baseline_pred_D2887_C": val_d2887,
                    "error_C": val_d2887 - val_d86,  # pred - actual
                    "d86_test2_repeatability_C": val_d86_rep,
                    "d86_repeatability_diff_C": val_d86_rep - val_d86 if pd.notna(val_d86_rep) else np.nan,
                })
                
    # -------------------------------------------------------------------------
    # 2. RP-1 Surrogate Model vs Fortin 2012 Measured ADC Curve
    # -------------------------------------------------------------------------
    surr_preds = compute_surrogate_ideal_distillation()
    fortin_rp1 = df[(df["source"] == "fortin2012") & (df["fuel_name"].str.contains("RP-1 4572"))]
    
    if len(fortin_rp1) > 0:
        f_row = fortin_rp1.iloc[0]
        for cut in ["T5", "T10", "T20", "T50"]:
            meas_val = f_row[f"{cut}_C"]
            pred_val = surr_preds.get(f"{cut}_C", np.nan)
            if pd.notna(meas_val) and pd.notna(pred_val):
                rows.append({
                    "evaluation_type": "RP1_surrogate_ideal_vs_ADC_actual",
                    "fuel_id": "fortin_rp1_4572",
                    "fuel_name": "RP-1 4572 (Fortin et al. 2012)",
                    "recovery_point": cut,
                    "measured_D86_C": meas_val,
                    "baseline_pred_D2887_C": pred_val,
                    "error_C": pred_val - meas_val,
                    "d86_test2_repeatability_C": np.nan,
                    "d86_repeatability_diff_C": np.nan,
                })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_CSV, index=False)
    print(f"Saved distillation baseline evaluation to {OUT_CSV}")
    
    # -------------------------------------------------------------------------
    # Summary Statistics
    # -------------------------------------------------------------------------
    print("\n=== DISTILLATION BASELINE ACCURACY (ASTM D2887 SimDist vs ASTM D86 Batch) ===")
    d86_eval = out_df[out_df["evaluation_type"] == "D2887_ideal_vs_D86_actual"].copy()
    
    for cut in cuts:
        sub = d86_eval[d86_eval["recovery_point"] == cut]
        err = sub["error_C"].to_numpy()
        rep_err = sub["d86_repeatability_diff_C"].dropna().to_numpy()
        print(f"Cut {cut:5s} (n={len(sub)}):")
        print(f"  Baseline Error (Pred D2887 - True D86): Mean={err.mean():+6.2f} C, MAE={np.abs(err).mean():5.2f} C, RMSE={np.sqrt((err**2).mean()):5.2f} C")
        if len(rep_err) > 0:
            print(f"  Experimental Repeatability (Test2 - Test1): Mean={rep_err.mean():+6.2f} C, MAE={np.abs(rep_err).mean():5.2f} C")

    print("\n=== RP-1 SURROGATE IDEAL BOILING CURVE VS ACTUAL DISTILLATION ===")
    rp1_eval = out_df[out_df["evaluation_type"] == "RP1_surrogate_ideal_vs_ADC_actual"]
    for _, r in rp1_eval.iterrows():
        print(f"Cut {r['recovery_point']:4s}: Measured={r['measured_D86_C']:.1f} C, Surrogate Ideal={r['baseline_pred_D2887_C']:.1f} C, Error={r['error_C']:+.1f} C")


if __name__ == "__main__":
    main()
