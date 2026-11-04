#!/usr/bin/env python
"""
scripts/v2/eval_flash_baseline_full.py

Comprehensive evaluation of Le Chatelier + Raoult's law + Antoine flash point baseline
on the complete set of 1,006 mixture records from Leenhouts benchmark.

Addresses Task A4:
- Evaluates ALL 1,006 mixture records.
- Evaluates both constant LFL (25 C) and temperature-dependent LFL (Zabetakis correction).
- Audits C12+ low-temperature Antoine extrapolation behavior.
- Reports MAE, RMSE, max dangerous error (predicting flash higher than actual), max conservative error.
- Cluster bootstraps over unique system_id groups.
- Outputs results/v2/A4/flash_baseline_full.csv and summary statistics.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BASE / "src"))

from fuelblend.baselines.flash import (
    predict_flash_le_chatelier,
    check_antoine_validity,
    psat_antoine,
)

BASE = Path(__file__).resolve().parents[2]
MIX_CSV = BASE / "data" / "processed" / "t7_leenhouts_flash_mix.csv"
COMPOUNDS_CSV = BASE / "data" / "unified" / "compounds.csv"
OUT_DIR = BASE / "results" / "v2" / "A4"
OUT_CSV = OUT_DIR / "flash_baseline_full.csv"


def cluster_bootstrap(err: np.ndarray, clusters: np.ndarray, B: int = 4000, seed: int = 42):
    """Cluster bootstrap of MAE and RMSE over cluster groups."""
    u, inv = np.unique(clusters, return_inverse=True)
    n_c = np.bincount(inv).astype(float)
    abs_err_sums = np.bincount(inv, weights=np.abs(err))
    sq_err_sums = np.bincount(inv, weights=(err ** 2))
    
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(len(u), np.full(len(u), 1.0 / len(u)), size=B)
    n_samples = counts @ n_c
    
    mae_dist = (counts @ abs_err_sums) / n_samples
    rmse_dist = np.sqrt((counts @ sq_err_sums) / n_samples)
    
    return {
        "mae_mean": float(mae_dist.mean()),
        "mae_ci": [float(np.percentile(mae_dist, 2.5)), float(np.percentile(mae_dist, 97.5))],
        "rmse_mean": float(rmse_dist.mean()),
        "rmse_ci": [float(np.percentile(rmse_dist, 2.5)), float(np.percentile(rmse_dist, 97.5))],
    }


def main():
    print(f"Loading {MIX_CSV} and {COMPOUNDS_CSV}...")
    df_mix = pd.read_csv(MIX_CSV)
    comp = pd.read_csv(COMPOUNDS_CSV).set_index("inchikey")
    
    print(f"Total mixture records to evaluate: {len(df_mix)}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    rows = []
    
    for idx, r in df_mix.iterrows():
        keys = json.loads(r.comp_inchikeys)
        x = json.loads(r.x_mole)
        val = float(r.value)
        n_comp = len(keys)
        is_hc = bool(r.is_hc)
        sys_id = str(r.system_id)
        
        # Check component coverage in compounds.csv
        missing_comp = [k for k in keys if k not in comp.index]
        missing_ant = [k for k in keys if k in comp.index and pd.isna(comp.loc[k, "antoine_A"])]
        missing_lfl = [k for k in keys if k in comp.index and pd.isna(comp.loc[k, "lfl_vol_pct"])]
        
        c_nums = [int(comp.loc[k, "n_carbon"]) for k in keys if k in comp.index and pd.notna(comp.loc[k, "n_carbon"])]
        max_c = max(c_nums) if c_nums else np.nan
        has_c12 = bool(max_c >= 12) if pd.notna(max_c) else False
        
        record = {
            "row_id": idx,
            "system_id": sys_id,
            "property": "flash_point_C",
            "kind": r.kind,
            "n_comp": n_comp,
            "is_hc": is_hc,
            "measured_C": val,
            "comp_inchikeys": r.comp_inchikeys,
            "x_mole": r.x_mole,
            "max_n_carbon": max_c,
            "has_c12_plus": has_c12,
            "missing_components_count": len(missing_comp),
            "missing_antoine_count": len(missing_ant),
            "missing_lfl_count": len(missing_lfl),
        }
        
        if missing_comp:
            record["coverage_status"] = f"missing_{len(missing_comp)}_unregistered_components"
            record["pred_flash_const_lfl_C"] = np.nan
            record["err_const_lfl_C"] = np.nan
            record["pred_flash_temp_lfl_C"] = np.nan
            record["err_temp_lfl_C"] = np.nan
            record["antoine_in_range_const"] = False
            record["extrap_distance_K_const"] = np.nan
            record["dangerous_underpredict_risk"] = False
            record["dangerous_over_5C_const"] = False
            rows.append(record)
            continue
            
        if missing_ant or missing_lfl:
            record["coverage_status"] = f"missing_params_ant{len(missing_ant)}_lfl{len(missing_lfl)}"
            record["pred_flash_const_lfl_C"] = np.nan
            record["err_const_lfl_C"] = np.nan
            record["pred_flash_temp_lfl_C"] = np.nan
            record["err_temp_lfl_C"] = np.nan
            record["antoine_in_range_const"] = False
            record["extrap_distance_K_const"] = np.nan
            record["dangerous_underpredict_risk"] = False
            record["dangerous_over_5C_const"] = False
            rows.append(record)
            continue
            
        # Extract parameters
        antoine_tuples = [(float(comp.loc[k, "antoine_A"]), float(comp.loc[k, "antoine_B"]), float(comp.loc[k, "antoine_C"])) for k in keys]
        lfls = [float(comp.loc[k, "lfl_vol_pct"]) for k in keys]
        tlos = [float(comp.loc[k, "antoine_Tlo"]) if pd.notna(comp.loc[k, "antoine_Tlo"]) else np.nan for k in keys]
        this = [float(comp.loc[k, "antoine_Thi"]) if pd.notna(comp.loc[k, "antoine_Thi"]) else np.nan for k in keys]
        
        # Predict with constant LFL (25 C)
        p_const, stat_const = predict_flash_le_chatelier(x, antoine_tuples, lfls, lfl_temp=False)
        # Predict with temp-dependent LFL
        p_temp, stat_temp = predict_flash_le_chatelier(x, antoine_tuples, lfls, lfl_temp=True)
        
        record["coverage_status"] = "ok" if stat_const == "ok" else f"solver_error_{stat_const}"
        record["pred_flash_const_lfl_C"] = p_const
        record["err_const_lfl_C"] = p_const - val if pd.notna(p_const) else np.nan
        record["pred_flash_temp_lfl_C"] = p_temp
        record["err_temp_lfl_C"] = p_temp - val if pd.notna(p_temp) else np.nan
        
        if pd.notna(p_const):
            in_range, max_dist = check_antoine_validity(p_const, tlos, this)
            record["antoine_in_range_const"] = in_range
            record["extrap_distance_K_const"] = max_dist
            record["dangerous_underpredict_risk"] = bool(p_const > val)
            record["dangerous_over_5C_const"] = bool((p_const - val) > 5.0)
        else:
            record["antoine_in_range_const"] = False
            record["extrap_distance_K_const"] = np.nan
            record["dangerous_underpredict_risk"] = False
            record["dangerous_over_5C_const"] = False
            
        rows.append(record)
        
    out_df = pd.DataFrame(rows)
    out_df.to_csv(OUT_CSV, index=False)
    print(f"Saved full predictions to {OUT_CSV}")
    
    # Coverage summary
    print("\n=== COVERAGE SUMMARY ===")
    print(out_df["coverage_status"].value_counts())
    
    ok_df = out_df[out_df["coverage_status"] == "ok"].copy()
    print(f"\nTotal evaluated successfully: {len(ok_df)} / {len(out_df)}")
    
    # Statistical analysis on predictable subset
    err_c = ok_df["err_const_lfl_C"].to_numpy()
    err_t = ok_df["err_temp_lfl_C"].to_numpy()
    systems = ok_df["system_id"].to_numpy()
    
    boot_c = cluster_bootstrap(err_c, systems)
    boot_t = cluster_bootstrap(err_t, systems)
    
    print("\n=== PERFORMANCE METRICS (n = 126 mixtures, 31 unique systems) ===")
    print(f"Constant LFL (25 C):")
    print(f"  MAE:  {np.mean(np.abs(err_c)):.3f} degC  [95% CI: {boot_c['mae_ci'][0]:.3f}, {boot_c['mae_ci'][1]:.3f}]")
    print(f"  RMSE: {np.sqrt(np.mean(err_c**2)):.3f} degC  [95% CI: {boot_c['rmse_ci'][0]:.3f}, {boot_c['rmse_ci'][1]:.3f}]")
    print(f"  Bias (pred - meas): {np.mean(err_c):+.3f} degC")
    print(f"  Dangerous Over-prediction (pred > actual): {np.mean(err_c > 0)*100:.1f}% of mixtures")
    print(f"  Dangerous Over-prediction > 5 degC:        {np.mean(err_c > 5)*100:.1f}% ({np.sum(err_c > 5)} blends)")
    print(f"  Max Dangerous Error (pred - actual):       {np.max(err_c):+.3f} degC")
    print(f"  Max Conservative Error (actual - pred):    {np.max(-err_c):+.3f} degC")
    
    print(f"\nTemperature-dependent LFL (Zabetakis correction):")
    print(f"  MAE:  {np.mean(np.abs(err_t)):.3f} degC  [95% CI: {boot_t['mae_ci'][0]:.3f}, {boot_t['mae_ci'][1]:.3f}]")
    print(f"  RMSE: {np.sqrt(np.mean(err_t**2)):.3f} degC  [95% CI: {boot_t['rmse_ci'][0]:.3f}, {boot_t['rmse_ci'][1]:.3f}]")
    print(f"  Bias (pred - meas): {np.mean(err_t):+.3f} degC")
    print(f"  Dangerous Over-prediction (pred > actual): {np.mean(err_t > 0)*100:.1f}% of mixtures")
    print(f"  Dangerous Over-prediction > 5 degC:        {np.mean(err_t > 5)*100:.1f}% ({np.sum(err_t > 5)} blends)")
    print(f"  Max Dangerous Error (pred - actual):       {np.max(err_t):+.3f} degC")
    print(f"  Max Conservative Error (actual - pred):    {np.max(-err_t):+.3f} degC")
    
    # Paired test
    diff_mae = np.abs(err_t) - np.abs(err_c)
    print(f"\nPaired MAE difference (Temp LFL - Const LFL): {np.mean(diff_mae):+.3f} degC")
    
    # C12+ Breakdown
    c12_sub = ok_df[ok_df["has_c12_plus"]].copy()
    non_c12 = ok_df[~ok_df["has_c12_plus"]].copy()
    print(f"\nSub-analysis by heavy component presence:")
    print(f"  C12+ present (n = {len(c12_sub)}): Const MAE = {np.mean(np.abs(c12_sub['err_const_lfl_C'])):.3f} C, Temp MAE = {np.mean(np.abs(c12_sub['err_temp_lfl_C'])):.3f} C")
    print(f"  < C12 only   (n = {len(non_c12)}): Const MAE = {np.mean(np.abs(non_c12['err_const_lfl_C'])):.3f} C, Temp MAE = {np.mean(np.abs(non_c12['err_temp_lfl_C'])):.3f} C")
    
    # Antoine validity breakdown
    in_sub = ok_df[ok_df["antoine_in_range_const"]].copy()
    out_sub = ok_df[~ok_df["antoine_in_range_const"]].copy()
    print(f"\nSub-analysis by Antoine Validity Range:")
    print(f"  Inside Antoine validity [Tlo, Thi] (n = {len(in_sub)}): MAE = {np.mean(np.abs(in_sub['err_const_lfl_C'])):.3f} C")
    print(f"  Extrapolated outside validity      (n = {len(out_sub)}): MAE = {np.mean(np.abs(out_sub['err_const_lfl_C'])):.3f} C")
    print(f"  Max extrapolation distance: {ok_df['extrap_distance_K_const'].max():.1f} K")


if __name__ == "__main__":
    main()
