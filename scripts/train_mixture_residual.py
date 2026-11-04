import os
import sys
import json
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

sys.path.insert(0, os.path.abspath("."))
from src.fuelblend.features.mixture_features import extract_mixture_features

def main():
    print("=== Training Mixture Residual Models (WP5) ===")
    
    # 1. Load compounds metadata
    compounds_df = pd.read_csv("data/unified/compounds.csv")
    compounds_dict = {}
    for _, r in compounds_df.iterrows():
        compounds_dict[r['inchikey']] = {
            'mw': float(r['mw']) if pd.notna(r['mw']) else 100.0,
            'n_carbon': float(r['n_carbon']) if pd.notna(r['n_carbon']) else 7.0,
            'cls': str(r['cls']) if pd.notna(r['cls']) else 'other',
        }
        
    # 2. Load model table and baseline predictions
    model_df = pd.read_csv("data/processed/model_table.csv")
    locked = json.load(open("splits/locked.json"))
    lk_sys = set(locked.get("multi_system_ids_locked", []))
    
    # Exactly mirror apply_splits filtering from baselines_eval.py
    mt_clean = model_df[~model_df['system_id'].isin(lk_sys)]
    target_mt = mt_clean[mt_clean['property'].isin([
        'viscosity_dyn_mPas', 'viscosity_kin_cSt', 'density_kg_m3', 'flash_point_C', 'vapor_pressure_kPa'
    ])].reset_index(drop=True)
    
    preds_df = pd.read_csv("data/runs/baseline_preds.csv")
    assert len(target_mt) == len(preds_df), f"Row count mismatch: {len(target_mt)} vs {len(preds_df)}"
    
    # Transfer composition columns directly by exact row alignment
    preds_df['comp_inchikeys'] = target_mt['comp_inchikeys']
    preds_df['w_mass'] = target_mt['w_mass']
    preds_df['x_mole'] = target_mt['x_mole']
    preds_df['quality_weight'] = target_mt['quality_weight']
    
    print(f"Verified 1:1 alignment across all {len(preds_df)} prediction records.")
    
    seeds = [11, 23, 37, 41, 59]
    
    # -------------------------------------------------------------
    # Viscosity Residual Model: Target = ln(value) - ln(gn_mass)
    # -------------------------------------------------------------
    print("\n--- Fitting Viscosity Residual Model ---")
    visc_mask = (
        (preds_df['property'] == 'viscosity_dyn_mPas') &
        (preds_df['conflict'] == False) &
        (preds_df['gn_mass'].notna()) &
        (preds_df['gn_mass'] > 0) &
        (preds_df['value'] > 0) &
        (preds_df['comp_inchikeys'].notna())
    )
    v_data = preds_df[visc_mask].copy()
    print(f"Total valid viscosity mixture rows: {len(v_data)}")
    
    feat_list = []
    for _, row in v_data.iterrows():
        feat_list.append(extract_mixture_features(row, compounds_dict))
    X_v = pd.DataFrame(feat_list)
    
    y_target_log = np.log(v_data['value'].values)
    y_base_log = np.log(v_data['gn_mass'].values)
    residual_log = y_target_log - y_base_log
    
    train_mask = (v_data['split'] == 'S1').values
    test_mask = (v_data['split'] == 'S2').values
    print(f"Viscosity S1 train rows: {train_mask.sum()}, S2 test rows: {test_mask.sum()}")
    
    v_models = []
    for s in seeds:
        gbr = HistGradientBoostingRegressor(
            max_iter=150,
            max_depth=6,
            min_samples_leaf=20,
            learning_rate=0.05,
            random_state=s
        )
        sample_w = v_data['quality_weight'].values[train_mask]
        sample_w = np.nan_to_num(sample_w, nan=1.0)
        gbr.fit(X_v[train_mask], residual_log[train_mask], sample_weight=sample_w)
        v_models.append(gbr)
        
    s1_res_pred = np.mean([m.predict(X_v[train_mask]) for m in v_models], axis=0)
    s2_res_pred = np.mean([m.predict(X_v[test_mask]) for m in v_models], axis=0)
    
    s1_base_rmse = np.sqrt(np.mean((residual_log[train_mask])**2))
    s1_model_rmse = np.sqrt(np.mean((residual_log[train_mask] - s1_res_pred)**2))
    
    s2_base_rmse = np.sqrt(np.mean((residual_log[test_mask])**2))
    s2_model_rmse = np.sqrt(np.mean((residual_log[test_mask] - s2_res_pred)**2))
    
    print(f"[Viscosity S1] Baseline ln-RMSE: {s1_base_rmse:.4f} -> Residual Model ln-RMSE: {s1_model_rmse:.4f}")
    print(f"[Viscosity S2] Baseline ln-RMSE: {s2_base_rmse:.4f} -> Residual Model ln-RMSE: {s2_model_rmse:.4f}")
    
    # -------------------------------------------------------------
    # Density Residual Model: Target = value - dens_idvol
    # -------------------------------------------------------------
    print("\n--- Fitting Density Residual Model ---")
    dens_mask = (
        (preds_df['property'] == 'density_kg_m3') &
        (preds_df['conflict'] == False) &
        (preds_df['dens_idvol'].notna()) &
        (preds_df['value'].notna()) &
        (preds_df['comp_inchikeys'].notna())
    )
    d_data = preds_df[dens_mask].copy()
    print(f"Total valid density mixture rows: {len(d_data)}")
    
    d_feat_list = []
    for _, row in d_data.iterrows():
        d_feat_list.append(extract_mixture_features(row, compounds_dict))
    X_d = pd.DataFrame(d_feat_list)
    
    y_target_d = d_data['value'].values
    y_base_d = d_data['dens_idvol'].values
    residual_d = y_target_d - y_base_d
    
    d_train_mask = (d_data['split'] == 'S1').values
    d_test_mask = (d_data['split'] == 'S2').values
    print(f"Density S1 train rows: {d_train_mask.sum()}, S2 test rows: {d_test_mask.sum()}")
    
    d_models = []
    for s in seeds:
        gbr = HistGradientBoostingRegressor(
            max_iter=150,
            max_depth=5,
            min_samples_leaf=20,
            learning_rate=0.05,
            random_state=s
        )
        sample_w = d_data['quality_weight'].values[d_train_mask]
        sample_w = np.nan_to_num(sample_w, nan=1.0)
        gbr.fit(X_d[d_train_mask], residual_d[d_train_mask], sample_weight=sample_w)
        d_models.append(gbr)
        
    d_s1_pred = np.mean([m.predict(X_d[d_train_mask]) for m in d_models], axis=0)
    d_s2_pred = np.mean([m.predict(X_d[d_test_mask]) for m in d_models], axis=0)
    
    d_s1_base_rmse = np.sqrt(np.mean((residual_d[d_train_mask])**2))
    d_s1_model_rmse = np.sqrt(np.mean((residual_d[d_train_mask] - d_s1_pred)**2))
    
    d_s2_base_rmse = np.sqrt(np.mean((residual_d[d_test_mask])**2))
    d_s2_model_rmse = np.sqrt(np.mean((residual_d[d_test_mask] - d_s2_pred)**2))
    
    print(f"[Density S1] Baseline RMSE (kg/m3): {d_s1_base_rmse:.3f} -> Residual Model RMSE: {d_s1_model_rmse:.3f}")
    print(f"[Density S2] Baseline RMSE (kg/m3): {d_s2_base_rmse:.3f} -> Residual Model RMSE: {d_s2_model_rmse:.3f}")

    # Write summary log
    log_content = f"""# RESIDUAL_EVAL_LOG.md - WP5 Mixture Residual Models

Evaluated on 5-seed ensemble of HistGradientBoostingRegressor with sample weights.
Residual definition:
- Viscosity: ln(measured) - ln(gn_mass baseline)
- Density: measured - ideal volume baseline (kg/m3)

## Viscosity (viscosity_dyn_mPas, ln-RMSE)
- S1 Train (n={train_mask.sum()}): Baseline = {s1_base_rmse:.4f} -> Residual Model = {s1_model_rmse:.4f} ({(1 - s1_model_rmse/s1_base_rmse)*100:.1f} % error reduction)
- S2 Test (n={test_mask.sum()}): Baseline = {s2_base_rmse:.4f} -> Residual Model = {s2_model_rmse:.4f} ({(1 - s2_model_rmse/s2_base_rmse)*100:.1f} % error reduction)

## Density (density_kg_m3, RMSE kg/m3)
- S1 Train (n={d_train_mask.sum()}): Baseline = {d_s1_base_rmse:.3f} -> Residual Model = {d_s1_model_rmse:.3f} ({(1 - d_s1_model_rmse/d_s1_base_rmse)*100:.1f} % error reduction)
- S2 Test (n={d_test_mask.sum()}): Baseline = {d_s2_base_rmse:.3f} -> Residual Model = {d_s2_model_rmse:.3f} ({(1 - d_s2_model_rmse/d_s2_base_rmse)*100:.1f} % error reduction)
"""
    with open("data/runs/RESIDUAL_EVAL_LOG.md", "w") as f:
        f.write(log_content)
    print("\nWrote data/runs/RESIDUAL_EVAL_LOG.md successfully.")

if __name__ == "__main__":
    main()
