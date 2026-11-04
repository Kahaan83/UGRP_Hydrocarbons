import pandas as pd, numpy as np, json, joblib
from pathlib import Path

BASE = Path('.')
preds = pd.read_csv('data/runs/baseline_preds.csv')
mt = pd.read_csv('data/processed/model_table.csv')
locked = json.load(open('splits/locked.json'))
lk_sys = set(locked.get('multi_system_ids_locked', []))

import sys; sys.path.insert(0, 'src')
from fuelblend.features.mixture_features import extract_mixture_features

comp_df = pd.read_csv('data/unified/compounds.csv')
cd = {str(r.inchikey): {'mw': float(r.mw) if pd.notna(r.mw) else 100.0,
                         'n_carbon': float(r.n_carbon) if pd.notna(r.n_carbon) else 7.0,
                         'cls': str(r.cls) if pd.notna(r.cls) else 'other'}
      for _, r in comp_df.iterrows()}

mt_clean = mt[~mt['system_id'].isin(lk_sys)].copy()
mt_clean['T_K_r4'] = mt_clean['T_K'].round(4).astype(str)
mt_lookup = mt_clean[['system_id','property','T_K_r4','comp_inchikeys','w_mass','x_mole','quality_weight']].drop_duplicates(subset=['system_id','property','T_K_r4'])
preds['T_K_r4'] = pd.to_numeric(preds['T_K'], errors='coerce').round(4).astype(str)
merged = preds.merge(mt_lookup, on=['system_id','property','T_K_r4'], how='inner')
print(f'Merged: {len(merged)} rows')

for prop, base_col, log_mode in [('viscosity_dyn_mPas', 'gn_mass', True), ('density_kg_m3', 'dens_idvol', False)]:
    mask = ((merged['property']==prop) & (merged['conflict']==False) &
            merged[base_col].notna() & merged['comp_inchikeys'].notna())
    if log_mode:
        mask &= (merged[base_col] > 0) & (merged['value'] > 0)
    data = merged[mask].copy()
    feats = [extract_mixture_features(r, cd) for _, r in data.iterrows()]
    X = pd.DataFrame(feats)
    if log_mode:
        residual = np.log(data['value'].values) - np.log(data[base_col].values)
    else:
        residual = data['value'].values - data[base_col].values
    s2 = (data['split'] == 'S2').values
    s2_preds = np.mean([joblib.load(f'runs/mixture_models/{prop}__seed_{s}.joblib').predict(X[s2]) for s in [11,23,37,41,59]], axis=0)
    rmse_base = np.sqrt(np.mean(residual[s2]**2))
    rmse_mod = np.sqrt(np.mean((residual[s2] - s2_preds)**2))
    print(f'{prop} S2: baseline={rmse_base:.4f} residual={rmse_mod:.4f} (n_s2={s2.sum()})')
