import json
import numpy as np
import pandas as pd

def extract_mixture_features(row, compounds_dict):
    """
    Extract blend-level tabular features from composition vectors.
    Permutation-invariant weighted sums and standard deviations over components.
    """
    try:
        inchis = json.loads(row['comp_inchikeys'])
        w = json.loads(row['w_mass'])
        x = json.loads(row['x_mole'])
    except Exception:
        return {}

    w = np.array(w, dtype=float)
    x = np.array(x, dtype=float)
    
    # Normalize fractions if slight numerical drift
    if w.sum() > 0:
        w = w / w.sum()
    if x.sum() > 0:
        x = x / x.sum()

    n_comp = len(inchis)
    
    mws = []
    n_carbons = []
    classes = []
    
    for k in inchis:
        info = compounds_dict.get(k, {})
        mws.append(info.get('mw', np.nan))
        n_carbons.append(info.get('n_carbon', np.nan))
        classes.append(info.get('cls', 'other'))

    mws = np.array(mws, dtype=float)
    n_carbons = np.array(n_carbons, dtype=float)
    
    # Mass-weighted statistics
    mw_mean = np.nansum(w * mws)
    mw_std = np.sqrt(np.nansum(w * (mws - mw_mean)**2)) if n_comp > 1 else 0.0
    
    c_mean = np.nansum(w * n_carbons)
    c_std = np.sqrt(np.nansum(w * (n_carbons - c_mean)**2)) if n_comp > 1 else 0.0
    
    # Class fractions (mass-based)
    class_types = ['aromatic', 'alkene_alkyne', 'cycloalkane', 'iso_paraffin', 'n_paraffin']
    class_w = {f'frac_{cls}': 0.0 for cls in class_types}
    for cls, weight in zip(classes, w):
        if cls in class_types:
            class_w[f'frac_{cls}'] += weight
            
    # Entropy of mixing (ideal proxy)
    entropy_mix = -np.sum(x * np.log(np.clip(x, 1e-12, 1.0)))

    feats = {
        'n_comp': float(n_comp),
        'mw_mean': mw_mean,
        'mw_std': mw_std,
        'c_mean': c_mean,
        'c_std': c_std,
        'entropy_mix': entropy_mix,
        'T_K': float(row.get('T_K', 298.15)),
        'inv_T': 1000.0 / float(row.get('T_K', 298.15)) if float(row.get('T_K', 298.15)) > 0 else 0.0,
    }
    feats.update(class_w)
    return feats
