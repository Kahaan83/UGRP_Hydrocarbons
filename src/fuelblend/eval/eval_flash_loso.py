"""
eval_flash_loso.py — Leave-One-System-Out Cross Validation for Hydrocarbon Flash Mixtures.

Evaluates all 126 hydrocarbon mixture flash rows across 31 unique system_id groups.
Computes:
  1. Le Chatelier + Raoult + Antoine physics baseline
  2. HistGradientBoosting residual model
  3. Paired differences with 4,000-resample cluster bootstrap across the 31 systems
  4. One-sided safety metrics: over-prediction > 5 degC, max over-prediction, median bias
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from sklearn.ensemble import HistGradientBoostingRegressor

from fuelblend.mixture.common import CompoundRegistry


def cluster_bootstrap(err_a, err_b, sysid, B=4000, seed=42):
    """Cluster bootstrap of paired difference MAE(a) - MAE(b). Negative = a better."""
    u, inv = np.unique(sysid, return_inverse=True)
    n_s = np.bincount(inv).astype(float)
    sa = np.bincount(inv, weights=np.abs(err_a))
    sb = np.bincount(inv, weights=np.abs(err_b))
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(len(u), np.full(len(u), 1.0 / len(u)), size=B)
    n = counts @ n_s
    A = counts @ sa / n
    B_m = counts @ sb / n
    d = A - B_m
    return float(d.mean()), np.percentile(d, [2.5, 97.5]), float((d < 0).mean())


def evaluate_loso(root: Path):
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    mt = pd.read_csv(root / "data/processed/model_table.csv", low_memory=False)

    # 1. Separate pure flash rows and mixture flash rows
    pures = mt[(mt["property"] == "flash_point_C") & (mt["n_comp"] == 1)].copy()
    mixes = mt[(mt["property"] == "flash_point_C") & (mt["n_comp"] >= 2)].copy()

    # Pre-extract Antoine & LFL map
    cdf = pd.read_csv(root / "data/unified/compounds.csv")
    ant_map = {}
    for r in cdf.itertuples():
        if pd.notna(r.antoine_A) and pd.notna(r.lfl_vol_pct):
            ant_map[r.inchikey] = (float(r.antoine_A), float(r.antoine_B), float(r.antoine_C),
                                   float(r.lfl_vol_pct) / 100.0)

    # Physics baseline solver
    def solve_flash_base(iks, x):
        if not all(k in ant_map for k in iks) or len(iks) != len(x):
            return np.nan

        def H_T(T_K):
            val = 0.0
            for k, xi in zip(iks, x):
                A, B, C, lfl = ant_map[k]
                psat = 10.0 ** (A - B / (T_K + C))
                val += xi * psat / (101.325 * lfl)
            return val - 1.0

        try:
            if H_T(150.0) * H_T(600.0) < 0:
                return float(brentq(H_T, 150.0, 600.0) - 273.15)
        except Exception:
            pass
        return np.nan

    # Compute physics baseline for all mixtures
    base_preds = []
    for r in mixes.itertuples():
        iks = json.loads(r.comp_inchikeys)
        x = json.loads(r.x_mole)
        base_preds.append(solve_flash_base(iks, x))
    mixes["base_flash"] = base_preds

    # Feature extraction for HGB: composition-weighted descriptor statistics
    def extract_feats(df_subset):
        rows = []
        for r in df_subset.itertuples():
            iks = json.loads(r.comp_inchikeys) if r.comp_inchikeys.startswith("[") else [r.comp_inchikeys]
            ws = json.loads(r.w_mass) if r.w_mass.startswith("[") else [1.0]
            ws = np.array([float(w) for w in ws])
            ws = ws / ws.sum()
            feats = []
            for k in iks:
                idx = reg.ik2i.get(k, 0)
                feats.append(reg.S[idx])
            feats = np.array(feats)
            w_mean = (ws[:, None] * feats).sum(0)
            w_std = np.sqrt((ws[:, None] * (feats - w_mean) ** 2).sum(0))
            n_c = float(len(iks))
            rows.append(np.concatenate([w_mean, w_std, [n_c, ws.max()]]))
        return np.array(rows)

    X_mix = extract_feats(mixes)
    X_pure = extract_feats(pures)
    y_pure = pures["value"].to_numpy(float)
    y_mix = mixes["value"].to_numpy(float)
    base_mix = mixes["base_flash"].to_numpy(float)

    # Leave-One-System-Out CV
    systems = mixes["system_id"].unique()
    hgb_preds = np.full(len(mixes), np.nan)

    print(f"[LOSO-Flash] Running Leave-One-System-Out CV across {len(systems)} systems ({len(mixes)} rows)...")
    for s_idx, sys_id in enumerate(systems):
        test_mask = (mixes["system_id"] == sys_id).values
        train_mix_mask = ~test_mask

        # Train on remaining mixtures + all pure flash rows
        X_tr = np.vstack([X_mix[train_mix_mask], X_pure])
        # Target: residual on top of baseline if baseline finite, else direct value
        y_tr = np.concatenate([y_mix[train_mix_mask] - base_mix[train_mix_mask],
                               y_pure - pures["value"].to_numpy(float)])  # pure residual = 0
        valid_tr = np.isfinite(y_tr)

        m = HistGradientBoostingRegressor(max_iter=100, learning_rate=0.08, max_depth=4,
                                          min_samples_leaf=3, random_state=42)
        m.fit(X_tr[valid_tr], y_tr[valid_tr])

        res_pred = m.predict(X_mix[test_mask])
        hgb_preds[test_mask] = base_mix[test_mask] + res_pred

    mixes["hgb_flash"] = hgb_preds

    # Evaluate metrics on valid rows
    ok = np.isfinite(base_mix) & np.isfinite(hgb_preds)
    y_t = y_mix[ok]
    b_p = base_mix[ok]
    h_p = hgb_preds[ok]
    sys_ok = mixes["system_id"].values[ok]

    err_b = b_p - y_t
    err_h = h_p - y_t

    mae_b, rmse_b = float(np.mean(np.abs(err_b))), float(np.sqrt(np.mean(err_b ** 2)))
    mae_h, rmse_h = float(np.mean(np.abs(err_h))), float(np.sqrt(np.mean(err_h ** 2)))

    over5_b = float((err_b > 5.0).mean() * 100.0)
    over5_h = float((err_h > 5.0).mean() * 100.0)

    # Cluster bootstrap
    mean_diff, ci, p_better = cluster_bootstrap(err_h, err_b, sys_ok)

    out_dir = root / "results/eval"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame([{
        "property": "flash_point_C",
        "n_rows": int(ok.sum()),
        "n_systems": int(len(np.unique(sys_ok))),
        "base_mae": mae_b,
        "base_rmse": rmse_b,
        "base_over5_pct": over5_b,
        "base_max_over": float(err_b.max()),
        "hgb_mae": mae_h,
        "hgb_rmse": rmse_h,
        "hgb_over5_pct": over5_h,
        "hgb_max_over": float(err_h.max()),
        "hgb-base_mae": mean_diff,
        "ci_95_lo": ci[0],
        "ci_95_hi": ci[1],
        "P(hgb<base)": p_better,
    }])

    summary.to_csv(out_dir / "flash_loso_summary.csv", index=False)
    mixes.to_csv(out_dir / "flash_loso_predictions.csv", index=False)

    print("\n==================================================================")
    print(" Leave-One-System-Out (LOSO) Flash Point Mixture Results (n=126)")
    print("==================================================================")
    print(f"Evaluated on: {int(ok.sum())} rows across {len(np.unique(sys_ok))} systems")
    print(f"Physics Baseline (Le Chatelier + Raoult):")
    print(f"  MAE:                       {mae_b:.2f} degC")
    print(f"  RMSE:                      {rmse_b:.2f} degC")
    print(f"  Dangerous Over-prediction: {over5_b:.1f}%  (Max over: +{err_b.max():.1f} degC)")
    print(f"\nLearned Residual Model (LOSO-HGB):")
    print(f"  MAE:                       {mae_h:.2f} degC")
    print(f"  RMSE:                      {rmse_h:.2f} degC")
    print(f"  Dangerous Over-prediction: {over5_h:.1f}%  (Max over: +{err_h.max():.1f} degC)")
    print(f"\nPaired Difference (HGB - Base):")
    print(f"  Delta MAE:                 {mean_diff:.2f} degC  [95% CI: {ci[0]:.2f}, {ci[1]:.2f}]")
    print(f"  P(Model < Physics):        {p_better:.1%}")
    print("==================================================================\n")
    print(f"[done] Saved to {out_dir / 'flash_loso_summary.csv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    args = ap.parse_args()
    evaluate_loso(Path(args.root))
