"""
eval_grouped_cv.py — 5-Fold Grouped Cross-Validation across All 510 Mixture Systems.

Replaces the single S2 split with rigorous 5-fold CV grouped by system_id.
Every mixture observation is evaluated as an out-of-fold test sample.
Computes:
  1. Viscosity (Refutas physics baseline vs HGB residual)
  2. Density (Ideal-volume physics baseline vs HGB residual)
  3. Flash Point (Le Chatelier physics baseline vs HGB residual)
  4. 4,000-resample cluster bootstrap across all 510 systems for 95% confidence intervals.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider, RegistryTensors
from fuelblend.mixture.t8_deepsets import build_data


def cluster_bootstrap_metric(err_a, err_b, sysid, is_rmse=True, B=4000, seed=42):
    """Cluster bootstrap of paired difference metric(a) - metric(b). Negative = a better."""
    f = np.square if is_rmse else np.abs
    u, inv = np.unique(sysid, return_inverse=True)
    n_s = np.bincount(inv).astype(float)
    sa = np.bincount(inv, weights=f(err_a))
    sb = np.bincount(inv, weights=f(err_b))
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(len(u), np.full(len(u), 1.0 / len(u)), size=B)
    n = counts @ n_s
    A = counts @ sa / n
    B_m = counts @ sb / n
    if is_rmse:
        A, B_m = np.sqrt(A), np.sqrt(B_m)
    d = A - B_m
    return float(d.mean()), np.percentile(d, [2.5, 97.5]), float((d < 0).mean())


def run_grouped_cv(root: Path, n_splits=5, seed=42):
    import torch
    dev = torch.device("cpu")
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    rt = RegistryTensors(reg, dev)
    prov = PurePropertyProvider(reg, root / "data/cache/pure_grid.npz", verbose=False)
    data = build_data(root, reg, prov, rt, dev)

    meta = data.meta
    is_mix = meta["n_comp"].to_numpy() >= 2
    lab = data.lab
    y = data.y
    base = data.base_raw
    sys_ids = meta["system_id"].to_numpy()
    pidx = meta["prop"].to_numpy()

    # Mixture-level feature extraction for tabular HGB baseline
    idx_arr, w_arr = data.idx, data.w
    safe_idx = np.where(idx_arr >= 0, idx_arr, 0)
    Sw = reg.S[safe_idx]
    wm = w_arr[..., None]
    feat_mean = (wm * Sw).sum(1)
    feat_std = np.sqrt((wm * (Sw - feat_mean[:, None, :]) ** 2).sum(1))
    ln_p, rho_p = np.nan_to_num(data.ln_eta_i), np.nan_to_num(data.rho_i, nan=750.0)
    ln_m, rho_m = (w_arr * ln_p).sum(1), (w_arr * rho_p).sum(1)
    T = meta["T_K"].to_numpy(float)
    X = np.column_stack([feat_mean, feat_std, base, T, data.ctx, meta["n_comp"].to_numpy(float),
                         w_arr.max(1), ln_m, rho_m])

    # Out-of-fold prediction array
    oof_hgb = np.full((len(y), 3), np.nan)

    # GroupKFold over mixture systems
    unique_mix_sys = np.unique(sys_ids[is_mix])
    print(f"[Grouped-CV] Setting up {n_splits}-fold CV over {len(unique_mix_sys)} unique mixture systems...")

    gkf = GroupKFold(n_splits=n_splits)
    mix_indices = np.where(is_mix)[0]
    mix_sys = sys_ids[mix_indices]

    for fold, (train_s_idx, test_s_idx) in enumerate(gkf.split(mix_indices, groups=mix_sys), 1):
        test_rows = mix_indices[test_s_idx]
        train_rows = np.setdiff1d(np.arange(len(y)), test_rows)  # includes pures in training

        print(f"  Fold {fold}/{n_splits}: {len(test_rows)} test mixture rows, {len(train_rows)} train rows...")
        for p in range(3):
            # Train mask for property p
            p_tr = train_rows[lab[train_rows, p]]
            if len(p_tr) < 20:
                continue
            X_tr = X[p_tr].copy()
            X_te = X[test_rows].copy()
            # Handle all-nan columns (e.g. ambient T for flash)
            all_nan = np.isnan(X_tr).all(axis=0)
            X_tr[:, all_nan] = 0.0
            X_te[:, all_nan] = 0.0

            target = y[p_tr, p] - base[p_tr, p]
            m = HistGradientBoostingRegressor(max_iter=250, learning_rate=0.06, max_leaf_nodes=12,
                                              min_samples_leaf=5, l2_regularization=1.0, random_state=seed + fold)
            m.fit(X_tr, target)
            oof_hgb[test_rows, p] = base[test_rows, p] + m.predict(X_te)

    # Evaluation summary
    out_dir = root / "results/eval"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    prop_names = ["viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"]
    short_names = ["visc", "dens", "flash"]

    print("\n=========================================================================================")
    print(f" 5-Fold Grouped Cross-Validation Results Across All 510 Mixture Systems")
    print("=========================================================================================")

    for p in range(3):
        mask = is_mix & lab[:, p] & np.isfinite(oof_hgb[:, p]) & (pidx == p)
        if mask.sum() == 0:
            continue
        r = np.where(mask)[0]
        y_true = y[r, p]
        y_b = base[r, p]
        y_h = oof_hgb[r, p]
        sids = sys_ids[r]

        err_b = y_b - y_true
        err_h = y_h - y_true

        is_rmse = p != 2
        rmse_b = float(np.sqrt(np.mean(err_b ** 2)))
        rmse_h = float(np.sqrt(np.mean(err_h ** 2)))
        mae_b = float(np.mean(np.abs(err_b)))
        mae_h = float(np.mean(np.abs(err_h)))

        diff_mean, ci, p_better = cluster_bootstrap_metric(err_h, err_b, sids, is_rmse=is_rmse)

        row = {
            "ci_method": "4000-resample system-level cluster bootstrap (system_id)",
            "property": short_names[p],
            "n_obs": int(len(r)),
            "n_systems": int(len(np.unique(sids))),
            "base_rmse": rmse_b,
            "base_mae": mae_b,
            "model_rmse": rmse_h,
            "model_mae": mae_h,
            "diff_metric": diff_mean,
            "ci_95_lo": ci[0],
            "ci_95_hi": ci[1],
            "P(model<base)": p_better,
        }
        if p == 0:
            row["base_mape_pct"] = float(np.mean(np.abs(np.exp(err_b) - 1.0)) * 100.0)
            row["model_mape_pct"] = float(np.mean(np.abs(np.exp(err_h) - 1.0)) * 100.0)
            print(f"Kinematic Viscosity (n={len(r)} rows across {row['n_systems']} systems):")
            print(f"  Physics Baseline (Refutas): log-RMSE = {rmse_b:.4f} (MAPE = {row['base_mape_pct']:.2f}%)")
            print(f"  Learned Model (5-Fold CV):  log-RMSE = {rmse_h:.4f} (MAPE = {row['model_mape_pct']:.2f}%)")
            print(f"  Bootstrap Delta log-RMSE:   {diff_mean:.4f}  [95% CI: {ci[0]:.4f}, {ci[1]:.4f}]  P={p_better:.1%}")
        elif p == 1:
            print(f"\nLiquid Density (n={len(r)} rows across {row['n_systems']} systems):")
            print(f"  Physics Baseline (Ideal Vol): RMSE = {rmse_b:.4f} kg/m3 (MAE = {mae_b:.4f} kg/m3)")
            print(f"  Learned Model (5-Fold CV):    RMSE = {rmse_h:.4f} kg/m3 (MAE = {mae_h:.4f} kg/m3)")
            print(f"  Bootstrap Delta RMSE:         {diff_mean:.4f} kg/m3  [95% CI: {ci[0]:.4f}, {ci[1]:.4f}]  P={p_better:.1%}")
        else:
            over5_b = float((err_b > 5.0).mean() * 100.0)
            over5_h = float((err_h > 5.0).mean() * 100.0)
            row["base_over5_pct"] = over5_b
            row["model_over5_pct"] = over5_h
            print(f"\nFlash Point (n={len(r)} rows across {row['n_systems']} systems):")
            print(f"  Physics Baseline (Le Chatelier): MAE = {mae_b:.2f} degC (Over > 5C = {over5_b:.1f}%)")
            print(f"  Learned Model (5-Fold CV):       MAE = {mae_h:.2f} degC (Over > 5C = {over5_h:.1f}%)")
            print(f"  Bootstrap Delta MAE:             {diff_mean:.2f} degC  [95% CI: {ci[0]:.2f}, {ci[1]:.2f}]  P={p_better:.1%}")
        summary_rows.append(row)

    print("=========================================================================================\n")
    df_res = pd.DataFrame(summary_rows)
    df_res.to_csv(out_dir / "grouped_cv_summary.csv", index=False)
    print(f"[done] Saved 5-fold grouped CV summary to {out_dir / 'grouped_cv_summary.csv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--splits", type=int, default=5)
    args = ap.parse_args()
    run_grouped_cv(Path(args.root), n_splits=args.splits)
