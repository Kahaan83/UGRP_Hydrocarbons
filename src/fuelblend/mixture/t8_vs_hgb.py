"""
T8 vs HistGradientBoosting residual model (WP5-style) - paired comparison on identical rows
===========================================================================================
Install at fuelblend_data/src/fuelblend/mixture/t8_vs_hgb.py

    python -m fuelblend.mixture.t8_vs_hgb --root . --pred results/t8/<tag>_predictions.csv --device cpu

What it does
  1. Rebuilds the same MixData as T8 (same rows, same splits, same physics baselines).
  2. Trains one HGB residual model per property on the train splits (y - baseline ~ mixture features),
     predicts every held-out row.
  3. Aligns with the DeepSets predictions CSV (asserts y_true matches row by row).
  4. Reports baseline / HGB / DeepSets error on mixtures (n_comp >= 2) and a CLUSTER bootstrap over system_id
     of the paired metric difference (negative = first model better), so correlated rows of one system
     do not inflate significance.
Metric: ln-viscosity RMSE, density RMSE, flash MAE (the headline metrics used elsewhere).

NOTE: this HGB uses mixture-level aggregates of the registry descriptors, not your original WP5 features.
For an exact WP5 comparison, export your WP5 held-out predictions in the same row order and replace
`hgb_predict` (the rest of the script is model-agnostic).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.ensemble import HistGradientBoostingRegressor

from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider, RegistryTensors
from fuelblend.mixture.t8_deepsets import SHORT, build_data


def mixture_features(data, reg):
    idx, w = data.idx, data.w
    safe = np.where(idx >= 0, idx, 0)
    Sw = reg.S[safe]                                        # (R,N,D)
    wm = w[..., None]
    mean = (wm * Sw).sum(1)
    std = np.sqrt((wm * (Sw - mean[:, None, :]) ** 2).sum(1))
    ln, rho = np.nan_to_num(data.ln_eta_i), np.nan_to_num(data.rho_i, nan=750.0)
    ln_m, rho_m = (w * ln).sum(1), (w * rho).sum(1)
    ln_s = np.sqrt((w * (ln - ln_m[:, None]) ** 2).sum(1))
    rho_s = np.sqrt((w * (rho - rho_m[:, None]) ** 2).sum(1))
    ent = -(w * np.log(np.clip(w, 1e-12, None))).sum(1)
    T = data.meta["T_K"].to_numpy(float)                    # NaN for flash rows: HGB handles NaN natively
    return np.column_stack([mean, std, data.base_raw, T, data.ctx, data.meta["n_comp"].to_numpy(float),
                            w.max(1), ent, ln_m, ln_s, rho_m, rho_s])


def hgb_predict(X, data, train_rows_by_prop, pred_rows, seed=0):
    out = np.full((len(data.y), 3), np.nan)
    for p, rows in train_rows_by_prop.items():
        if len(rows) < 20:
            continue
        X_tr = X[rows].copy()
        X_pr = X[pred_rows].copy()
        all_nan = np.isnan(X_tr).all(axis=0)
        X_tr[:, all_nan] = 0.0
        X_pr[:, all_nan] = 0.0
        m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=8,
                                          min_samples_leaf=max(3, min(20, len(rows) // 10)), l2_regularization=1.0,
                                          early_stopping=False, random_state=seed)
        m.fit(X_tr, data.y[rows, p] - data.base_raw[rows, p])
        out[pred_rows, p] = data.base_raw[pred_rows, p] + m.predict(X_pr)
    return out


def cluster_bootstrap(err_a, err_b, sysid, p, B=4000, seed=0):
    """paired difference metric(a) - metric(b); negative = a better. Resamples whole systems."""
    f = np.abs if p == 2 else np.square
    u, inv = np.unique(sysid, return_inverse=True)
    n_s = np.bincount(inv).astype(float)
    sa, sb = np.bincount(inv, weights=f(err_a)), np.bincount(inv, weights=f(err_b))
    counts = np.random.default_rng(seed).multinomial(len(u), np.full(len(u), 1.0 / len(u)), size=B)
    n = counts @ n_s
    A, Bm = counts @ sa / n, counts @ sb / n
    if p != 2:
        A, Bm = np.sqrt(A), np.sqrt(Bm)
    d = A - Bm
    return float(d.mean()), np.percentile(d, [2.5, 97.5]), float((d < 0).mean())


def metric(p, e):
    return float(np.abs(e).mean()) if p == 2 else float(np.sqrt(np.mean(np.square(e))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--pred", required=True, help="T8 predictions CSV")
    ap.add_argument("--train-splits", nargs="+", default=["S1"])
    ap.add_argument("--out", default="results/t8/hgb_comparison.csv")
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()

    root, dev = Path(a.root), torch.device(a.device)
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    rt = RegistryTensors(reg, dev)
    data = build_data(root, reg, PurePropertyProvider(reg, root / "data/cache/pure_grid.npz"), rt, dev)
    meta, lab = data.meta, data.lab
    split = meta["split"].to_numpy()
    tr = [s.upper() for s in a.train_splits]
    X = mixture_features(data, reg)
    held = np.where(~np.isin(split, tr) & (split != "NONE"))[0]
    hgb = hgb_predict(X, data, {p: np.where(np.isin(split, tr) & lab[:, p])[0] for p in range(3)}, held)

    ds = pd.read_csv(a.pred)
    rows_out = []
    for sp in sorted(set(split[held])):
        rows = np.where(split == sp)[0]
        d = ds[ds["split"] == sp].reset_index(drop=True)
        pidx = meta["prop"].to_numpy()[rows]
        if len(d) != len(rows) or not np.allclose(d["y_true"].to_numpy(), data.y[rows, pidx], equal_nan=True):
            raise ValueError(f"Row alignment failed for {sp}: predictions CSV does not match rebuilt data "
                             f"(len {len(d)} vs {len(rows)}). Rebuild with the same code/splits.")
        for p in range(3):
            m = (pidx == p) & lab[rows, p] & (meta["n_comp"].to_numpy()[rows] >= 2) & np.isfinite(hgb[rows, p])
            if m.sum() < 5:
                continue
            r = rows[m]
            y = data.y[r, p]
            e = {"base": data.base_raw[r, p] - y, "hgb": hgb[r, p] - y, "ds": d["deepsets"].to_numpy()[m] - y}
            sid = meta["system_id"].to_numpy()[r]
            row = {"split": sp, "prop": SHORT[p], "n": int(m.sum()), "systems": int(len(np.unique(sid))),
                   **{f"{k}_{'mae' if p == 2 else 'rmse'}": metric(p, v) for k, v in e.items()}}
            for a_, b_ in (("ds", "hgb"), ("ds", "base"), ("hgb", "base")):
                mean, ci, p_better = cluster_bootstrap(e[a_], e[b_], sid, p)
                row.update({f"{a_}-{b_}": mean, f"{a_}-{b_}_lo": ci[0], f"{a_}-{b_}_hi": ci[1],
                            f"P({a_}<{b_})": p_better})
            rows_out.append(row)
    res = pd.DataFrame(rows_out)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out, index=False)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.float_format", "{:.4f}".format)
    print(res.to_string(index=False))
    print("\nNegative differences mean the first model is better. If a 95% interval contains 0 the difference "
          "is not distinguishable from noise at this number of held-out systems.")


if __name__ == "__main__":
    main()
