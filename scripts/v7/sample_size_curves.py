"""v7 Task 8: sample-size learning curves. Usage: sample_size_curves.py <gbm_trackb|dmpnn_trackb|freezing> [--quick]"""
import sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/v6/freezing")
from pathlib import Path
import numpy as np, pandas as pd
from fuelblend.v7 import runner
from fuelblend.v7.sample_size import group_subsample, slope_ci, FRACS

R = Path(__file__).resolve().parents[2]
WHICH = sys.argv[1]; QUICK = "--quick" in sys.argv
FR = (0.4, 1.0) if QUICK else FRACS
REPS = 1 if QUICK else 3
rows = []


def subsample(grp_tr, f, rep):
    return group_subsample(grp_tr, f, np.random.default_rng(1000 + rep))


if WHICH in ("gbm_trackb", "dmpnn_trackb"):
    comp, folds, tasks = runner.load_track_b_pool()
    g = comp.grp.values
    Yraw = comp[tasks].values.astype(float)
    Mk = ~np.isnan(Yraw)
    if WHICH == "gbm_trackb":
        import lightgbm as lgb
        from fuelblend.features.rdkit_desc import rdkit_2d
        X = rdkit_2d(comp.smiles.tolist())
    else:
        from fuelblend.ml import dmpnn as D
        feats = [D.featurize(s) for s in comp.smiles]
        NM, EP, PAT = (1, 8, 3) if QUICK else (2, 300, 50)
    t0 = time.time()
    for k in range(5):
        te = np.where(folds == k)[0]
        tr_all = np.where(folds != k)[0]
        for rep in range(REPS):
            for f in FR:
                sub = tr_all[subsample(g[tr_all], f, rep)]
                if WHICH == "gbm_trackb":
                    for j, t in enumerate(tasks):
                        tr = sub[Mk[sub, j]]
                        tt = te[Mk[te, j]]
                        if len(tr) < 5 or len(tt) == 0:
                            continue
                        m = lgb.LGBMRegressor(num_leaves=15, learning_rate=0.05, n_estimators=300, min_child_samples=5, subsample=0.8,
                                              subsample_freq=1, colsample_bytree=0.8, verbose=-1, n_jobs=4, random_state=rep).fit(X.iloc[tr], Yraw[tr, j])
                        rows.append(dict(frac=f, repeat=rep, fold=k, target=t, n_train=len(tr), n_train_groups=len(np.unique(g[tr])),
                                         mae=float(np.abs(m.predict(X.iloc[tt]) - Yraw[tt, j]).mean())))
                else:
                    ug = np.unique(g[sub])
                    rg = np.random.default_rng(k * 7 + rep)
                    vg = set(rg.choice(ug, max(1, int(.15 * len(ug))), replace=False))
                    va = np.array([i for i in sub if g[i] in vg])
                    tr2 = np.array([i for i in sub if g[i] not in vg])
                    if len(tr2) < 4 or len(va) < 1:
                        continue
                    Ym = np.where(Mk[tr2], Yraw[tr2], np.nan)
                    mu = np.nan_to_num(np.nanmean(Ym, 0))
                    sd = np.nanstd(Ym, 0) + 1e-9
                    sd = np.where(np.isnan(sd), 1, sd)
                    Z = np.nan_to_num((Yraw - mu) / sd)
                    ms = D.fit_ensemble([feats[i] for i in tr2], Z[tr2], Mk[tr2], [feats[i] for i in va], Z[va], Mk[va], n_models=NM,
                                        epochs=EP, patience=PAT, seed=k * 10 + rep, hidden=200, depth=4, dropout=0.1)
                    p = D.predict(ms, [feats[i] for i in te])[0] * sd + mu
                    for j, t in enumerate(tasks):
                        mm = Mk[te, j]
                        if mm.sum() == 0 or Mk[tr2, j].sum() < 1:
                            continue
                        rows.append(dict(frac=f, repeat=rep, fold=k, target=t, n_train=int(Mk[tr2, j].sum()),
                                         n_train_groups=len(np.unique(g[tr2])), mae=float(np.abs(p[mm, j] - Yraw[te, j][mm]).mean())))
                print(WHICH, "fold", k, "rep", rep, "frac", f, f"{time.time() - t0:.0f}s", flush=True)
    OUTS = {WHICH: pd.DataFrame(rows)}
else:   # freezing LOGO, two models
    import step3_models as S
    ds = pd.read_csv(R / "results/v6/freezing/features_preds.csv", dtype={"fuel_id": str})
    dev = ds[ds.role == "dev"].reset_index(drop=True)
    y = dev.freezing_C.to_numpy()
    fits = {"freezing_ridge": S.fp_ridge, "freezing_lgbm": S.make_lgbm("pred_e")}
    out = {k: [] for k in fits}
    for i in range(len(dev)):
        tr_all = np.array([j for j in range(len(dev)) if j != i])
        grp = dev.entity.values[tr_all]
        for rep in range(REPS):
            for f in FR:
                sub = tr_all[subsample(grp, f, rep)]
                for nm, fn in fits.items():
                    try:
                        p = fn(dev.iloc[sub], dev.iloc[[i]])[0]
                    except Exception:
                        continue
                    out[nm].append(dict(frac=f, repeat=rep, fold=i, target="freezing_C", n_train=len(sub), n_train_groups=len(sub),
                                        mae=float(abs(p - y[i]))))
    OUTS = {k: pd.DataFrame(v) for k, v in out.items()}

for name, df in OUTS.items():
    od = R / "results/v7/sample_size" / (name + ("_quick" if QUICK else ""))
    od.mkdir(parents=True, exist_ok=True)
    df.to_csv(od / "curve.csv", index=False)
    top = sorted(df.frac.unique())[-2:]
    sl = []
    for t, d in df.groupby("target"):
        piv = d[d.frac.isin(top)].groupby(["repeat", "fold", "frac"]).mae.mean().unstack("frac").dropna()
        if len(piv) < 3:
            continue
        n = d[d.frac.isin(top)].groupby("frac").n_train.mean().loc[top].values
        lo, hi, pt = slope_ci(n, piv[top].values)
        sl.append(dict(target=t, frac_lo=top[0], frac_hi=top[1], n_lo=n[0], n_hi=n[1], mae_lo=piv[top[0]].mean(), mae_hi=piv[top[1]].mean(),
                       slope=pt, slope_lo=lo, slope_hi=hi, flattened=bool(lo <= 0 <= hi)))
    pd.DataFrame(sl).to_csv(od / "slopes.csv", index=False)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        for t, d in df.groupby("target"):
            m = d.groupby("frac").mae.mean()
            ax.plot(m.index, m.values / m.values[-1], marker="o", label=t)
        ax.set_xlabel("fraction of training groups"); ax.set_ylabel("MAE / MAE(frac=1)"); ax.legend(fontsize=6)
        fig.tight_layout(); fig.savefig(od / "curve.png", dpi=120)
    except Exception as e:
        print("plot", e)
    print(name)
    print(pd.DataFrame(sl).round(3).to_string())
