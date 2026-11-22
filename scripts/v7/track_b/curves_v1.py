"""v7 Task 5: epoch curves of the v1 D-MPNN configs, v1-cap vs long-cap. Dev rows only (ysi_holdout masked). See PREDECLARATION.md."""
import sys, json, time
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from fuelblend.ml import dmpnn as D
from fuelblend.ml.curves import CurveLogger
from fuelblend.v7 import runner

R = Path(__file__).resolve().parents[3]
OUT = R / "results/v7/track_b/curves_v1"
QUICK = "--quick" in sys.argv
LONG_ONLY = "--long-only" in sys.argv   # post-hoc rerun that also saves OOF predictions (Task 5 first run did not)
if QUICK: OUT = OUT.parent / "curves_v1_quick"
if LONG_ONLY: OUT = OUT.parent / "curves_v1_long"
OUT.mkdir(parents=True, exist_ok=True)
CFGS = [dict(hidden=200, depth=4, dropout=0.1), dict(hidden=300, depth=3, dropout=0.2)]
REG = {"v1cap": dict(epochs=150, patience=25), "long": dict(epochs=1000, patience=100)}
NM = 2 if QUICK else 5
if LONG_ONLY: REG = {"long": REG["long"]}
if QUICK and not LONG_ONLY: REG = {"v1cap": dict(epochs=6, patience=3), "long": dict(epochs=12, patience=6)}

comp, folds, tasks = runner.load_track_b_pool()
Yraw = comp[tasks].values.astype(float); Mask = ~np.isnan(Yraw)
feats = [D.featurize(s) for s in comp.smiles]
g = comp.grp.values
print(len(comp), "compounds", len(tasks), "tasks", D.device_info(), flush=True)
for p in (OUT / "curves.csv",):
    if p.exists(): p.unlink()
lg = CurveLogger(OUT / "curves.csv")

oof = {r: np.full(Yraw.shape, np.nan) for r in REG}; const = np.full(Yraw.shape, np.nan)
best_rows, chosen = [], []
t0 = time.time()
for k in range(5):
    te = np.where(folds == k)[0]; tr = np.where(folds != k)[0]
    ug = np.unique(g[tr]); rng = np.random.default_rng(k)
    vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
    va = np.array([i for i in tr if g[i] in vg]); tr2 = np.array([i for i in tr if g[i] not in vg])
    Ym = np.where(Mask[tr2], Yraw[tr2], np.nan)
    mu = np.nan_to_num(np.nanmean(Ym, 0)); sd = np.nanstd(Ym, 0) + 1e-9; sd = np.where(np.isnan(sd), 1, sd)
    Z = np.nan_to_num((Yraw - mu) / sd)
    for j in range(len(tasks)):
        trm = tr[Mask[tr, j]]; const[te, j] = np.median(Yraw[trm, j]) if len(trm) else np.nan
    for rn, rc in REG.items():
        best = None
        for ci, c in enumerate(CFGS):
            ms = D.fit_ensemble([feats[i] for i in tr2], Z[tr2], Mask[tr2], [feats[i] for i in va], Z[va], Mask[va],
                                n_models=NM, seed=k, logger=lg.child(model=f"{rn}_cfg{ci}", fold=k), **rc, **c)
            pv, _ = D.predict(ms, [feats[i] for i in va])
            vl = float((((pv - Z[va]) ** 2) * Mask[va]).sum() / max(Mask[va].sum(), 1))
            for s_, m in enumerate(ms): best_rows.append(dict(regime=rn, fold=k, cfg=ci, seed=s_, best_epoch=m.best_epoch_))
            if best is None or vl < best[0]: best = (vl, ms, ci)
        p, _ = D.predict(best[1], [feats[i] for i in te])
        oof[rn][te] = p * sd + mu; chosen.append(dict(regime=rn, fold=k, cfg=best[2]))
        print("fold", k, rn, "cfg", best[2], f"{time.time()-t0:.0f}s", flush=True)
lg.close()

rows = []
for j, t in enumerate(tasks):
    m = Mask[:, j]; y = Yraw[m, j]; r = dict(target=t, n=int(m.sum()), mae_const=np.abs(const[m, j] - y).mean())
    for rn in REG: r[f"mae_{rn}"] = np.abs(oof[rn][m, j] - y).mean()
    rows.append(r)
oo = comp[["inchikey", "smiles", "grp"]].copy()
for j, t in enumerate(tasks):
    oo[t + "_obs"] = Yraw[:, j]
    for rn in REG: oo[f"{t}_{rn}"] = oof[rn][:, j]
oo.to_csv(OUT / "oof.csv", index=False)
res = pd.DataFrame(rows); res.to_csv(OUT / "cv_summary.csv", index=False)
be = pd.DataFrame(best_rows); be.to_csv(OUT / "best_epoch.csv", index=False)
cv = pd.read_csv(OUT / "curves.csv")
spe = cv.groupby("model").sec.median().to_dict()
meta = dict(median_best_epoch=be.groupby("regime").best_epoch.median().to_dict(), sec_per_epoch_median=spe,
            chosen=chosen, n_models=NM, regimes=REG, n_locked_ysi=comp.attrs.get("n_locked_ysi"), total_s=time.time() - t0, **D.device_info())
(OUT / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    for mdl, d in cv.groupby("model"): ax.plot(d.groupby("epoch").val_loss.median(), label=mdl)
    ax.set_xlabel("epoch"); ax.set_ylabel("median val loss (z)"); ax.set_yscale("log"); ax.legend(fontsize=7); fig.tight_layout()
    fig.savefig(OUT / "curves.png", dpi=130)
except Exception as e: print("plot failed", e)
print(res.round(3).to_string()); print(meta["median_best_epoch"])
