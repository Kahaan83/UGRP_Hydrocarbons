"""v7 Task 10: Track B D-MPNN retrain with nested Optuna on the GPU. See PREDECLARATION.md.
Usage: FUELBLEND_DEVICE=cuda python scripts/v7/track_b/retrain.py [--quick]"""
import sys, json, time, os, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from fuelblend.ml import dmpnn as D
from fuelblend.ml.curves import CurveLogger
from fuelblend.v7 import runner, hpo

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/track_b" / ("retrain_quick" if QUICK else "retrain")
OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs"))
NT, TO, NM = (2, 300, 2) if QUICK else (60, 5400, 5)
space = hpo.load_space("dmpnn")
if QUICK:
    space = {**space, "fixed": {"epochs": 6, "patience": 3}}
comp, folds, tasks = runner.load_track_b_pool()
g = comp.grp.values
Yraw = comp[tasks].values.astype(float)
Mk = ~np.isnan(Yraw)
feats = [D.featurize(s) for s in comp.smiles]
lg = CurveLogger(OUT / "curves.csv", tb_dir=RUNS / "track_b" / ("quick" if QUICK else "retrain"))
oom = {"n": 0}


def scale(tr):
    Ym = np.where(Mk[tr], Yraw[tr], np.nan)
    mu = np.nan_to_num(np.nanmean(Ym, 0))
    sd = np.nanstd(Ym, 0) + 1e-9
    sd = np.where(np.isnan(sd) | (Mk[tr].sum(0) < 2) | (sd < 1e-6), 1.0, sd)
    return mu, sd


def _fit(tr, va, params, n_models, seed, logger=None, on_epoch=None, bs=None, oom_retry=0):
    if oom_retry:
        oom["n"] += 1
    mu, sd = scale(tr)
    Z = np.nan_to_num((Yraw - mu) / sd)
    p = {k: v for k, v in params.items() if k != "bs"}
    ms = D.fit_ensemble([feats[i] for i in tr], Z[tr], Mk[tr], [feats[i] for i in va], Z[va], Mk[va], n_models=n_models, seed=seed,
                        bs=int(bs if bs is not None else params.get("bs", 32)), logger=logger, on_epoch=on_epoch, **p)
    return ms, mu, sd, Z


def fit_score(params, tr, va, logger, on_epoch=None):
    ms, mu, sd, Z = runner.fit_with_oom_retry(lambda **kw: _fit(tr, va, params, 1, 0, on_epoch=on_epoch, **kw), bs=params.get("bs", 32))
    pv, _ = D.predict(ms, [feats[i] for i in va])
    return float((((pv - Z[va]) ** 2) * Mk[va]).sum() / max(Mk[va].sum(), 1))


def final_fit(params, tr, te):
    ug = np.unique(g[tr])
    k = int(folds[te[0]])
    vg = set(np.random.default_rng(k).choice(ug, max(1, int(.15 * len(ug))), replace=False))
    va = np.array([i for i in tr if g[i] in vg])
    tr2 = np.array([i for i in tr if g[i] not in vg])
    ms, mu, sd, Z = runner.fit_with_oom_retry(lambda **kw: _fit(tr2, va, params, NM, k, logger=lg.child(model="final", fold=k), **kw),
                                              bs=params.get("bs", 32))
    return D.predict(ms, [feats[i] for i in te])[0] * sd + mu


t0 = time.time()
res = hpo.nested_cv(comp, g, folds, fit_score=fit_score, space=space, n_trials=NT, timeout_s=TO,
                    study_prefix="v7_trackb_dmpnn" + ("_quick" if QUICK else ""), final_fit=final_fit, out_dir=OUT / "preds", logger=None)
lg.close()
res.to_csv(OUT / "best_params.csv", index=False)
oof = np.full(Yraw.shape, np.nan)
for k in range(5):
    oof[folds == k] = np.load(OUT / "preds" / f"pred_fold{k}.npy")

v6 = pd.read_csv(R / "results/v6/gnn/oof.csv").set_index("inchikey").reindex(comp.inchikey)
tr_ = pd.read_csv(R / "results/v7/track_b/descriptor_gbm/oof.csv").set_index("inchikey").reindex(comp.inchikey)
const = np.full(Yraw.shape, np.nan)
for k in range(5):
    for j in range(len(tasks)):
        trm = (folds != k) & Mk[:, j]
        const[folds == k, j] = np.median(Yraw[trm, j])
rows, ps_v1, ps_tr = [], [], []


def boot_p(d, gg):
    ug = np.unique(gg)
    S = np.array([d[gg == u].sum() for u in ug])
    N = np.array([(gg == u).sum() for u in ug], float)
    idx = np.random.default_rng(0).integers(0, len(ug), (4000, len(ug)))
    est = S[idx].sum(1) / N[idx].sum(1)
    return est, 2 * min((est <= 0).mean(), (est >= 0).mean())


for j, t in enumerate(tasks):
    y = Yraw[:, j]
    m = Mk[:, j] & ~np.isnan(oof[:, j]) & ~np.isnan(v6[t + "_oof"].values) & ~np.isnan(tr_[t + "_gbm"].values)
    e7 = np.abs(oof[m, j] - y[m]); e6 = np.abs(v6[t + "_oof"].values[m] - y[m])
    eg = np.abs(tr_[t + "_gbm"].values[m] - y[m]); er = np.abs(tr_[t + "_rf"].values[m] - y[m]); ec = np.abs(const[m, j] - y[m])
    et = eg if eg.mean() <= er.mean() else er
    est1, p1 = boot_p(e6 - e7, g[m])
    est2, p2 = boot_p(et - e7, g[m])
    ps_v1.append(p1); ps_tr.append(p2)
    rows.append(dict(target=t, n=int(m.sum()), mae_gnn=e7.mean(), mae_const=ec.mean(), skill=1 - e7.mean() / ec.mean(),
                     rmse_gnn=float(np.sqrt(((oof[m, j] - y[m]) ** 2).mean())), over_pred=float(np.maximum(oof[m, j] - y[m], 0).mean()),
                     mae_v1=e6.mean(), mae_gbm=eg.mean(), mae_rf=er.mean(), gain_vs_v1=(e6 - e7).mean(), v1_lo=np.percentile(est1, 5),
                     v1_hi=np.percentile(est1, 95), p_v1=p1, gain_vs_trees=(et - e7).mean(), tr_lo=np.percentile(est2, 5),
                     tr_hi=np.percentile(est2, 95), p_trees=p2))


def holm(ps):
    o = np.argsort(ps)
    adj = np.empty(len(ps))
    run = 0
    for r_, i in enumerate(o):
        run = max(run, min(1, (len(ps) - r_) * ps[i]))
        adj[i] = run
    return adj


S = pd.DataFrame(rows)
S["p_holm_v1"] = holm(ps_v1)
S["p_holm_trees"] = holm(ps_tr)
S["beats_v1"] = (S.p_holm_v1 < 0.10) & (S.gain_vs_v1 > 0)
S["earns_place_vs_trees"] = (S.p_holm_trees < 0.10) & (S.gain_vs_trees > 0)
oo = comp[["inchikey", "smiles", "is_hc", "grp"]].copy()
for j, t in enumerate(tasks):
    oo[t + "_obs"] = Yraw[:, j]
    oo[t + "_oof"] = oof[:, j]
runner.write_results(OUT, summary=S, oof=oo, meta=dict(task="10_trackb_retrain", n_trials=NT, timeout_s=TO, n_models=NM, oom_retries=oom["n"],
                                                     total_s=time.time() - t0, n_locked_ysi=comp.attrs.get("n_locked_ysi")))
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cv = pd.read_csv(OUT / "curves.csv")
    fig, ax = plt.subplots(figsize=(7, 4))
    for f_, d in cv.groupby("fold"):
        ax.plot(d.groupby("epoch").val_loss.median(), label=f"fold {f_}")
    ax.set_yscale("log"); ax.set_xlabel("epoch"); ax.set_ylabel("median val loss (z)"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "curves.png", dpi=130)
except Exception as e:
    print("plot", e)
print(S[["target", "n", "mae_gnn", "mae_v1", "mae_gbm", "mae_rf", "gain_vs_v1", "p_holm_v1", "gain_vs_trees", "p_holm_trees"]].round(3).to_string())
