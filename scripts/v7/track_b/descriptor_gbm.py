"""v7 Task 6: RDKit-descriptor LightGBM and RF comparators for Track B (CPU). See PREDECLARATION.md."""
import sys, json
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd, optuna, lightgbm as lgb
from sklearn.ensemble import RandomForestRegressor
from fuelblend.features.rdkit_desc import rdkit_2d
from fuelblend.v7 import runner

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/track_b" / ("descriptor_gbm_quick" if QUICK else "descriptor_gbm")
OUT.mkdir(parents=True, exist_ok=True)
NT = 3 if QUICK else 30
optuna.logging.set_verbosity(optuna.logging.WARNING)

comp, folds, tasks = runner.load_track_b_pool()
X = rdkit_2d(comp.smiles.tolist()); g = comp.grp.values
pred = {m: np.full((len(comp), len(tasks)), np.nan) for m in ("gbm", "rf")}

def lgb_fit(Xtr, ytr, Xva, yva, p):
    m = lgb.LGBMRegressor(num_leaves=p["num_leaves"], min_child_samples=p["mcs"], learning_rate=p["lr"], n_estimators=p["n"],
                          verbose=-1, n_jobs=4, random_state=0)
    m.fit(Xtr, ytr, eval_set=[(Xva, yva)], callbacks=[lgb.early_stopping(30, verbose=False)])
    return m

for j, t in enumerate(tasks):
    y = comp[t].values; has = ~np.isnan(y)
    for k in range(5):
        tr = np.where((folds != k) & has)[0]; te = np.where(folds == k)[0]
        if len(tr) < 20: continue
        ug = np.unique(g[tr]); rng = np.random.default_rng(k)
        vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
        va = np.array([i for i in tr if g[i] in vg]); tr2 = np.array([i for i in tr if g[i] not in vg])
        keep = X.iloc[tr].nunique() > 1; cols = X.columns[keep.values]
        Xa = X[cols]
        def obj(tri):
            p = dict(num_leaves=tri.suggest_int("num_leaves", 4, 64, log=True), mcs=tri.suggest_int("mcs", 3, 40),
                     lr=tri.suggest_float("lr", 0.01, 0.2, log=True), n=tri.suggest_int("n", 100, 800))
            m = lgb_fit(Xa.iloc[tr2], y[tr2], Xa.iloc[va], y[va], p)
            return float(np.abs(m.predict(Xa.iloc[va]) - y[va]).mean())
        st = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(seed=k)); st.optimize(obj, n_trials=NT)
        bp = st.best_params; bp = dict(num_leaves=bp["num_leaves"], mcs=bp["mcs"], lr=bp["lr"], n=bp["n"])
        m = lgb_fit(Xa.iloc[tr2], y[tr2], Xa.iloc[va], y[va], bp)   # train on tr2, early stop on va, as tuned
        pred["gbm"][te, j] = m.predict(Xa.iloc[te])
        med = Xa.iloc[tr].median()
        rf = RandomForestRegressor(n_estimators=100 if QUICK else 500, n_jobs=4, random_state=k).fit(Xa.iloc[tr].fillna(med), y[tr])
        pred["rf"][te, j] = rf.predict(Xa.iloc[te].fillna(med))
    print(t, "done", flush=True)

oof6 = pd.read_csv(R / "results/v6/gnn/oof.csv").set_index("inchikey").reindex(comp.inchikey)
rows, ps = [], []
for j, t in enumerate(tasks):
    y = comp[t].values; m = ~np.isnan(y) & ~np.isnan(pred["gbm"][:, j]) & ~np.isnan(oof6[t + "_oof"].values)
    eg, er, en = (np.abs(pred["gbm"][m, j] - y[m]), np.abs(pred["rf"][m, j] - y[m]), np.abs(oof6[t + "_oof"].values[m] - y[m]))
    best = eg if eg.mean() <= er.mean() else er
    d = best - en; gg = g[m]; ug = np.unique(gg)
    S = np.array([d[gg == u].sum() for u in ug]); N = np.array([(gg == u).sum() for u in ug], float)
    idx = np.random.default_rng(0).integers(0, len(ug), (4000, len(ug))); est = S[idx].sum(1) / N[idx].sum(1)
    p = 2 * min((est <= 0).mean(), (est >= 0).mean())
    rows.append(dict(target=t, n=int(m.sum()), mae_gnn=en.mean(), mae_gbm=eg.mean(), mae_rf=er.mean(), gain=d.mean(),
                     gain_gnn_vs_gbm_lo=np.percentile(est, 5), gain_gnn_vs_gbm_hi=np.percentile(est, 95), p_boot=p)); ps.append(p)
res = pd.DataFrame(rows)
o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
for r_, i in enumerate(o): run = max(run, min(1, (len(ps) - r_) * ps[i])); adj[i] = run
res["p_holm"] = adj; res["gnn_earns_place"] = (res.p_holm < 0.10) & (res.gain > 0)
res.to_csv(OUT / "cv_summary.csv", index=False)
oo = comp[["inchikey", "smiles", "grp"]].copy()
for j, t in enumerate(tasks): oo[t + "_obs"] = comp[t].values; oo[t + "_gbm"] = pred["gbm"][:, j]; oo[t + "_rf"] = pred["rf"][:, j]
oo.to_csv(OUT / "oof.csv", index=False)
(OUT / "meta.json").write_text(json.dumps(dict(n_trials=NT, n_locked_ysi=comp.attrs.get("n_locked_ysi")), default=str))
print(res.round(3).to_string())
