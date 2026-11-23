"""v11c: nested tuning of the random forest used for the pure-compound models (docs/v11c_predeclaration.md).
Same outer scaffold folds and inner group hold-out as v7 descriptor_gbm.py; 20 Optuna trials over max_features, min_samples_leaf, max_depth.
Compared with the saved v7 RF out-of-fold predictions; Holm over 10 properties; adopt if gain CI above 0 and Holm p < 0.05."""
import sys, json
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd, optuna
from sklearn.ensemble import RandomForestRegressor
from fuelblend.features.rdkit_desc import rdkit_2d
from fuelblend.v7 import runner
R = Path(__file__).resolve().parents[2]
OUT = R / "results/v11/rf_pure"; OUT.mkdir(parents=True, exist_ok=True)
optuna.logging.set_verbosity(optuna.logging.WARNING)
comp, folds, tasks = runner.load_track_b_pool()
X = rdkit_2d(comp.smiles.tolist()); g = comp.grp.values
base = pd.read_csv(R / "results/v7/track_b/descriptor_gbm/oof.csv")
assert (base.inchikey.values == comp.inchikey.values).all()
NJ, NT = 6, 20
def rf(p, Xa, y, seed):
    return RandomForestRegressor(n_estimators=300, max_features=p["mf"], min_samples_leaf=p["msl"], max_depth=p["md"], n_jobs=NJ, random_state=seed).fit(Xa, y)
rows, oof = [], comp[["inchikey", "grp"]].copy()
for t in tasks:
    y = comp[t].values; has = ~np.isnan(y); pr = np.full(len(comp), np.nan)
    for k in range(5):
        tr = np.where((folds != k) & has)[0]; te = np.where(folds == k)[0]
        if len(tr) < 20: continue
        ug = np.unique(g[tr]); rng = np.random.default_rng(k)
        vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
        va = np.array([i for i in tr if g[i] in vg]); tr2 = np.array([i for i in tr if g[i] not in vg])
        keep = X.iloc[tr].nunique() > 1; Xa = X[X.columns[keep.values]]
        med = Xa.iloc[tr].median(); Xf = Xa.fillna(med)
        def obj(tri):
            p = dict(mf=tri.suggest_float("mf", 0.1, 1.0), msl=tri.suggest_int("msl", 1, 8), md=tri.suggest_categorical("md", [None, 8, 14, 22]))
            return float(np.abs(rf(p, Xf.iloc[tr2], y[tr2], k).predict(Xf.iloc[va]) - y[va]).mean())
        st = optuna.create_study(direction="minimize", sampler=optuna.samplers.TPESampler(seed=k)); st.optimize(obj, n_trials=NT)
        bp = st.best_params
        pr[te] = rf(dict(mf=bp["mf"], msl=bp["msl"], md=bp["md"]), Xf.iloc[tr], y[tr], k).predict(Xf.iloc[te])
    oof[t + "_tuned"] = pr
    m = has & ~np.isnan(pr) & base[t + "_rf"].notna().values
    d = np.abs(base[t + "_rf"].values[m] - y[m]) - np.abs(pr[m] - y[m]); gg = g[m]; u = np.unique(gg)
    S = np.array([d[gg == a].sum() for a in u]); N = np.array([(gg == a).sum() for a in u], float)
    idx = np.random.default_rng(42).integers(0, len(u), (2000, len(u))); est = S[idx].sum(1) / N[idx].sum(1)
    rows.append(dict(property=t, n=int(m.sum()), mae_rf_v7=float(np.abs(base[t + "_rf"].values[m] - y[m]).mean()), mae_rf_tuned=float(np.abs(pr[m] - y[m]).mean()),
                     gain=float(d.mean()), lo=float(np.percentile(est, 2.5)), hi=float(np.percentile(est, 97.5)), p=float((np.sum(est <= 0) + 1) / 2001)))
    print(rows[-1], flush=True)
res = pd.DataFrame(rows); o = res.p.rank(method="first").astype(int)
res["p_holm"] = [min(1.0, p * (len(res) - r + 1)) for p, r in zip(res.p, o)]
res["adopt"] = (res.lo > 0) & (res.p_holm < 0.05)
res.to_csv(OUT / "summary.csv", index=False); oof.to_csv(OUT / "oof.csv", index=False)
print(res.round(4).to_string())
