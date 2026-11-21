"""Track B step 2: grouped 5-fold CV of the multi-task D-MPNN vs constant and existing estimators (Tm: v5 fusion_estimate,
YSI: v4 ridge GC). Group = generic Murcko scaffold (acyclics bucketed by carbon count). Config chosen per outer fold on an
inner validation split (group-held-out). Oxygenates only carry ron/mon labels (see build_pool.py)."""
import sys, time, json
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold
from fuelblend.ml import dmpnn as D
from fuelblend.ml import ysi_model as Y
from fuelblend.thermo import fusion_estimate as F

R = Path(__file__).resolve().parents[3]
V2 = "--v2" in sys.argv   # post-hoc optimisation variant (wider grid, longer patience); writes to results/v6/gnn_v2
OUT = R / ("results/v6/gnn_v2" if V2 else "results/v6/gnn")
POOL = R / "results/v6/gnn"
import os
if V2 and os.environ.get("FUELBLEND_DEVICE", "cpu") != "cpu": OUT = R / "results/v6/gnn_v2_gpu"   # keep the CPU and GPU runs apart
QUICK = "--quick" in sys.argv
CFGS = [dict(hidden=200, depth=4, dropout=0.1), dict(hidden=300, depth=3, dropout=0.2)]
if V2:
    CFGS += [dict(hidden=300, depth=3, dropout=0.2, lr=5e-4, wd=1e-4, patience=40),
             dict(hidden=300, depth=4, dropout=0.2, lr=1e-3, wd=1e-4, patience=40)]
OUT.mkdir(parents=True, exist_ok=True)
NM, EP = (2, 8) if QUICK else (5, 150)

def group(s):
    m = Chem.MolFromSmiles(s)
    if m.GetRingInfo().NumRings() == 0: return f"acyc_C{sum(a.GetAtomicNum()==6 for a in m.GetAtoms())}"
    sc = MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m))
    return Chem.MolToSmiles(sc)

pool = pd.read_csv(POOL / "pool.csv")
comp = pool.drop_duplicates("inchikey")[["inchikey", "smiles", "is_hc"]].reset_index(drop=True)
comp["grp"] = comp.smiles.map(group)
tasks = sorted(pool.target.unique())
Yraw = pd.pivot_table(pool, index="inchikey", columns="target", values="value").reindex(comp.inchikey)[tasks].values
Mask = ~np.isnan(Yraw)
feats = [D.featurize(s) for s in comp.smiles]
print(len(comp), "compounds", len(tasks), "tasks", comp.grp.nunique(), "groups", flush=True)

def fit_fold(tr, te, seed):
    # inner validation: 15% of train groups held out
    g = comp.grp.values
    ug = np.unique(g[tr]); rng = np.random.default_rng(seed); vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
    va = np.array([i for i in tr if g[i] in vg]); tr2 = np.array([i for i in tr if g[i] not in vg])
    mu = np.nanmean(np.where(Mask[tr2], Yraw[tr2], np.nan), 0); sd = np.nanstd(np.where(Mask[tr2], Yraw[tr2], np.nan), 0) + 1e-9
    mu = np.nan_to_num(mu); sd = np.where(np.isnan(sd), 1, sd)
    Z = np.nan_to_num((Yraw - mu) / sd)
    best = None
    for c in CFGS:
        ms = D.fit_ensemble([feats[i] for i in tr2], Z[tr2], Mask[tr2], [feats[i] for i in va], Z[va], Mask[va],
                            n_models=NM, epochs=EP, seed=seed, **c)
        pv, _ = D.predict(ms, [feats[i] for i in va])
        vl = float((((pv - Z[va]) ** 2) * Mask[va]).sum() / max(Mask[va].sum(), 1))
        if best is None or vl < best[0]: best = (vl, ms, c)
    p, s = D.predict(best[1], [feats[i] for i in te])
    return p * sd + mu, s * sd, best[2], mu

oof = np.full(Yraw.shape, np.nan); oofsd = np.full(Yraw.shape, np.nan); const = np.full(Yraw.shape, np.nan)
base = {t: np.full(len(comp), np.nan) for t in ("tm_K", "ysi")}
chosen = []
t0 = time.time()
for k, (tr, te) in enumerate(GroupKFold(5).split(comp, groups=comp.grp)):
    p, s, c, mu = fit_fold(tr, te, k)
    oof[te], oofsd[te] = p, s
    for j, t in enumerate(tasks):
        trm = tr[Mask[tr, j]]
        const[te, j] = np.median(Yraw[trm, j]) if len(trm) else np.nan
    # existing estimators
    j = tasks.index("tm_K"); trm = tr[Mask[tr, j]]
    try: base["tm_K"][te] = np.asarray(F.fit_predict(F.features(comp.smiles[trm].tolist()), Yraw[trm, j], F.features(comp.smiles[te].tolist())))
    except Exception as e: print("tm baseline failed", e)
    j = tasks.index("ysi"); trm = tr[Mask[tr, j]]
    try:
        Xtr = np.vstack([Y.features(s_) for s_ in comp.smiles[trm]]); m = Y.fit_ridge(Xtr, Yraw[trm, j], comp.grp.values[trm])
        base["ysi"][te] = m.predict(np.vstack([Y.features(s_) for s_ in comp.smiles[te]]))
    except Exception as e: print("ysi baseline failed", e)
    chosen.append(c); print("fold", k, c, f"{time.time()-t0:.0f}s", flush=True)

rows = []
rng = np.random.default_rng(0); ug = comp.grp.unique()
for j, t in enumerate(tasks):
    m = Mask[:, j] & ~np.isnan(oof[:, j]); y = Yraw[m, j]; e = np.abs(oof[m, j] - y); ec = np.abs(const[m, j] - y)
    r = dict(target=t, n=int(m.sum()), mae_gnn=e.mean(), mae_const=ec.mean(), skill=1 - e.mean() / ec.mean(),
             rmse_gnn=np.sqrt(((oof[m, j] - y) ** 2).mean()), over_pred=np.maximum(oof[m, j] - y, 0).mean())
    gg = comp.grp.values[m]
    d = ec - e; bs = []
    idx = {g_: np.where(gg == g_)[0] for g_ in np.unique(gg)}; keys = list(idx)
    for _ in range(500):
        pick = rng.choice(len(keys), len(keys)); ii = np.concatenate([idx[keys[q]] for q in pick]); bs.append(d[ii].mean())
    r["gain_lo"], r["gain_hi"] = np.percentile(bs, [5, 95])
    if t in base:
        b = base[t][m]; ok = ~np.isnan(b); r["mae_existing"] = np.abs(b[ok] - y[ok]).mean(); r["mae_gnn_same"] = e[ok].mean()
        r["n_existing"] = int(ok.sum())
    rows.append(r)
res = pd.DataFrame(rows); res.to_csv(OUT / ("cv_summary_quick.csv" if QUICK else "cv_summary.csv"), index=False)
oo = comp[["inchikey", "smiles", "is_hc", "grp"]].copy()
for j, t in enumerate(tasks): oo[t + "_obs"] = Yraw[:, j]; oo[t + "_oof"] = oof[:, j]; oo[t + "_sd"] = oofsd[:, j]
oo.to_csv(OUT / ("oof_quick.csv" if QUICK else "oof.csv"), index=False)
json.dump(chosen, open(OUT / ("chosen_quick.json" if QUICK else "chosen_cfg.json"), "w"))
print(res.round(3).to_string())
