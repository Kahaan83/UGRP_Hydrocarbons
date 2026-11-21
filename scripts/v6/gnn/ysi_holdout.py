"""Track B YSI locked holdout (declared in results/v6/gnn/ysi_holdout/README.md).
  python scripts/v6/gnn/ysi_holdout.py --stage lock    write lock.json (no model is fitted)
  python scripts/v6/gnn/ysi_holdout.py --stage score   fit on non-holdout compounds, score the holdout ONCE
"""
import sys, json, hashlib, math, argparse
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
sys.path.insert(0, "src")
R = Path(__file__).resolve().parents[3]
OUT = R / "results/v6/gnn/ysi_holdout"
ap = argparse.ArgumentParser(); ap.add_argument("--stage", required=True, choices=["lock", "score"]); args = ap.parse_args()
ROW_MAX = 0.25


def group(s):
    m = Chem.MolFromSmiles(s)
    if m.GetRingInfo().NumRings() == 0: return f"acyc_C{sum(a.GetAtomicNum()==6 for a in m.GetAtoms())}"
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))


pool = pd.read_csv(R / "results/v6/gnn/pool.csv")
comp = pool.drop_duplicates("inchikey")[["inchikey", "smiles", "is_hc"]].reset_index(drop=True)
comp["grp"] = comp.smiles.map(group)
ysi = pool[pool.target == "ysi"].set_index("inchikey").value
yc = comp[comp.inchikey.isin(ysi.index)]

if args.stage == "lock":
    cnt = yc.groupby("grp").size()
    order = sorted(cnt.index, key=lambda g: hashlib.md5(f"v6ysi:{g}".encode()).hexdigest())
    locked, n = [], 0
    for g in order:
        if n + cnt[g] <= ROW_MAX * len(yc):
            locked.append(g); n += int(cnt[g])
    hold = yc[yc.grp.isin(locked)]
    lock = dict(rule="md5('v6ysi:'+scaffold group) ascending; lock a group while locked YSI compounds <= 25 % of YSI compounds",
                n_ysi=len(yc), n_groups=len(cnt), locked_groups=locked, n_locked=len(hold), locked_inchikeys=sorted(hold.inchikey),
                pool_sha256=hashlib.sha256((R / "results/v6/gnn/pool.csv").read_bytes()).hexdigest())
    json.dump(lock, open(OUT / "lock.json", "w"), indent=1)
    print(len(yc), "YSI compounds,", len(cnt), "groups; locked", len(hold), "compounds in", len(locked), "groups")
    sys.exit()

# ---------------------------------------------------------------- score (once)
import torch
from fuelblend.ml import dmpnn as D, ysi_model as Y
flag = OUT / "scored.flag"
assert not flag.exists(), "holdout already scored once"
lock = json.load(open(OUT / "lock.json"))
assert hashlib.sha256((R / "results/v6/gnn/pool.csv").read_bytes()).hexdigest() == lock["pool_sha256"]
hold_ik = set(lock["locked_inchikeys"])
is_h = comp.inchikey.isin(hold_ik).to_numpy()
tr_all = np.where(~is_h)[0]; te = np.where(is_h)[0]
tasks = sorted(pool.target.unique())
Yraw = pd.pivot_table(pool, index="inchikey", columns="target", values="value").reindex(comp.inchikey)[tasks].values
Mask = ~np.isnan(Yraw)
feats = [D.featurize(s) for s in comp.smiles]
CFGS = [dict(hidden=200, depth=4, dropout=0.1), dict(hidden=300, depth=3, dropout=0.2)]
g = comp.grp.values
rng = np.random.default_rng(0); ug = np.unique(g[tr_all]); vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
va = np.array([i for i in tr_all if g[i] in vg]); tr = np.array([i for i in tr_all if g[i] not in vg])
mu = np.nan_to_num(np.nanmean(np.where(Mask[tr], Yraw[tr], np.nan), 0)); sd = np.nanstd(np.where(Mask[tr], Yraw[tr], np.nan), 0) + 1e-9
sd = np.where(np.isnan(sd), 1, sd)
Z = np.nan_to_num((Yraw - mu) / sd)
best = None
for c in CFGS:
    ms = D.fit_ensemble([feats[i] for i in tr], Z[tr], Mask[tr], [feats[i] for i in va], Z[va], Mask[va], n_models=5, epochs=150, seed=0, **c)
    pv, _ = D.predict(ms, [feats[i] for i in va])
    vl = float((((pv - Z[va]) ** 2) * Mask[va]).sum() / max(Mask[va].sum(), 1))
    print("cfg", c, "val", round(vl, 4), flush=True)
    if best is None or vl < best[0]: best = (vl, ms, c)
p, s = D.predict(best[1], [feats[i] for i in te])
j = tasks.index("ysi")
pred = p[:, j] * sd[j] + mu[j]
y = Yraw[te, j]
ytr = Yraw[tr_all, j]; trm = tr_all[Mask[tr_all, j]]
const = np.full(len(te), np.median(Yraw[trm, j]))
Xtr = np.vstack([Y.features(s_) for s_ in comp.smiles[trm]])
ridge = Y.fit_ridge(Xtr, Yraw[trm, j], g[trm]).predict(np.vstack([Y.features(s_) for s_ in comp.smiles[te]]))
grp = g[te]
brng = np.random.default_rng(0)
def boot(d):
    ugs = np.unique(grp); idx = {q: np.where(grp == q)[0] for q in ugs}
    b = [d[np.concatenate([idx[ugs[i]] for i in brng.integers(0, len(ugs), len(ugs))])].mean() for _ in range(2000)]
    return [float(np.percentile(b, 5)), float(np.percentile(b, 95))]
e = np.abs(pred - y); ec = np.abs(const - y); er = np.abs(ridge - y)
res = dict(n=int(len(te)), groups=int(len(np.unique(grp))), config=best[2], mae_gnn=float(e.mean()), mae_const=float(ec.mean()), mae_ridge=float(er.mean()),
           gain_vs_const_ci90=boot(ec - e), gain_vs_ridge_ci90=boot(er - e), over_pred_mean=float(np.maximum(pred - y, 0).mean()))
res["gate_pass"] = bool(res["gain_vs_const_ci90"][0] > 0 and res["gain_vs_ridge_ci90"][0] > 0)
# conformal-style absolute 90 % margin from the existing Track B out-of-fold YSI errors on non-holdout compounds (disclosed: those folds saw holdout compounds in training)
oo = pd.read_csv(R / "results/v6/gnn/oof.csv")
oo = oo[oo.ysi_obs.notna() & ~oo.inchikey.isin(hold_ik)]
q = float(np.quantile(np.abs(oo.ysi_oof - oo.ysi_obs), 0.9))
res.update(margin90_abs_dev_oof=q, holdout_coverage=float((e <= q).mean()), n_margin=int(len(oo)))
pd.DataFrame(dict(inchikey=comp.inchikey[te].values, smiles=comp.smiles[te].values, grp=grp, ysi=y, gnn=pred, const=const, ridge=ridge)).to_csv(OUT / "holdout_pred.csv", index=False)
json.dump(res, open(OUT / "readings.json", "w"), indent=1)
flag.write_text("scored once")
print(json.dumps(res, indent=1))
