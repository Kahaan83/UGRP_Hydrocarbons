"""Blend properties from pure-compound predictions (linear mole-fraction rule). See results/v7/final/blend_from_pure/PREDECLARATION.md."""
import sys, json
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from fuelblend.api import predict_pure as PP
from fuelblend.v7 import guard

R = Path(__file__).resolve().parents[3]
OUT = R / "results/v7/final/blend_from_pure"
LOCKS = ("mixture_lock", "mixture_holdout51", "splits_locked")
locked = set().union(*[{guard.canon(i) for i in guard.locked_ids(n)} for n in LOCKS])

oof = pd.read_csv(R / "results/v7/track_b/curves_v1_long/oof.csv").set_index("inchikey")
comp = pd.read_csv(R / "data/unified/compounds.csv").set_index("inchikey")
rows = pd.read_csv(R / "results/v6/mixture/rows.csv"); rows = rows[rows.role == "dev"]
um = pd.read_csv(R / "data/unified/measurements_dedup.csv", low_memory=False)
um = um[(um.property == "dcn") & um.kind.isin(["binary", "multi"])]

T = []
for head, prop, task in (("leenhouts", "mon", "mon_leenhouts"), ("kaust", "mon", "mon_kaust"), ("kaust", "ron", "ron_kaust")):
    d = rows[(rows["head"] == head) & (rows.prop == prop)]
    T.append((task, prop, [(r.system_id, json.loads(r.comps), np.array(json.loads(r.x), float), float(r.y)) for r in d.itertuples()]))
T.append(("dcn_leenhouts", "dcn", [(r.system_id, json.loads(r.comp_inchikeys), np.array(json.loads(r.x_mole), float), float(r.value)) for r in um.itertuples()]))

cache = {}
def pure_pred(ik, prop):
    k = (ik, prop)
    if k not in cache:
        if ik in oof.index and not np.isnan(oof.loc[ik, f"{prop}_long"]):
            cache[k] = (float(oof.loc[ik, f"{prop}_long"]), "oof")
        elif ik in comp.index:
            v = PP.predict_pure(comp.loc[ik, "smiles"], props=[prop])[comp.loc[ik, "smiles"]][prop]
            cache[k] = (v["value"], "final") if v["status"] == "ok" else (np.nan, "fail")
        else:
            cache[k] = (np.nan, "fail")
    return cache[k]

def pure_obs(ik, prop):
    return float(oof.loc[ik, f"{prop}_obs"]) if ik in oof.index else np.nan

out, summ, ps = [], [], []
for task, prop, items in T:
    keep = [it for it in items if guard.canon(it[0]) not in locked]
    n_locked = len(items) - len(keep)
    sysid = np.array([k[0] for k in keep]); y = np.array([k[3] for k in keep])
    pred = np.full(len(keep), np.nan); ref = np.full(len(keep), np.nan); src = []
    for i, (_, iks, x, _) in enumerate(keep):
        pp = [pure_pred(ik, prop) for ik in iks]
        if any(np.isnan(p[0]) for p in pp) or abs(x.sum() - 1) > 1e-3:
            src.append("dropped"); continue
        xs = x / x.sum()
        pred[i] = float(np.dot(xs, [p[0] for p in pp])); src.append("oof_only" if all(p[1] == "oof" for p in pp) else "uses_final")
        ob = np.array([pure_obs(ik, prop) for ik in iks])
        if not np.isnan(ob).any(): ref[i] = float(np.dot(xs, ob))
    ok = ~np.isnan(pred)
    c0 = np.full(len(keep), np.nan)
    for s in np.unique(sysid[ok]):
        c0[sysid == s] = np.median(y[ok & (sysid != s)])
    ok &= ~np.isnan(c0)
    e_m, e_c = np.abs(pred - y)[ok], np.abs(c0 - y)[ok]
    g = sysid[ok]; ug = np.unique(g)
    d = e_c - e_m
    S = np.array([d[g == u].sum() for u in ug]); N = np.array([(g == u).sum() for u in ug], float)
    idx = np.random.default_rng(0).integers(0, len(ug), (4000, len(ug))); est = S[idx].sum(1) / N[idx].sum(1)
    p = 2 * min((est <= 0).mean(), (est >= 0).mean()); ps.append(p)
    rk = ok & ~np.isnan(ref)
    summ.append(dict(task=task, n_rows=int(ok.sum()), n_systems=int(len(ug)), n_dropped=int(len(keep) - ok.sum()), n_locked_removed=n_locked,
                     n_rows_oof_only=int(sum(1 for s_, o in zip(src, ok) if o and s_ == "oof_only")),
                     mae_model=e_m.mean(), mae_const=e_c.mean(), gain=d.mean(), lo=np.percentile(est, 5), hi=np.percentile(est, 95), p_boot=p,
                     n_ref=int(rk.sum()), mae_model_on_ref=float(np.abs(pred - y)[rk].mean()) if rk.any() else np.nan,
                     mae_measured_pure_blend=float(np.abs(ref - y)[rk].mean()) if rk.any() else np.nan))
    out.append(pd.DataFrame(dict(task=task, system_id=sysid, y=y, pred=pred, const=c0, ref_measured_blend=ref, src=src)))
o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
for r_, i in enumerate(o): run = max(run, min(1, (len(ps) - r_) * ps[i])); adj[i] = run
S = pd.DataFrame(summ); S["p_holm"] = adj; S["usable_exploratory"] = (S.p_holm < 0.10) & (S.gain > 0) & (S.lo > 0)
S.to_csv(OUT / "summary.csv", index=False, float_format="%.5g"); pd.concat(out).to_csv(OUT / "rows.csv", index=False, float_format="%.5g")
print(S.round(3).to_string())
