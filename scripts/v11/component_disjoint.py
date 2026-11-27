"""Exploratory component-disjoint (cold-start molecule) test of the viscosity (single-regime default, T>=283.15 K) and density (tuned)
pair corrections. Post hoc: pools development and locked systems; training excludes every system containing any held-out compound.
Hold-out units: (a) random compound folds (5 folds, seeds 0-2); (b) whole chemical classes. Does not write to results/v5/BASE."""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/v11")); sys.path.insert(0, str(ROOT / "src")); warnings.simplefilter("ignore")
import nested_tune as N
from fuelblend.ml import viscosity_pair as vp

OUT = ROOT / "results/v11/component_disjoint"
ctx = vp.context(); reg = ctx.reg; CLS = np.asarray(ctx.cls)
MARGIN = {"visc": 0.1343, "dens": 4.5424}
dp = json.load(open(ROOT / "results/v11/nested/dens_final_params.json"))
PARAMS = {"visc": ({**vp.LGB_PARAMS, "n_jobs": 4}, 2), "dens": ({**dp["lgb"]}, dp["weight_exponent"])}

def comps(S, which):
    if which == "visc":
        return [set(map(int, t)) for t in S.R.ii]
    return [set(reg.ik2i[k] for k in t) for t in S.R.iks]

def run(which):
    S = N.Setup(which); R = S.R; sysarr = S.sysarr
    rowc = comps(S, which)
    sc = {}
    for s, c in zip(sysarr, rowc): sc[s] = c
    systems = sorted(sc); allc = sorted(set().union(*sc.values()))
    mask = (R.T_K.to_numpy() >= S.D.T_WARM) if which == "visc" else np.ones(len(R), bool)
    mask &= np.array([len(c) >= 1 for c in rowc])
    prm, p = PARAMS[which]
    def evaluate(held):
        test_sys = {s for s in systems if sc[s] & held}
        train_sys = set(systems) - test_sys
        rows = np.where(np.isin(sysarr, list(test_sys)))[0]
        corr = S.fit_predict(prm, p, train_sys, rows)
        err = S.err_from(corr)
        mk = np.zeros(len(R), bool); mk[rows] = True; mk &= mask
        base = S.D.macro(R, S.base_err, mk); ml = S.D.macro(R, err, mk)
        worst = pd.Series(np.abs(err[mk])).groupby(sysarr[mk]).max()
        nun = pd.Series({s: len(sc[s] & held) for s in ml.index}); ncomp = pd.Series({s: len(sc[s]) for s in ml.index})
        return pd.DataFrame(dict(base=base, ml=ml, covered=(worst.reindex(ml.index) <= MARGIN[which]).astype(int), n_unseen=nun, n_comp=ncomp)), len(train_sys)
    res = []
    for seed in (0, 1, 2):
        rng = np.random.default_rng(seed); fold = dict(zip(allc, rng.integers(0, 5, len(allc))))
        parts = []
        for f in range(5):
            d, ntr = evaluate({c for c in allc if fold[c] == f}); d["seed"] = seed; d["fold"] = f; d["n_train_sys"] = ntr; parts.append(d)
            print(which, "seed", seed, "fold", f, len(d), flush=True)
        res.append(pd.concat(parts))
    cls_res = {}
    for k in sorted(set(CLS)):
        held = {c for c in allc if CLS[c] == k}
        if not held: continue
        d, ntr = evaluate(held); d["n_train_sys"] = ntr; cls_res[str(k)] = d; print(which, "class", k, len(d), ntr, flush=True)
    return pd.concat(res), cls_res

def summarise(d, rng):
    g = d.groupby(level=0).mean(numeric_only=True)  # one row per system (average over appearances)
    gain = (g.base - g.ml).to_numpy()
    bs = np.array([gain[rng.integers(0, len(gain), len(gain))].mean() for _ in range(5000)])
    from scipy.stats import beta
    k, n = int(g.covered.round().sum()), len(g)
    return dict(n_systems=n, base_mae=float(g.base.mean()), ml_mae=float(g.ml.mean()), gain=float(gain.mean()),
                ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                cov=f"{k}/{n}", cov_frac=k / n, cp=[float(beta.ppf(.025, k, n - k + 1)) if k else 0.0, float(beta.ppf(.975, k + 1, n - k)) if k < n else 1.0])

out = {}
for which in ("visc", "dens"):
    res, cls_res = run(which)
    rng = np.random.default_rng(1)
    o = {"seeds": {int(s): summarise(res[res.seed == s], rng) for s in (0, 1, 2)}}
    r0 = res[res.seed == 0]
    o["seed0_all_comps_unseen"] = summarise(r0[r0.n_unseen == r0.n_comp], rng) if (r0.n_unseen == r0.n_comp).any() else None
    o["seed0_one_comp_unseen"] = summarise(r0[r0.n_unseen < r0.n_comp], rng) if (r0.n_unseen < r0.n_comp).any() else None
    o["class_holdout"] = {k: summarise(v, rng) | dict(n_train_sys=int(v.n_train_sys.iloc[0])) for k, v in cls_res.items()}
    out[which] = o
    res.to_csv(OUT / f"{which}_folds.csv")
(OUT / "summary.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
