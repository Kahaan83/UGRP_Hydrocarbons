"""Task 13 adapter: mixture_gbr = the `desc_gbm` comparator of v6 Track C (`scripts/v6/mixture/run_cv.py: gbm_residual`, sklearn GradientBoostingRegressor
n_estimators 45, depth 3, lr 0.08) on mole-weighted RDKit descriptors + physics baseline (+ T for Boehm), residual target, GroupKFold(min(5, groups)) on dev rows.
Gate = per (head, prop) match of the refit out-of-fold prediction to `results/v6/mixture/dev_oof.csv` column `pred_desc_gbm` (relative 1e-6 on the MAE of each unit)."""
import sys, json
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/external_ingest")
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[3]


def _descriptors():
    """`descriptors` from the v6 script, imported by path (its argparse runs at import, so argv is set for the call)."""
    import importlib.util
    old = sys.argv
    sys.argv = ["run_cv", "--stage", "summary"]
    try:
        spec = importlib.util.spec_from_file_location("v6_mix_run_cv", ROOT / "scripts/v6/mixture/run_cv.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules["v6_mix_run_cv"] = mod
        spec.loader.exec_module(mod)
    finally:
        sys.argv = old
    return mod.descriptors


def run(nt, to, quick, out, tb, T):
    desc = _descriptors()
    rows = pd.read_csv(ROOT / "results/v6/mixture/rows.csv")
    ref = pd.read_csv(ROOT / "results/v6/mixture/dev_oof.csv")
    from _trees import nested_lgb, boot_gain
    cache = {}

    def D(smi):
        if smi not in cache:
            cache[smi] = desc(smi)
        return cache[smi]

    units, gate_units, summ, ok = [], [], {}, True
    for (head, prop), sub in rows.groupby(["head", "prop"]):
        dev = sub[sub.role == "dev"].reset_index(drop=True)
        r = ref[(ref["head"] == head) & (ref["prop"] == prop)].reset_index(drop=True)
        assert len(r) == len(dev) and np.allclose(r.y, dev.y), f"dev rows differ from v6 for {head}/{prop}"
        sm = [json.loads(s) for s in dev.smiles]
        xs = [np.array(json.loads(s), float) for s in dev.x]
        Dm = np.stack([(x[:, None] * np.stack([D(s) for s in ss])).sum(0) for ss, x in zip(sm, xs)])
        base, y, g = dev.base.to_numpy(float), dev.y.to_numpy(float), dev.group.to_numpy()
        useT = head == "boehm"
        F = np.hstack([Dm, base[:, None]] + ([dev.T_K.to_numpy(float)[:, None]] if useT else []))
        folds = np.full(len(dev), -1)
        for k, (_, te) in enumerate(GroupKFold(min(5, len(np.unique(g)))).split(np.arange(len(dev)), groups=g)):
            folds[te] = k
        p_orig = np.full(len(dev), np.nan)
        for k in np.unique(folds):
            tr, te = np.where(folds != k)[0], np.where(folds == k)[0]
            m = GradientBoostingRegressor(n_estimators=45, max_depth=3, learning_rate=0.08, random_state=42)
            m.fit(F[tr], y[tr] - base[tr])
            p_orig[te] = base[te] + m.predict(F[te])
        mae, base_mae = float(np.abs(p_orig - y).mean()), float(np.abs(r.pred_desc_gbm.to_numpy() - y).mean())
        rel = abs(mae - base_mae) / base_mae
        gate_units.append(dict(head=head, prop=prop, repro_mae=mae, base_mae=base_mae, rel_diff=rel, passed=bool(rel < 1e-6)))
        ok &= rel < 1e-6
        units.append((head, prop, F, y, base, g, folds, p_orig))
    gate = dict(model="mixture_gbr", units=gate_units, rel_diff=max(u["rel_diff"] for u in gate_units), passed=bool(ok))
    if not ok:
        return gate, None
    nsig, bps = 0, {}
    for head, prop, F, y, base, g, folds, p_orig in units:
        X = pd.DataFrame(F, columns=[f"f{i}" for i in range(F.shape[1])])
        oof, res = nested_lgb(X, y - base, np.ones(len(y)), np.ones(len(y), bool), g, folds, cat=None,
                              prefix=f"v7_mixture_gbr_{head}_{prop}" + ("_quick" if quick else ""), out_dir=out, nt=nt, to=to, tb=tb, quick=quick, n_jobs=1)
        p_new = base + oof
        e_new, e_orig = np.abs(p_new - y), np.abs(p_orig - y)
        gain, lo, hi = boot_gain(e_orig, e_new, g)
        nsig += int(lo > 0)
        k = f"{head}_{prop}"
        summ.update({f"{k}_n": len(y), f"{k}_n_groups": int(len(np.unique(g))), f"{k}_mae_v7": float(e_new.mean()), f"{k}_mae_refit_original": float(e_orig.mean()),
                     f"{k}_mae_physics_baseline": float(np.abs(base - y).mean()), f"{k}_gain": gain, f"{k}_gain_lo": lo, f"{k}_gain_hi": hi,
                     })
        bps[k] = res.best_params.tolist()
    summ = dict(model="mixture_gbr", metric="dev OOF MAE per (head, prop); gain = refit original minus v7, group bootstrap 90 % CI", units_with_gain_ci_above_0=nsig,
                n_units=len(units), best_params=bps, **summ)
    return gate, summ
