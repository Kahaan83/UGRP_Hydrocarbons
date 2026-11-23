"""v11: tune LightGBM settings of the viscosity pair residual (dev systems only). Protocol: docs/v11_predeclaration.md."""
import sys, json, itertools, warnings, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
OUT = ROOT / "results/v11/tune_viscosity"
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("dmv", ROOT / "scripts/v5/base/delta_ml_viscosity.py")
D = importlib.util.module_from_spec(spec)
sys.modules["dmv"] = D
spec.loader.exec_module(D)

from fuelblend.ml import viscosity_pair as vp  # noqa: E402
import lightgbm as lgb  # noqa: E402
from sklearn.model_selection import GroupKFold  # noqa: E402

GRID = dict(n_estimators=[200, 400, 800], learning_rate=[0.02, 0.05, 0.1], num_leaves=[7, 15, 31],
            min_child_samples=[20, 40, 80], colsample_bytree=[0.6, 0.8, 1.0], reg_lambda=[1.0, 5.0, 20.0])


def oof_err(R, P, params, folds):
    """Per-system macro |error| of corrected ln viscosity (ML) on dev warm rows; folds = list of (train_rows, test_rows) index arrays."""
    n = len(R)
    sysarr = R.system_id.to_numpy()
    resid_ref = R.ln_pred.to_numpy() - R.ln_meas.to_numpy()
    oof = np.zeros(n)
    for tr, te in folds:
        Ptr = D.train_set(R, P, set(sysarr[tr]))
        g, w = D.targets(R, Ptr)
        m = lgb.LGBMRegressor(**params)
        m.fit(Ptr[vp.FEATURES], g, sample_weight=w, categorical_feature=["class_pair"])
        Pte = P[P.row.isin(te)]
        oof[te] = D.corr_by_row(Pte, m.booster_.predict(Pte[vp.FEATURES]), n)[te]
    return resid_ref + oof


def main():
    ctx = vp.context()
    R = D.regenerate_rows(ctx)
    P = D.build_pairs(R, ctx)
    sysarr = R.system_id.to_numpy()
    dev_rows = np.where(R.role.to_numpy() == "dev")[0]
    warm = R.T_K.to_numpy() >= D.T_WARM
    mask = np.zeros(len(R), bool)
    mask[dev_rows] = True
    mask &= warm
    base = dict(vp.LGB_PARAMS)

    def v5_folds():
        return [(dev_rows[a], dev_rows[b]) for a, b in GroupKFold(n_splits=5).split(dev_rows, groups=sysarr[dev_rows])]

    def shuffled_folds(seed):
        us = np.array(sorted(set(sysarr[dev_rows])))
        rng = np.random.default_rng(seed)
        fold_of = dict(zip(us, rng.permutation(len(us)) % 5))
        f = np.array([fold_of[s] for s in sysarr[dev_rows]])
        return [(dev_rows[f != k], dev_rows[f == k]) for k in range(5)]

    rng = np.random.default_rng(11)
    keys = list(GRID)
    cfgs = [dict(zip(keys, [GRID[k][rng.integers(len(GRID[k]))] for k in keys])) for _ in range(30)]
    rows = []
    folds = v5_folds()
    for i, c in enumerate([{}] + cfgs):
        p = {**base, **c, "n_jobs": 4}
        e = D.macro(R, oof_err(R, P, p, folds), mask)
        rows.append(dict(cfg=i, **c, macro_mae=float(e.mean())))
        print(i, c, round(float(e.mean()), 5), flush=True)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "search.csv", index=False)
    best = T.iloc[1:].sort_values("macro_mae").iloc[0]
    bp = {k: type(GRID[k][0])(best[k]) for k in keys}
    print("best", bp, flush=True)
    diffs = []
    for s in (1, 2, 3):
        f = shuffled_folds(s)
        ed = D.macro(R, oof_err(R, P, {**base, "n_jobs": 4}, f), mask)
        et = D.macro(R, oof_err(R, P, {**base, **bp, "n_jobs": 4}, f), mask)
        diffs.append((ed - et.reindex(ed.index)))
        print("confirm seed", s, float(ed.mean()), float(et.mean()), flush=True)
    d = pd.concat(diffs).groupby(level=0).mean()
    brng = np.random.default_rng(42)
    bs = np.array([d.to_numpy()[brng.integers(0, len(d), len(d))].mean() for _ in range(2000)])
    res = dict(best=bp, n_systems=int(len(d)), mean_gain=float(d.mean()), ci_lo=float(np.percentile(bs, 2.5)),
               ci_hi=float(np.percentile(bs, 97.5)), v5_search_mae=float(T.macro_mae.iloc[0]), best_search_mae=float(best.macro_mae))
    res["adopt"] = bool(res["ci_lo"] > 0 and res["mean_gain"] >= 0.002)
    (OUT / "summary.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
