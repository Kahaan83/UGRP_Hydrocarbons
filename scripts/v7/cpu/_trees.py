"""Shared v7 helpers for Task 13 (CPU tree models): nested Optuna LightGBM with curve logging, and file-hash guards."""
import sys, hashlib, os
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np
import lightgbm as lgb
import optuna
from fuelblend.ml.curves import CurveLogger, lgb_curve_callback
from fuelblend.v7 import hpo

ROOT = Path(__file__).resolve().parents[3]


def tree_hash(*rels):
    """sha256 over all files under the given repo-relative dirs (sorted), used to assert nothing under src/ or BASE changed."""
    h = hashlib.sha256()
    for rel in rels:
        for p in sorted((ROOT / rel).rglob("*")):
            if p.is_file() and "__pycache__" not in p.parts:
                h.update(str(p.relative_to(ROOT)).encode())
                h.update(p.read_bytes())
    return h.hexdigest()


def lgb_params(p, n_jobs=4, seed=42):
    return dict(objective="regression", learning_rate=p["learning_rate"], num_leaves=int(p["num_leaves"]),
                min_data_in_leaf=int(p["min_child_samples"]), bagging_fraction=p["subsample"], bagging_freq=1,
                feature_fraction=p["colsample"], lambda_l2=p["lambda_l2"], seed=seed, deterministic=True, force_row_wise=True,
                num_threads=n_jobs, verbose=-1)


def _val_groups(grp, idx, seed):
    ug = np.unique(grp[idx])
    if len(ug) < 2:
        return np.zeros(len(idx), bool)
    vg = set(np.random.default_rng(seed).choice(ug, max(1, int(.15 * len(ug))), replace=False))
    return np.array([grp[i] in vg for i in idx])


def fit_es(X, y, w, grp, idx, params, cat, seed, logger=None, trial=None, early=100, n_est=5000):
    """lgb.train on rows idx with a 15 % group-held-out split for early stopping; returns (booster, best_iter, train/val loss at best)."""
    v = _val_groups(grp, idx, seed)
    if v.all() or not v.any():
        v = np.zeros(len(idx), bool)
        v[::7] = True
    tr, va = idx[~v], idx[v]
    dtr = lgb.Dataset(X.iloc[tr], y[tr], weight=w[tr], categorical_feature=cat or "auto", free_raw_data=False)
    dva = lgb.Dataset(X.iloc[va], y[va], weight=w[va], reference=dtr)
    cbs = [lgb.early_stopping(early, verbose=False)]
    if logger is not None:
        cbs.append(lgb_curve_callback(logger))
    if trial is not None:
        def prune(env):
            vl = float(env.evaluation_result_list[0][2])
            trial.report(vl, env.iteration)
            if trial.should_prune():
                raise optuna.TrialPruned()
        prune.order = 40
        cbs.append(prune)
    b = lgb.train({**params, "metric": "l2"}, dtr, num_boost_round=n_est, valid_sets=[dva], callbacks=cbs)
    return b


def nested_lgb(X, y, w, elig, grp, folds, *, cat, prefix, out_dir, nt, to, tb=None, quick=False, n_jobs=4, extra=None):
    """Nested CV over prediction units X (rows). Trainable rows are `elig`; outer folds `folds` (one id per row); inner GroupKFold(3) by hpo.
    Returns (oof prediction per row, per-fold search summary, curves path). Final fit per outer fold uses the best params with early stopping."""
    space = hpo.load_space("lgbm")
    if quick:
        space = {**space, "fixed": {"n_estimators": 60, "early_stopping_rounds": 10}}
    n_est, early = space["fixed"]["n_estimators"], space["fixed"]["early_stopping_rounds"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lg = CurveLogger(out_dir / f"curves_{prefix}.csv", tb_dir=tb)

    def fit_score(params, tr, va, logger, on_epoch=None):
        p = {**lgb_params(params, n_jobs), **(extra or {})}
        tre, vae = tr[elig[tr]], va[elig[va]]
        dtr = lgb.Dataset(X.iloc[tre], y[tre], weight=w[tre], categorical_feature=cat or "auto", free_raw_data=False)
        dva = lgb.Dataset(X.iloc[vae], y[vae], weight=w[vae], reference=dtr)
        cbs = [lgb.early_stopping(early, verbose=False)]
        if on_epoch is not None:
            def rep(env):
                on_epoch(env.iteration + 1, float(env.evaluation_result_list[0][2]))
            rep.order = 40
            cbs.append(rep)
        b = lgb.train({**p, "metric": "l2"}, dtr, num_boost_round=n_est, valid_sets=[dva], callbacks=cbs)
        return float(b.best_score["valid_0"]["l2"])

    def final_fit(params, tr, te):
        k = int(folds[te[0]])
        b = fit_es(X, y, w, grp, tr[elig[tr]], {**lgb_params(params, n_jobs), **(extra or {})}, cat, k, logger=lg.child(model="final", fold=k), early=early, n_est=n_est)
        return b.predict(X.iloc[te], num_iteration=b.best_iteration)

    res = hpo.nested_cv(None, grp, folds, fit_score=fit_score, space=space, n_trials=nt, timeout_s=to, study_prefix=prefix,
                        final_fit=final_fit, out_dir=out_dir / "preds" / prefix)
    lg.close()
    oof = np.full(len(y), np.nan)
    for k in np.unique(folds):
        oof[folds == k] = np.load(out_dir / "preds" / prefix / f"pred_fold{k}.npy")
    return oof, res


def boot_gain(e_old, e_new, groups, n=4000, seed=0):
    """Per-group MAE gain (old - new), mean and 90 % percentile interval over group resamples."""
    ug = np.unique(groups)
    a = np.array([np.mean(e_old[groups == u]) for u in ug])
    b = np.array([np.mean(e_new[groups == u]) for u in ug])
    d = a - b
    idx = np.random.default_rng(seed).integers(0, len(ug), (n, len(ug)))
    est = d[idx].mean(1)
    return float(d.mean()), float(np.percentile(est, 5)), float(np.percentile(est, 95))
