"""Shared v7 helpers for the MolPool heads (Task 11): nested Optuna on the residual over a physics baseline."""
import sys, os
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np
from fuelblend.ml import molpool as MP
from fuelblend.ml.curves import CurveLogger
from fuelblend.v7 import hpo

DEVICE = "cuda" if os.environ.get("FUELBLEND_DEVICE") == "cuda" else None


def _val_pos(groups, idx, seed):
    """Positions (within idx) of a 15 % group-held-out early-stopping split; None if it cannot be formed."""
    ug = np.unique(groups[idx])
    if len(ug) < 2:
        return None
    vg = set(np.random.default_rng(seed).choice(ug, max(1, int(.15 * len(ug))), replace=False))
    v = np.array([j for j, i in enumerate(idx) if groups[i] in vg])
    return v if 0 < len(v) < len(idx) else None


def fit_resid(H, X, A, r, groups, tr, te, params, seeds, seed, logger=None, on_epoch=None):
    """Fit on rows tr (standardised residual), predict rows te; returns (prediction in residual units, scale)."""
    s = float(r[tr].std())
    s = s if s > 1e-9 else 1.0
    A_ = None
    if A is not None:
        A_ = (A - A[tr].mean(0)) / (A[tr].std(0) + 1e-9)
    p = dict(params)
    ep, pat = p.pop("epochs"), p.pop("patience")
    pr = MP.fit_predict(H[tr], X[tr], r[tr] / s, H[te], X[te], None if A_ is None else A_[tr], None if A_ is None else A_[te],
                        seeds=seeds, epochs=ep, patience=pat, val_idx=_val_pos(groups, tr, seed), device=DEVICE,
                        logger=logger, on_epoch=on_epoch, **p)
    return pr * s


def tune(H, X, A, r, groups, folds, *, prefix, out_dir, nt, to, tb, quick=False, n_final=5, on_fold=None):
    """Nested CV: returns OOF residual predictions (one per row) and the per-fold search summary."""
    space = hpo.load_space("molpool")
    if quick:
        space = {**space, "fixed": {"epochs": 6, "patience": 3}}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lg = CurveLogger(out_dir / f"curves_{prefix}.csv", tb_dir=tb)

    def fit_score(params, tr, va, logger, on_epoch=None):
        pr = fit_resid(H, X, A, r, groups, tr, va, params, (0,), 0, on_epoch=on_epoch)
        s = float(r[tr].std())
        s = s if s > 1e-9 else 1.0
        return float(np.mean((pr / s - r[va] / s) ** 2))

    def final_fit(params, tr, te):
        k = int(folds[te[0]])
        return fit_resid(H, X, A, r, groups, tr, te, params, tuple(range(n_final)), k, logger=lg.child(model="final", fold=k))

    res = hpo.nested_cv(None, groups, folds, fit_score=fit_score, space=space, n_trials=nt, timeout_s=to, study_prefix=prefix,
                        final_fit=final_fit, out_dir=out_dir / "preds" / prefix, on_fold=on_fold)
    lg.close()
    oof = np.full(len(r), np.nan)
    for k in np.unique(folds):
        oof[folds == k] = np.load(out_dir / "preds" / prefix / f"pred_fold{k}.npy")
    return oof, res


def boot_p(d, g, n=4000, seed=0):
    """Mean of d with a 90 % group-bootstrap CI and a two-sided bootstrap p (d = comparator error minus candidate error)."""
    ug = np.unique(g)
    S = np.array([d[g == u].sum() for u in ug])
    N = np.array([(g == u).sum() for u in ug], float)
    idx = np.random.default_rng(seed).integers(0, len(ug), (n, len(ug)))
    est = S[idx].sum(1) / N[idx].sum(1)
    return float(d.mean()), float(np.percentile(est, 5)), float(np.percentile(est, 95)), float(2 * min((est <= 0).mean(), (est >= 0).mean()))


def holm(ps):
    o = np.argsort(ps)
    adj, run = np.empty(len(ps)), 0.0
    for r_, i in enumerate(o):
        run = max(run, min(1.0, (len(ps) - r_) * ps[i]))
        adj[i] = run
    return adj
