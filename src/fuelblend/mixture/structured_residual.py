"""v4 Task 11: physics-structured viscosity residual, plus fuel-grouped nested cross-validation helpers.

Model (ln = natural log of kinematic viscosity in cSt; physics baseline = raw Refutas):

    delta ln(nu)(z, T) = a0(z) + a1(z) * 1000/T

z is a vector of fuel descriptors (class mass fractions, mean/std carbon number, ...).  a0 and a1 are each a
ridge (affine) or GP function of z.  The T-dependence is therefore exactly linear in 1000/T *by construction*:
the correction cannot develop curvature or blow up when extrapolated to -40 C, unlike an unconstrained net.

  kind="ridge": design [1, t, zs, zs*t] (t = 1000/T centred, zs = standardised z); the intercept and the pure-t
                column (the "recalibration" terms) are unpenalised, the z-dependent columns carry an L2 penalty `alpha`.
  kind="gp":    stage 1 = OLS on [1, t] (recalibration), stage 2 = zero-mean GP on the remaining residual with
                kernel  s2 * RBF(z, z'; length_scale) * (1 + ts*ts')  i.e. a0 and a1 are independent GPs sharing the
                same z-kernel.  The posterior mean is exactly linear in t.  Predictive std (observation) is returned
                on request; stage-1 coefficient uncertainty is ignored (documented limitation).

Nested-CV helpers: `select_by_inner_cv` / `nested_loo` take a `fit_predict(train_idx, test_idx, hp)` callable and
group rows by fuel.  Hyper-parameter selection for an outer held-out fuel only ever receives the training-fold
row indices (ledger ruling R6), and the held-out fuel's targets are never read before the final scoring.
"""
from __future__ import annotations

import numpy as np

__all__ = ["StructuredResidual", "group_kfold", "select_by_inner_cv", "nested_loo"]


class StructuredResidual:
    def __init__(self, features: list[str], kind: str = "ridge", alpha: float = 1.0,
                 length_scale: float = 4.0, noise_ratio: float = 0.3):
        if kind not in ("ridge", "gp"):
            raise ValueError("kind must be 'ridge' or 'gp'")
        self.features, self.kind = list(features), kind
        self.alpha, self.length_scale, self.noise_ratio = float(alpha), float(length_scale), float(noise_ratio)

    # ------------------------------------------------------------------ helpers
    def _mat(self, X) -> np.ndarray:
        if hasattr(X, "loc") and hasattr(X, "columns"):
            X = X[self.features].to_numpy(float)
        X = np.asarray(X, float)
        if X.ndim != 2 or X.shape[1] != len(self.features):
            raise ValueError(f"X must have {len(self.features)} columns ({self.features})")
        return X

    def _zs(self, X):
        return (self._mat(X) - self.mu_) / self.sd_

    def _t(self, T_K):
        return 1000.0 / np.asarray(T_K, float) - self.tm_

    # ------------------------------------------------------------------ fit / predict
    def fit(self, X, T_K, y_resid):
        X = self._mat(X)
        y = np.asarray(y_resid, float)
        T_K = np.asarray(T_K, float)
        if not (len(X) == len(T_K) == len(y)) or not np.all(np.isfinite(y)):
            raise ValueError("X, T_K, y_resid must be aligned and finite")
        self.mu_ = X.mean(0)
        sd = X.std(0)
        self.sd_ = np.where(sd > 1e-12, sd, 1.0)
        self.tm_ = float((1000.0 / T_K).mean())
        zs, t = self._zs(X), self._t(T_K)
        if self.kind == "ridge":
            D = np.column_stack([np.ones(len(y)), t, zs, zs * t[:, None]])
            P = np.full(D.shape[1], self.alpha)
            P[:2] = 0.0
            self.beta_ = np.linalg.solve(D.T @ D + np.diag(P), D.T @ y)
        else:
            B = np.column_stack([np.ones(len(y)), t])
            self.c_ = np.linalg.lstsq(B, y, rcond=None)[0]
            r = y - B @ self.c_
            self.s2_ = max(float(r.var()), 1e-10)
            self.ts_ = max(float(t.std()), 1e-12)
            self.Z_, self.tt_ = zs, t / self.ts_
            K = self.s2_ * self._rbf(zs, zs) * (1.0 + np.outer(self.tt_, self.tt_))
            self.noise_ = self.noise_ratio * self.s2_
            self.L_ = np.linalg.cholesky(K + self.noise_ * np.eye(len(y)))
            self.w_ = np.linalg.solve(self.L_.T, np.linalg.solve(self.L_, r))
        return self

    def _rbf(self, A, B):
        d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
        return np.exp(-0.5 * d2 / self.length_scale ** 2)

    def predict(self, X, T_K, return_std: bool = False):
        X = self._mat(X)
        T_K = np.broadcast_to(np.asarray(T_K, float), (len(X),))
        zs, t = self._zs(X), self._t(T_K)
        if self.kind == "ridge":
            D = np.column_stack([np.ones(len(t)), t, zs, zs * t[:, None]])
            mean = D @ self.beta_
            return (mean, None) if return_std else mean
        ts = t / self.ts_
        Ks = self._rbf(zs, self.Z_)                                   # (m, n)
        kstar = self.s2_ * Ks * (1.0 + np.outer(ts, self.tt_))
        mean = self.c_[0] + self.c_[1] * t + kstar @ self.w_
        if not return_std:
            return mean
        v = np.linalg.solve(self.L_, kstar.T)
        var = self.s2_ * (1.0 + ts ** 2) - (v ** 2).sum(0) + self.noise_
        return mean, np.sqrt(np.maximum(var, 1e-12))

    def a0_a1(self, X):
        """Return (a0(z), a1(z)) in the 1000/T variable (not centred): delta = a0 + a1 * 1000/T exactly."""
        X = self._mat(X)
        # evaluate at invT = tm+1 and tm: slope = difference, intercept = value at invT=0
        f1 = self.predict(X, 1000.0 / (self.tm_ + 1.0))
        f0 = self.predict(X, 1000.0 / self.tm_)
        a1 = f1 - f0
        a0 = f0 - a1 * self.tm_
        return a0, a1


# ----------------------------------------------------------------------------- nested CV (grouped by fuel)
def group_kfold(fuel_ids, n_splits: int = 5, seed: int = 0):
    """Deterministic grouped K-fold over unique fuel ids: list of (train_pos, val_pos) index arrays into `fuel_ids`."""
    fuel_ids = np.asarray(fuel_ids)
    u = np.unique(fuel_ids)
    k = min(n_splits, len(u))
    if k < 2:
        raise ValueError("need at least 2 fuels for an inner split")
    perm = np.random.default_rng(seed).permutation(len(u))
    folds = np.array_split(u[perm], k)
    out = []
    for f in folds:
        va = np.isin(fuel_ids, f)
        out.append((np.where(~va)[0], np.where(va)[0]))
    return out


def select_by_inner_cv(train_idx, fuel_ids, grid, fit_predict, y, n_splits: int = 5, seed: int = 0):
    """Choose the hyper-parameter in `grid` with the lowest grouped-CV RMSE using ONLY the rows in `train_idx`.

    fit_predict(tr, te, hp) must use only rows `tr` to predict rows `te`; every index it receives is a subset of
    `train_idx`.  Returns (best_hp, scores list aligned with grid)."""
    train_idx = np.asarray(train_idx)
    sub_f = np.asarray(fuel_ids)[train_idx]
    scores = []
    for hp in grid:
        se, n = 0.0, 0
        for tr_p, va_p in group_kfold(sub_f, n_splits, seed):
            tr, va = train_idx[tr_p], train_idx[va_p]
            p = fit_predict(tr, va, hp)
            e = np.asarray(p, float) - np.asarray(y, float)[va]
            se += float(np.sum(e ** 2)); n += len(e)
        scores.append(np.sqrt(se / n) if n else np.inf)
    scores = np.array(scores)
    scores[~np.isfinite(scores)] = np.inf
    return grid[int(np.argmin(scores))], scores.tolist()


def nested_loo(fuel_ids, y, grid, fit_predict, outer_fuels=None, n_splits: int = 5, seed: int = 0,
               universe=None):
    """Leave-one-fuel-out with inner-CV hyper-parameter selection inside each training fold.

    `universe` (optional index array) restricts which rows may ever be trained on / selected over (e.g. dev fuels
    only, keeping a locked set out of every fit).  Returns (pred aligned with fuel_ids -- NaN for fuels not in
    outer_fuels, {fuel: chosen_hp})."""
    fuel_ids = np.asarray(fuel_ids)
    uni = np.arange(len(fuel_ids)) if universe is None else np.asarray(universe)
    pred = np.full(len(fuel_ids), np.nan)
    chosen = {}
    for f in (np.unique(fuel_ids[uni]) if outer_fuels is None else list(outer_fuels)):
        te = uni[fuel_ids[uni] == f]
        tr = uni[fuel_ids[uni] != f]
        hp, _ = select_by_inner_cv(tr, fuel_ids, grid, fit_predict, y, n_splits, seed)
        pred[te] = fit_predict(tr, te, hp)
        chosen[f] = hp
    return pred, chosen
