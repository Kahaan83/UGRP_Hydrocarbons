"""v4 Task 11: physics-structured viscosity residual + nested CV (R6: selection sees ONLY training-fold rows)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.mixture.structured_residual import (  # noqa: E402
    StructuredResidual, nested_loo, select_by_inner_cv)

FEATS = ["f1", "f2", "f3", "f4"]


def _synthetic(n_fuels=14, seed=3):
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n_fuels, len(FEATS)))
    T_C = np.array([20.0, -20.0, -40.0])
    X, T, y, fid = [], [], [], []
    for f in range(n_fuels):
        for tc in T_C:
            tk = tc + 273.15
            it = 1000.0 / tk
            a0 = 0.3 * Z[f, 0] - 0.2
            a1 = 0.05 * Z[f, 1] + 0.1
            X.append(Z[f]); T.append(tk); fid.append(f"F{f:02d}")
            y.append(a0 + a1 * it + rng.normal(scale=0.02))
    return np.array(X), np.array(T), np.array(y), np.array(fid)


@pytest.mark.parametrize("kind", ["ridge", "gp"])
def test_structured_residual_is_linear_in_invT(kind):
    X, T, y, _ = _synthetic()
    m = StructuredResidual(FEATS, kind=kind).fit(X, T, y)
    rng = np.random.default_rng(0)
    z = rng.normal(size=(5, len(FEATS)))
    inv = np.linspace(3.5, 5.0, 9)                    # equally spaced in 1000/T  (T from ~286 K to 200 K)
    for zi in z:
        pred = m.predict(np.tile(zi, (len(inv), 1)), 1000.0 / inv)
        d2 = np.diff(pred, n=2)
        assert np.max(np.abs(d2)) < 1e-9, d2
    # extrapolation far outside the training window stays finite and linear (cannot blow up)
    far = np.linspace(6.0, 9.0, 7)
    pred = m.predict(np.tile(z[0], (len(far), 1)), 1000.0 / far)
    assert np.all(np.isfinite(pred)) and np.max(np.abs(np.diff(pred, n=2))) < 1e-9


def test_ridge_recovers_known_structure():
    X, T, y, _ = _synthetic(n_fuels=30)
    m = StructuredResidual(FEATS, kind="ridge", alpha=1e-6).fit(X, T, y)
    assert np.sqrt(np.mean((m.predict(X, T) - y) ** 2)) < 0.05


def _make_fit_predict(X, T, y, fid):
    def fit_predict(tr, te, hp):
        m = StructuredResidual(FEATS, kind="ridge", alpha=hp).fit(X[tr], T[tr], y[tr])
        return m.predict(X[te], T[te])
    return fit_predict


GRID = [1e-3, 1e-1, 1.0, 10.0, 100.0, 1e3]


def test_nested_cv_selects_inside_fold():
    """R6: for each outer fold the selection function sees only training-fold rows; poisoning the
    held-out fuel's targets (NaN) changes neither the selected hyper-parameter nor the prediction."""
    X, T, y, fid = _synthetic()
    held = "F05"
    calls = []

    base_fp = _make_fit_predict(X, T, y, fid)

    def spying_fp(tr, te, hp):
        calls.append((np.asarray(tr).copy(), np.asarray(te).copy(), hp))
        return base_fp(tr, te, hp)

    pred1, chosen1 = nested_loo(fid, y, GRID, spying_fp, outer_fuels=[held])
    held_rows = set(np.where(fid == held)[0])
    # every call but the last is a selection call: held-out rows appear in neither its train nor its val set
    for tr, te, _ in calls[:-1]:
        assert not (held_rows & set(tr)) and not (held_rows & set(te))
    tr, te, _ = calls[-1]                              # final fit on the training fold, evaluated on the held-out fuel
    assert not (held_rows & set(tr)) and set(te) == held_rows
    assert chosen1[held] in GRID

    y_poison = y.copy()
    y_poison[fid == held] = np.nan
    pred2, chosen2 = nested_loo(fid, y_poison, GRID, _make_fit_predict(X, T, y_poison, fid), outer_fuels=[held])
    assert chosen2 == chosen1
    assert np.all(np.isfinite(pred2[fid == held]))
    assert np.allclose(pred1[fid == held], pred2[fid == held])


def test_nested_cv_full_loo_returns_all_folds():
    X, T, y, fid = _synthetic()
    pred, chosen = nested_loo(fid, y, GRID, _make_fit_predict(X, T, y, fid))
    assert set(chosen) == set(fid) and np.all(np.isfinite(pred))
    # selection is allowed to agree across folds (documented, R6); each choice must come from the grid
    assert all(v in GRID for v in chosen.values())


def test_select_by_inner_cv_uses_only_given_rows():
    X, T, y, fid = _synthetic()
    keep = np.where(fid != "F00")[0]
    seen = set()

    def fp(tr, te, hp):
        seen.update(tr.tolist()); seen.update(te.tolist())
        return _make_fit_predict(X, T, y, fid)(tr, te, hp)

    select_by_inner_cv(keep, fid, GRID, fp, y)
    assert seen <= set(keep.tolist())


def test_nested_cv_universe_excludes_locked_rows():
    """A locked set passed outside `universe` is never trained on or selected over."""
    X, T, y, fid = _synthetic()
    locked = np.where(np.isin(fid, ["F00", "F01"]))[0]
    uni = np.setdiff1d(np.arange(len(fid)), locked)
    seen = set()

    def fp(tr, te, hp):
        seen.update(tr.tolist()); seen.update(te.tolist())
        return _make_fit_predict(X, T, y, fid)(tr, te, hp)

    pred, chosen = nested_loo(fid, y, GRID, fp, outer_fuels=["F05", "F06"], universe=uni)
    assert not (seen & set(locked.tolist())) and set(chosen) == {"F05", "F06"}
