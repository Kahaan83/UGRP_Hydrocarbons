"""Small tests for scripts/v5/s5_freezing_nalkane.py (F0 is a max over species; censored values excluded; LOGO leak-free)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("s5", ROOT / "scripts" / "v5" / "s5_freezing_nalkane.py")
s5 = importlib.util.module_from_spec(_spec)
sys.modules["s5"] = s5
_spec.loader.exec_module(s5)

CONSTS = {14: dict(tm_K=279.0, dh_J_mol=45000.0), 16: dict(tm_K=291.0, dh_J_mol=53000.0)}


def _single_T(n, mol, total):
    return 1.0 / (1.0 / CONSTS[n]["tm_K"] - s5.R_GAS / CONSTS[n]["dh_J_mol"] * np.log(mol / total))


def test_f0_is_max_over_species():
    mol = {14: 0.05, 16: 0.01}
    T, arg = s5.t_first(mol, 1.0, CONSTS)
    t14, t16 = _single_T(14, 0.05, 1.0), _single_T(16, 0.01, 1.0)
    assert np.isclose(T, max(t14, t16))
    assert arg == (14 if t14 > t16 else 16)
    # pure species at x = 1 freezes at its Tm
    T1, a1 = s5.t_first({16: 1.0}, 1.0, CONSTS)
    assert np.isclose(T1, 291.0) and a1 == 16


def test_f0_inert_solvent_and_missing_constants():
    T, _ = s5.t_first({14: 0.05}, 1.0, CONSTS)
    T_more_solvent, _ = s5.t_first({14: 0.05}, 2.0, CONSTS)       # same moles, more inert solvent -> lower x -> lower T
    assert T_more_solvent < T
    T_none, arg = s5.t_first({9: 0.3}, 1.0, CONSTS)               # no constants for C9: non-solid-forming
    assert np.isnan(T_none) and arg is None


def test_censored_qualifiers_excluded():
    long = pd.DataFrame(dict(fuel_id=[1, 1, 2, 2, 3], property="freezing_point_C", value=[-50.0, -60.0, -70.0, -40.0, -55.0],
                             qualifier=[np.nan, "<", "<", "<", np.nan], method=["D5972"] * 5))
    f = s5.freezing_ajftd(long).set_index("fuel_id")
    assert f.loc["1", "value"] == -50.0 and f.loc["1", "n_censored"] == 1 and not f.loc["1", "censored_only"]
    assert np.isnan(f.loc["2", "value"]) and f.loc["2", "censored_only"]
    assert f.loc["3", "value"] == -55.0
    bulk = pd.DataFrame(dict(fuel_id=["a", "b"], coa_freezing_point_C=[-60.0, -50.0], coa_freezing_point_C_qualifier=["<", np.nan]))
    v = s5.freezing_av3322(bulk).set_index("fuel_id")
    assert np.isnan(v.loc["a", "value"]) and v.loc["a", "censored_only"] and v.loc["b", "value"] == -50.0


def test_logo_never_uses_scored_group_value():
    rng = np.random.default_rng(1)
    n = 12
    y = rng.normal(-50, 5, n)
    groups = np.array([f"g{i // 2}" for i in range(n)])          # pairs of rows share a group
    tf = y * 0.8 + rng.normal(0, 1, n) - 5
    X = rng.normal(size=(n, 5))
    base = s5.logo_predict(y, groups, tf, X)
    for gname in np.unique(groups):
        m = groups == gname
        y2 = y.copy()
        y2[m] = 1e6                                               # corrupt only the scored group's target
        alt = s5.logo_predict(y2, groups, tf, X)
        for col in ("C", "F1", "P"):
            assert np.allclose(base.loc[m, col], alt.loc[m, col]), col
    # and it does depend on the other groups
    y3 = y.copy()
    y3[groups == "g0"] = 1e6
    alt = s5.logo_predict(y3, groups, tf, X)
    assert not np.allclose(base.loc[groups == "g3", "C"], alt.loc[groups == "g3", "C"])


def test_locked_rows_not_trained_on():
    y = np.array([-50.0, -52.0, -48.0, -10.0])
    groups = np.array(["a", "b", "c", "L"])
    train = np.array([True, True, True, False])
    r = s5.logo_predict(y, groups, np.array([-40.0, -42.0, -38.0, -10.0]), np.zeros((4, 5)), train_mask=train, pred_idx=np.array([3]))
    assert np.isclose(r.C[3], np.mean(y[:3]))
