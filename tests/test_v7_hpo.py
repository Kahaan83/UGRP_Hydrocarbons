import json

import numpy as np
import optuna

from fuelblend.v7 import hpo

SPACE = {"fixed": {"epochs": 3}, "params": {"a": {"type": "float", "low": -1, "high": 1}}}


def _data():
    g = np.repeat(np.arange(20), 3)
    folds = np.array([i % 4 for i in g])
    return g, folds


def _st(tmp_path):
    return f"sqlite:///{(tmp_path / 'o.db').as_posix()}"


def test_study_resume_no_duplicate_trials(tmp_path):
    g, folds = _data()

    def fs(p, tr, va, lg):
        return p["a"] ** 2
    kw = dict(fit_score=fs, space=SPACE, timeout_s=60, study_prefix="t", storage=_st(tmp_path))
    hpo.nested_cv(None, g, folds, n_trials=4, **kw)
    hpo.nested_cv(None, g, folds, n_trials=6, **kw)
    for k in range(4):
        assert len(optuna.load_study(study_name=f"t_fold{k}", storage=_st(tmp_path)).trials) == 6


def test_outer_never_in_inner(tmp_path):
    g, folds = _data()
    seen = []

    def fs(p, tr, va, lg):
        used = set(tr) | set(va)
        assert not (set(tr) & set(va))
        assert any(not (used & set(np.where(folds == k)[0])) for k in range(4))
        assert not (set(g[tr]) & set(g[va]))
        seen.append(1)
        return p["a"] ** 2
    hpo.nested_cv(None, g, folds, fit_score=fs, space=SPACE, n_trials=2, timeout_s=60, study_prefix="n",
                  storage=_st(tmp_path))
    assert len(seen) == 4 * 2 * 3


def test_outer_scored_once(tmp_path):
    g, folds = _data()
    calls = []

    def ff(bp, tr, te):
        calls.append(1)
        return np.zeros(len(te))

    def fs(p, tr, va, lg):
        return p["a"] ** 2
    kw = dict(fit_score=fs, space=SPACE, n_trials=2, timeout_s=60, study_prefix="s",
              storage=_st(tmp_path), final_fit=ff, out_dir=tmp_path / "out")
    r1 = hpo.nested_cv(None, g, folds, **kw)
    r2 = hpo.nested_cv(None, g, folds, **kw)
    assert len(calls) == 4 and list(r1.fold) == list(r2.fold)
    assert json.loads(r1.best_params[0])["epochs"] == 3


def test_prune_hook_reports_from_inner_fold_zero(tmp_path):
    g, folds = _data()
    got = []

    def fs(p, tr, va, lg, on_epoch=None):
        for ep in range(1, 4):
            if on_epoch:
                on_epoch(ep, 1.0)
                got.append(ep)
        return 1.0
    hpo.nested_cv(None, g, folds, fit_score=fs, space=SPACE, n_trials=1, timeout_s=60, study_prefix="p",
                  storage=_st(tmp_path))
    assert len(got) == 4 * 3


def test_spaces_yaml_loads():
    for n in ("dmpnn", "molpool", "lgbm"):
        sp = hpo.load_space(n)
        assert sp["params"] and sp["fixed"]
    t = optuna.create_study().ask()
    p = hpo.suggest(t, hpo.load_space("dmpnn"))
    assert p["epochs"] == 1000 and p["bs"] in (32, 64, 128)
