"""v5 delta_smoke_cert 1.2: D1655 alternative smoke limit (>= 18 mm AND naphthalenes <= 3 vol %), flag use_alt_limits."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from fuelblend.design import engine_v4 as E

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "tests" / "data" / "smoke_alt_baseline.json"
MKEY = "smoke.test.key"
MARGINS = {MKEY: {"applied_margin": 3.0, "unit": "mm", "direction": "over"}}
POOL = ["a", "b"]
W = np.array([0.5, 0.5])


def _spec(use_alt=True, verification=None):
    sp = {"smoke_point_mm": {"min": 25.0, "alt_min": 18.0, "alt_requires": {"naphthalenes_vol_pct_max": 3.0}}}
    if verification:
        sp["smoke_point_mm"]["verification"] = verification
    return {"_family": "jet_a", "use_alt_limits": use_alt, "constraints": sp}


def _models(smoke, naph=None, naph_unavailable=False, naph_present=True):
    def smoke_fn(pool, w, T):
        return smoke

    def naph_fn(pool, w, T):
        if naph_unavailable:
            raise E.ModelUnavailable("no class/ring data")
        return naph

    props = {"smoke_point": E.PropertyModel("t", smoke_fn, margin_keys=[MKEY])}
    if naph_present:
        props["naphthalenes_vol_pct"] = E.PropertyModel("naph", naph_fn, by_construction=True)
    return {"props": props, "classes": {}, "validated_classes": None}


def _smoke(spec, models):
    return E.evaluate_blend(spec, POOL, W, models, MARGINS)[0]


def test_alt_limit_used():
    r = _smoke(_spec(), _models(22.0, 0.0))
    assert r.status == E.STATUS_PASS and "alt_min" in r.reason


def test_alt_limit_naph_unknown():
    assert _smoke(_spec(), _models(22.0, naph_unavailable=True)).status == E.STATUS_UNCERTIFIED
    assert _smoke(_spec(), _models(22.0, naph_present=False)).status == E.STATUS_UNCERTIFIED


def test_alt_limit_naph_over():
    assert _smoke(_spec(), _models(22.0, 3.5)).status == E.STATUS_FAIL


def test_alt_limit_primary_pass():
    r = _smoke(_spec(), _models(30.0, 9.9))
    assert r.status == E.STATUS_PASS and "primary" in r.reason and "alt_min" not in r.reason


def test_alt_limit_alt_margin_fails():
    # 20 - 3 = 17 < 18: alt limit fails with margin; stays fail
    assert _smoke(_spec(), _models(20.0, 0.0)).status == E.STATUS_FAIL


def test_alt_limit_flag_off_ignores_alt():
    r = _smoke(_spec(use_alt=False), _models(22.0, 0.0))
    assert r.status == E.STATUS_FAIL and "alt_min" not in r.reason


def test_alt_limit_unverified_spec_never_passes():
    r = _smoke(_spec(verification="secondary_sources_only"), _models(22.0, 0.0))
    assert r.status != E.STATUS_PASS and "limit used: alt_min" not in r.reason


def test_naphthalenes_model_class_fraction():
    from fuelblend.design import models_v4 as M
    L = M.library()
    mod = M.build_models(L)["props"]["naphthalenes_vol_pct"]
    assert mod.by_construction
    pool = ["1-methylnaphthalene", "n-decane"]
    v = mod.fn(pool, np.array([0.5, 0.5]), None)
    assert 30.0 < v < 70.0
    assert mod.fn(["tetralin", "n-decane"], np.array([0.5, 0.5]), None) == 0.0   # one aromatic ring only
    assert mod.fn(["methylcyclohexane", "n-decane"], np.array([0.5, 0.5]), None) == 0.0


def test_alt_limit_off_unchanged():
    """Flag off: the v5pool recipes give exactly the results frozen before the change (tests/data/gen_smoke_alt_baseline.py)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("gen_base", ROOT / "tests" / "data" / "gen_smoke_alt_baseline.py")
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    now = json.loads(json.dumps(g.snapshot(), default=str))
    old = json.loads(BASE.read_text())
    assert now.keys() == old.keys()
    for did in old:
        assert now[did].keys() == old[did].keys()
        for k in old[did]:
            a, b = now[did][k], old[did][k]
            assert a["status"] == b["status"] and a["reason"] == b["reason"], (did, k)
            for f in ("pred", "slack", "margin"):
                x, y = a[f], b[f]
                assert (x is None and y is None) or (str(x) == "nan" and str(y) == "nan") or x == y, (did, k, f)


def test_alt_pass_record_self_consistent():
    r = _smoke(_spec(), _models(22.0, 0.0))
    assert r.status == E.STATUS_PASS and r.lo == 18.0 and r.applied_margin == 3.0
    assert r.slack == pytest.approx((r.prediction - r.applied_margin - r.lo) / E._scale(r.lo, None))


def test_naphthalenes_missing_or_nan_smiles_unavailable():
    from fuelblend.design import models_v4 as M
    L = M.library()
    fn = M.build_models(L)["props"]["naphthalenes_vol_pct"].fn
    w = np.array([0.5, 0.5])
    with pytest.raises(E.ModelUnavailable):
        fn(["n-decane", "not-in-library"], w, None)
    L.smiles["n-undecane"] = float("nan")
    try:
        with pytest.raises(E.ModelUnavailable):
            fn(["n-decane", "n-undecane"], w, None)
    finally:
        L.smiles["n-undecane"] = "CCCCCCCCCCC"
