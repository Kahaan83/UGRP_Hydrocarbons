"""v4 Task 17: design engine contract tests (specs from yaml only, applied_margin only, 'inf' -> uncertified, ...)."""
import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.design import engine_v4 as E  # noqa: E402

POOL = ["A", "B"]                       # A = n-paraffin, B = aromatic
CLASSES = {"A": "n_paraffin", "B": "aromatic"}


def _write(tmp, spec, margins):
    cdir = Path(tmp)
    (cdir / "specs_test.yaml").write_text(yaml.safe_dump({"families": {"fam": {"constraints": spec}}}))
    mp = cdir / "unc.yaml"
    mp.write_text(yaml.safe_dump({"margins": margins}))
    return cdir, mp


def _pm(name, fn, keys=(), **kw):
    return E.PropertyModel(name=name, fn=fn, margin_keys=keys, **kw)


def _models(props, alternates=None, validated=None):
    return {"props": props, "alternates": alternates or {}, "classes": CLASSES, "validated_classes": validated}


def _design(cdir, mp, models):
    return E.design("test:fam", POOL, models, n_starts=4, seed=1, margins_path=mp, configs_dir=cdir)


# ---------------------------------------------------------------------------------------------------- extrapolation
def test_design_flags_extrapolation(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}},
                      {"k": {"applied_margin": 1.0, "unit": "C", "direction": "over"}})
    # flash model validated only on n-paraffins; the pool's second member is an aromatic
    fn = lambda pool, w, T: 50.0 * float(w[1])                    # only aromatic-rich blends reach the limit
    models = _models({"flash": _pm("m", fn, ["k"], validated_classes={"n_paraffin"})})
    res = _design(cdir, mp, models)
    assert res.weights["B"] > E.EXTRAP_WT_THRESHOLD       # optimiser pushes into the absent class
    assert res.extrapolating is True
    assert res.passed is False
    assert any(b["constraint"] == "extrapolation" for b in res.binding_constraints)
    # the same design in-domain passes
    models2 = _models({"flash": _pm("m", fn, ["k"], validated_classes={"n_paraffin", "aromatic"})})
    res2 = _design(cdir, mp, models2)
    assert res2.extrapolating is False and res2.passed is True


def test_pool_with_no_validated_class_is_forced_extrapolation(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}},
                      {"k": {"applied_margin": 1.0, "unit": "C", "direction": "over"}})
    models = _models({"flash": _pm("m", lambda p, w, T: 90.0, ["k"], validated_classes={"cycloalkane"})})
    res = _design(cdir, mp, models)
    assert res.extrapolating and not res.passed


# ---------------------------------------------------------------------------------------------------- AST: no spec literals
def _spec_numbers():
    nums = set()
    for f in ("specs_astm.yaml", "specs_bis.yaml"):
        doc = yaml.safe_load((ROOT / "configs" / f).read_text(encoding="utf-8"))
        for fam in doc["families"].values():
            for c in fam.get("constraints", {}).values():
                for k, v in c.items():
                    if k in ("min", "max") or k.startswith(("alt_", "lo", "hi")):
                        if isinstance(v, (int, float)) and not isinstance(v, bool):
                            nums.add(float(v))
    return nums


def test_specs_read_from_yaml_only():
    """No numeric literal in engine_v4.py equals a spec threshold from configs/specs_*.yaml.

    Scope (R7): spec thresholds only.  Generic constants are allowed: 0/1/2, and values bound to a module-level
    UPPER_CASE name (named engineering constants such as EXTRAP_WT_THRESHOLD)."""
    nums = _spec_numbers()
    assert len(nums) > 20                                  # the yaml really was parsed
    tree = ast.parse((ROOT / "src/fuelblend/design/engine_v4.py").read_text(encoding="utf-8"))
    named = set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) and t.id.isupper() for t in node.targets):
            named.update(id(n) for n in ast.walk(node.value))
    offenders = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)) and not isinstance(n.value, bool):
            if id(n) in named or n.value in (0, 1, 2):
                continue
            if abs(float(n.value)) in nums:
                offenders.append((n.lineno, n.value))
    assert not offenders, f"spec-threshold literals hard-coded in engine_v4.py: {offenders}"
    spec = E.load_spec("astm:gasoline")
    assert spec["constraints"]["dvpe_kPa"]["min"] in nums


def test_uncertainty_margin_keys_in_spec_are_ignored(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0, "uncertainty_margin": 99.0}},
                      {"k": {"applied_margin": 2.0, "unit": "C", "direction": "over"}})
    models = _models({"flash": _pm("m", lambda p, w, T: 42.0, ["k"])})
    cr = _design(cdir, mp, models).constraints[0]
    assert cr.applied_margin == 2.0 and cr.status == E.STATUS_PASS      # 42-2=40 >= 38; the hand-typed 99 is unused


# ---------------------------------------------------------------------------------------------------- volatility minimum
def test_volatility_minimum_enforced(tmp_path):
    cdir, mp = _write(tmp_path, {"dvpe_kPa": {"min": 45.0, "max": 70.0}}, {})
    calls = []

    def low(pool, w, T):
        calls.append(1)
        return 30.0                                   # below the yaml minimum whatever the blend

    models = _models({"dvpe": _pm("dvpe_raoult", low, margin_missing_reason="DVPE units unconfirmed, excluded")})
    res = _design(cdir, mp, models)
    cr = res.constraints[0]
    assert calls                                       # model really evaluated
    assert cr.status == E.STATUS_FAIL                  # below the minimum -> FAIL even though no margin exists
    assert res.passed is False
    # inside the window the pass is still never certified (no margin exists)
    ok = _models({"dvpe": _pm("dvpe_raoult", lambda p, w, T: 55.0, margin_missing_reason="DVPE units unconfirmed")})
    res2 = _design(cdir, mp, ok)
    assert res2.constraints[0].status == E.STATUS_UNCERTIFIED and res2.passed is False
    assert "units unconfirmed" in res2.constraints[0].reason


def test_real_gasoline_spec_has_dvpe_minimum_in_yaml():
    assert E.load_spec("astm:gasoline")["constraints"]["dvpe_kPa"].get("min") is not None


# ---------------------------------------------------------------------------------------------------- robustness
def test_robustness_runs_each_model(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}},
                      {"k": {"applied_margin": 1.0, "unit": "C", "direction": "over"}})
    counts = {"p": 0, "a1": 0, "a2": 0}

    def mk(tag, val):
        def fn(pool, w, T):
            counts[tag] += 1
            return val
        return fn

    models = _models({"flash": _pm("primary", mk("p", 60.0), ["k"])},
                     alternates={"flash": [_pm("alt1", mk("a1", 41.0), ["k"]), _pm("alt2", mk("a2", 55.0))]})
    rows = E.robustness_table("test:fam", POOL, [0.5, 0.5], models, margins_path=mp, configs_dir=cdir)
    assert counts == {"p": 1, "a1": 1, "a2": 1}                         # each model function called exactly once
    by = {r["model"]: r for r in rows}
    assert set(by) == {"primary", "alt1", "alt2"}
    assert by["primary"]["prediction"] == 60.0 and by["alt1"]["prediction"] == 41.0   # raw predictions, no offsets
    assert by["alt1"]["status"] == E.STATUS_PASS                        # 41 - 1 = 40 >= 38
    assert by["alt2"]["status"] == E.STATUS_UNCERTIFIED                 # no margin for this form: none borrowed


# ---------------------------------------------------------------------------------------------------- margins
def test_inf_margin_yields_not_certified(tmp_path):
    cdir, mp = _write(tmp_path, {"freezing_point_C": {"max": -40.0}},
                      {"svl": {"value": "inf", "unit": "C", "direction": "under"}})   # generator writes no applied_margin
    models = _models({"freezing": _pm("svl", lambda p, w, T: -60.0, ["svl"])})
    res = _design(cdir, mp, models)
    cr = res.constraints[0]
    assert cr.status == E.STATUS_UNCERTIFIED and cr.applied_margin is None
    assert "no valid margin" in cr.reason and "svl" in cr.reason
    assert res.passed is False
    assert cr.prediction == -60.0                       # nominal value is shown but not certified
    for bad in ("inf", float("inf"), None):
        assert E.read_applied_margin({"applied_margin": bad, "unit": "C"})[0] is None


def test_applied_margin_read_not_raw_quantile(tmp_path):
    # raw conformal quantile -6.5 (negative: would ADD headroom), applied (floored) 0.0
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}},
                      {"d93": {"value": -6.5, "applied_margin": 0.0, "unit": "C", "direction": "over"}})
    models = _models({"flash": _pm("m", lambda p, w, T: 40.0, ["d93"])})
    cr = _design(cdir, mp, models).constraints[0]
    assert cr.applied_margin == 0.0
    assert cr.status == E.STATUS_PASS and cr.slack == pytest.approx(2.0 / 38.0)   # 40 - 0 vs 38, NOT 40 + 6.5
    assert E.read_applied_margin({"value": 3.0})[0] is None                       # value alone is never used


def test_flash_margin_is_larger_of_lab_and_d93(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}},
                      {"lab": {"applied_margin": 3.7, "unit": "C", "direction": "over"},
                       "d93": {"applied_margin": 0.0, "unit": "C", "direction": "over"}})
    models = _models({"flash": _pm("m", lambda p, w, T: 50.0, ["lab", "d93"])})
    assert _design(cdir, mp, models).constraints[0].applied_margin == 3.7


def test_wrong_direction_margin_is_uncertified(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}},
                      {"k": {"applied_margin": 3.0, "unit": "C", "direction": "under"}})   # protects max limits only
    models = _models({"flash": _pm("m", lambda p, w, T: 60.0, ["k"])})
    assert _design(cdir, mp, models).constraints[0].status == E.STATUS_UNCERTIFIED


def test_model_unavailable_is_uncertified_not_exception(tmp_path):
    cdir, mp = _write(tmp_path, {"flash_point_C": {"min": 38.0}}, {})

    def boom(pool, w, T):
        raise E.ModelUnavailable("unverified constant")

    cr = _design(cdir, mp, _models({"flash": _pm("m", boom)})).constraints[0]
    assert cr.status == E.STATUS_UNCERTIFIED and "unverified constant" in cr.reason


# ---------------------------------------------------------------------------------------------------- fix round 1
def test_margin_validated_on_jet_only_is_not_certified_for_diesel(tmp_path):
    (tmp_path / "specs_test.yaml").write_text(yaml.safe_dump({"families": {
        "jet_x": {"constraints": {"flash_point_C": {"min": 38.0}}},
        "diesel_x": {"constraints": {"flash_point_C": {"min": 38.0}}}}}))
    mp = tmp_path / "unc.yaml"
    mp.write_text(yaml.safe_dump({"margins": {"k": {"applied_margin": 1.0, "unit": "C", "direction": "over"}}}))
    pm = _pm("m", lambda p, w, T: 60.0, ["k"])
    pm.certified_kinds = {"jet"}
    models = _models({"flash": pm})
    jet = E.design("test:jet_x", POOL, models, n_starts=2, seed=1, margins_path=mp, configs_dir=tmp_path)
    die = E.design("test:diesel_x", POOL, models, n_starts=2, seed=1, margins_path=mp, configs_dir=tmp_path)
    assert jet.passed and jet.constraints[0].status == E.STATUS_PASS
    assert die.constraints[0].status == E.STATUS_SCOPE and die.passed is False
    assert E.STATUS_SCOPE not in E.OK_STATUSES


@pytest.mark.parametrize("flag", [{"verification": "UNVERIFIED"}, {"verification": "secondary_sources_only"},
                                  {"unverified": True}])
def test_unverified_spec_limit_never_certified(tmp_path, flag):
    cdir, mp = _write(tmp_path, {"flash_point_C": dict({"min": 38.0}, **flag),
                                 "aromatics_wt_pct": dict({"max": 90.0}, **flag)},
                      {"k": {"applied_margin": 1.0, "unit": "C", "direction": "over"}})
    models = _models({"flash": _pm("m", lambda p, w, T: 60.0, ["k"]),
                      "aromatics_wt": _pm("c", lambda p, w, T: 10.0, by_construction=True)})
    res = _design(cdir, mp, models)
    assert [c.status for c in res.constraints] == [E.STATUS_UNCERTIFIED, E.STATUS_UNCERTIFIED]
    assert all("never certified" in c.reason for c in res.constraints) and res.passed is False


def test_ridge_envelope_guard_fires_for_pure_decane():
    from fuelblend.design import models_v4 as M
    L = M.library()
    fns = M.make_fns(L)
    w = np.array([1.0])
    with pytest.raises(E.ModelUnavailable, match="envelope"):
        fns["visc_ridge"](["n-decane"], w, -20.0)
    unguarded = fns["visc_ridge_unguarded"](["n-decane"], w, -20.0)
    assert unguarded < 0.5                              # the nonphysical value the guard prevents (0.126 cSt)
    assert fns["visc_recal"](["n-decane"], w, -20.0) > 2.0
    chk = pd.read_csv(ROOT / "results/v4/D/viscosity_envelope_check_v4.csv")
    row = chk[(chk.pool == '["n-decane"]') & (chk.T_C == -20.0)].iloc[0]
    assert row.ridge_unguarded_cSt == pytest.approx(unguarded, rel=1e-6) and np.isnan(row.ridge_guarded_cSt)
