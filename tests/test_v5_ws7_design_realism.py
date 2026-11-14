"""delta_ws7_design_realism (WS7): feature-off identity (R-W7-4), blendstock expansion, D86 domain, margin fit, optimiser, counterfactual."""
import importlib.util
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.design import blendstock_v5 as B
from fuelblend.design import engine_v4 as E
from fuelblend.design import library_route as LR

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "data" / "ws7_feature_off_fixture.json"
MARGINS = ROOT / "configs" / "uncertainty_v5.yaml"
warnings.simplefilter("ignore")

_s = importlib.util.spec_from_file_location("ws7fix", ROOT / "scripts/v5/base/delta_ws7_fixture.py")
FX = importlib.util.module_from_spec(_s)
_s.loader.exec_module(FX)

_r = importlib.util.spec_from_file_location("ws7run", ROOT / "scripts/v5/base/delta_ws7_design_realism.py")
RUN = importlib.util.module_from_spec(_r)
_r.loader.exec_module(RUN)


def _same(a, b):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= 1e-12 * max(1.0, abs(a))


# ------------------------------------------------------------------------------------------------ R-W7-4 feature off
def test_specs_unchanged_with_d86_off():
    fx = json.load(open(FIX, encoding="utf-8"))["specs"]
    now = FX.spec_snapshot()
    assert set(now) == set(fx)
    for k in fx:
        assert now[k] == fx[k], k


def test_d86_flag_adds_four_unverified_constraints():
    for sp in ("astm:jet_a", "bis:jet_a1"):
        c = E.load_spec(sp, d86=True)["constraints"]
        for k, lim in (("distillation_T10_C", ("max", 205.0)), ("distillation_FBP_C", ("max", 300.0)),
                       ("distillation_T50_minus_T10_C", ("min", 15.0)), ("distillation_T90_minus_T10_C", ("min", 40.0))):
            assert c[k][lim[0]] == lim[1] and c[k]["verification"] == "secondary_sources_only"
            assert E.canonical_property(k)[0] is not None
        assert not any(k.startswith("distillation_") for k in E.load_spec(sp)["constraints"])


def test_pool_c_recipes_identical_with_feature_off():
    fx = json.load(open(FIX, encoding="utf-8"))["pool_c"]
    now = FX.pool_c_snapshot()
    for spec in fx:
        assert set(now[spec]) == set(fx[spec])
        for k, a in fx[spec].items():
            b = now[spec][k]
            assert a["status"] == b["status"] and a["reason"] == b["reason"], (spec, k)
            for f in ("prediction", "slack", "margin"):
                assert _same(a[f], b[f]), (spec, k, f, a[f], b[f])


# ------------------------------------------------------------------------------------------------ blendstock expansion
def test_single_blendstock_equals_rep_level_recipe():
    bp = B.BlendPool([B.load_blendstock("atj_posf11498")])
    models = B.build_models(bp)
    w = bp.expand([1.0])
    r1 = B.judge("astm:jet_a", models, w, MARGINS, True, True, True, d86=False, represented_mass=0.0)
    r2 = LR.evaluate_recipe("astm:jet_a", models, dict(zip(models["pool"], w)), margins_path=MARGINS, synthetic_blend=True,
                            alt_limits=True, recipe_margins=True)
    assert [c.key for c in r1.constraints] == [c.key for c in r2.constraints]
    for a, b in zip(r1.constraints, r2.constraints):
        assert a.status == b.status and a.reason == b.reason
        for f in ("prediction", "slack", "applied_margin"):
            x, y = getattr(a, f), getattr(b, f)
            assert (x is None and y is None) or (isinstance(x, float) and math.isnan(x) and math.isnan(y)) or x == y


def test_blendpool_merges_same_compound_and_prefers_table_name():
    pure = B.pure_component("n-decane")
    bs = B.Blendstock(name="x", label="x", source="t", fuel_id=None, smiles=["CCCCCCCCCC", "CCCCCCCCCCCC"], w=np.array([0.4, 0.6]))
    bp = B.BlendPool([bs, pure])
    assert bp.comp_ids[0] == "n-decane" and len(bp.comp_ids) == 2
    np.testing.assert_allclose(bp.expand([0.5, 0.5]), [0.5 * 0.4 + 0.5, 0.5 * 0.6])
    assert bp.represented_mass([0.5, 0.5]) == pytest.approx(0.5)


def test_represented_rule_downgrades_pass_but_not_construction_or_fail():
    mk = lambda k, st, prop: E.ConstraintResult(key=k, prop=prop, lo=None, hi=1.0, model="m", prediction=0.5, applied_margin=0.1,
                                                margin_keys=(), status=st, reason="", slack=0.4)
    res = E.DesignResult(weights={}, predictions={}, margins={}, passed=True, extrapolating=False, binding_constraints=[],
                         constraints=[mk("a", "pass", "density"), mk("b", "by_construction", "aromatics_wt"), mk("c", "fail", "flash")])
    B.apply_represented_rule(res, 0.3)
    st = [c.status for c in res.constraints]
    assert st == ["uncertified", "by_construction", "fail"] and not res.passed
    assert "represented_blendstock_margin_not_validated" in res.constraints[0].reason


# ------------------------------------------------------------------------------------------------ D86
def test_d86_narrow_boiling_is_outside_domain_never_pass():
    q = B.d86_quantities([172.0, 179.0, 169.0], [0.3, 0.6, 0.1])
    assert not q["monotone"] or q["T90_T10"] < 40
    dom = dict(tbp_spread_min=40.0, tbp_spread_max=120.0)
    assert B.domain_distance(q, dom) > 0


def test_tb_from_registry_antoine_close_to_literature():
    assert abs(B.tb_C("CCCCCCCCCC") - 174.1) < 5.0
    assert abs(B.tb_C("CC(C)CC(C)(C)C") - 99.2) < 1.0


def test_conformal_margin_rank_and_floor():
    r = np.arange(1, 20, dtype=float)          # n = 19 -> rank ceil(20 * 0.9) = 18
    assert B.conformal_margin(r) == 18.0
    assert B.conformal_margin(-np.arange(1, 20, dtype=float)) == 0.0
    assert math.isinf(B.conformal_margin([1.0, 2.0, 3.0]))


def test_d86_margin_gate_and_certifying_max():
    rng = np.random.default_rng(0)

    def frame(n_av, n_aj):
        rows = []
        for s, n in (("av3322", n_av), ("ajftd", n_aj)):
            for _ in range(n):
                p = rng.normal(180, 5)
                rows.append(dict(source=s, **{f"pred_{q}": p for q in ("T10", "T50", "T90", "FBP", "T50_T10", "T90_T10")},
                                 **{f"meas_{q}": p + rng.normal(3 if s == "ajftd" else 0, 2) for q in ("T10", "T50", "T90", "FBP", "T50_T10", "T90_T10")}))
        return pd.DataFrame(rows)
    f = B.fit_d86_margins(frame(25, 5))
    assert all(math.isinf(v["certifying_margin"]) for v in f.values())
    f = B.fit_d86_margins(frame(25, 20))
    v = f["T10"]
    assert v["gate_passed"] and v["certifying_margin"] == max([v["pooled"]] + list(v["per_source"].values()))
    assert f["T50_T10"]["direction"] == "over" and f["T10"]["direction"] == "under"


# ------------------------------------------------------------------------------------------------ optimiser
def test_optimiser_respects_hard_cap_and_maximises():
    # two members; objective prefers member 1, a hard cap keeps member 1 <= 0.3
    vec = lambda u: np.array([u[1] - 0.1, 1.0])
    out = B.optimise(vec, 2, n_starts=5, seed=42, hard=lambda u: np.array([0.3 - u[1]]))
    assert out["feasible"] and out["u"][1] == pytest.approx(0.3, abs=1e-4)


def test_optimiser_uncertainty_penalty_moves_optimum():
    vec = lambda u: np.array([u[0]])
    a = B.optimise(vec, 2, n_starts=3, seed=42)
    b = B.optimise(vec, 2, n_starts=3, seed=42, penalty=lambda u: 4.0 * u[0] ** 2, lam=1.0)
    assert a["u"][0] == pytest.approx(1.0, abs=1e-4) and b["u"][0] == pytest.approx(0.125, abs=1e-3)


# ------------------------------------------------------------------------------------------------ counterfactual port
def test_counterfactual_removes_only_spec_flag_clauses():
    c = pd.DataFrame([dict(constraint="a", status="uncertified", slack=0.1,
                           reason="spec limit verification: secondary_sources_only (an unverified limit is never certified)"),
                      dict(constraint="b", status="uncertified", slack=0.1, reason=B.REPRESENTED_REASON)]).set_index("constraint")
    ok, rem = RUN.counterfactual(c.iloc[:1])
    assert ok and rem == []
    ok, rem = RUN.counterfactual(c)
    assert not ok and "represented_blendstock" in rem[0]
