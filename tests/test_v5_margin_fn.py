"""v5 delta_smoke_cert 1.5: recipe-dependent margins (PropertyModel.margin_fn) in the optimiser, alt-aware _Sides/rejudge."""
import dataclasses
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.design import engine_v4 as E  # noqa: E402
from fuelblend.design import library_route as LR  # noqa: E402
from fuelblend.design import smoke_margin_v5 as SM  # noqa: E402

REAL_MARGINS = ROOT / "configs" / "uncertainty_v5.yaml"
V5_SMOKE_KEY = "smoke_point_mm.ajftd_v5_in_domain.bulk_regression"


def _margins_without_v5_smoke_key():
    """Pre-change margins file: the live yaml minus the real-fuel smoke key, so these tests keep exercising the lab-key (7.8 mm) path."""
    import tempfile
    doc = yaml.safe_load(REAL_MARGINS.read_text(encoding="utf-8"))
    doc["margins"].pop(V5_SMOKE_KEY, None)
    f = Path(tempfile.mkdtemp()) / "uncertainty_no_v5_smoke.yaml"
    f.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
    return f


MARGINS = _margins_without_v5_smoke_key()
DESIGNS = ROOT / "results" / "v5" / "BASE" / "delta_design_ml_visc" / "designs"
FOLDERS = ["ext_n30", "ext_n8", "guard4_v5pool", "v5pool"]
BASELINE = ROOT / "tests" / "data" / "margin_fn_baseline.json"


def _cases():
    out = []
    for f in FOLDERS:
        d = pd.read_csv(DESIGNS / f / "designs.csv")
        for r in d.itertuples():
            out.append(pytest.param(f, r.spec_name, json.loads(r.pool), json.loads(r.weights), id=f"{f}-{r.design_id}"))
    return out


_MODELS: dict = {}


def _models_for(pool):
    k = tuple(pool)
    if k not in _MODELS:
        _MODELS[k] = LR.build_route_models(list(pool))
    return _MODELS[k]


def _num(x):
    return float("nan") if x is None else float(x)


def _same(a, b, tol=1e-9):
    a, b = _num(a), _num(b)
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= tol


def _by_key(res):
    return {c.key: c for c in res.constraints}


# ------------------------------------------------------------------------------------------------ stubs (no data, fast)
POOL3 = ["toluene", "n_dodecane", "iso_octane"]
STUB_CLASSES = {"toluene": "aromatic", "n_dodecane": "n_paraffin", "iso_octane": "iso_paraffin"}


def _stub_cfg(tmp, arom_min=10.0, with_alt=True, verification=None):
    smoke = {"min": 25.0}
    if verification:
        smoke["verification"] = verification
    if with_alt:
        smoke.update(alt_min=18.0, alt_requires={"naphthalenes_vol_pct_max": 3.0})
    spec = {"aromatics_vol_pct": {"min": arom_min}, "smoke_point_mm": smoke}
    cdir = Path(tmp)
    (cdir / "specs_test.yaml").write_text(yaml.safe_dump({"families": {"jet_stub": {"constraints": spec}}}))
    mp = cdir / "unc.yaml"
    mp.write_text(yaml.safe_dump({"margins": {}}))
    return cdir, mp


def _stub_models(margin_fn=None, margin_keys=(), naph_by_construction=True):
    extra = {} if margin_fn is None else dict(margin_fn=margin_fn, margin_direction="over")
    smoke = E.PropertyModel("stub_smoke", lambda pool, w, T: 28.0 - 30.0 * float(w[0]), margin_keys=margin_keys, **extra)
    return {"props": {"aromatics_vol": E.PropertyModel("arom", lambda pool, w, T: 100.0 * float(w[0]), by_construction=True),
                      "naphthalenes_vol_pct": E.PropertyModel("naph", lambda pool, w, T: 0.0, by_construction=naph_by_construction),
                      "smoke_point": smoke},
            "alternates": {}, "classes": STUB_CLASSES, "validated_classes": None}


def _in_domain_3mm(pool, w, T_C):
    return 3.0 if float(w[0]) <= 0.6 else None


def _stub_design(tmp_path, models, alt=True, **kw):
    cdir, mp = _stub_cfg(tmp_path)
    return E.design("test:jet_stub", POOL3, models, n_starts=5, seed=42, margins_path=mp, configs_dir=cdir,
                    alt_limits=alt, **kw)


# ------------------------------------------------------------------------------------------------ tests
def test_optimiser_sees_smoke_margin(tmp_path):
    res = _stub_design(tmp_path, _stub_models(margin_fn=_in_domain_3mm))
    cr = {c.key: c for c in res.constraints}
    sm = cr["smoke_point_mm"]
    assert sm.applied_margin == 3.0
    assert sm.status == E.STATUS_PASS, sm.reason
    assert sm.prediction - 3.0 >= 18.0 - 1e-6                      # passes the 18 mm alternative limit WITH the margin
    assert cr["aromatics_vol_pct"].prediction >= 10.0 - 1e-6        # and the aromatics minimum
    assert res.passed
    # the alt-blind optimiser (alt flag off) cannot satisfy both: smoke primary 25 mm is unreachable at >= 10 vol % aromatics
    res_off = _stub_design(tmp_path, _stub_models(margin_fn=_in_domain_3mm), alt=False)
    assert not res_off.passed


def test_alt_met_on_unverified_spec_is_uncertified_never_fail(tmp_path):
    """Review fix A: a recipe that misses the primary 25 mm limit but meets the 18 mm alternative (with margin and naphthalenes) on a spec
    whose verification is 'secondary_sources_only' is UNCERTIFIED on every path (margin_fn, yaml margin), never FAIL; the same recipe on
    a verified spec passes; a recipe missing the alternative limit with its margin still fails."""
    w_met = np.array([0.2, 0.4, 0.4])          # smoke 22 mm: primary miss, alt 18 met with 3 mm margin (19)
    w_miss = np.array([0.4, 0.3, 0.3])         # smoke 16 mm: alt 18 missed even nominally
    w_marg = np.array([0.3, 0.35, 0.35])       # smoke 19 mm: alt met nominally, 19 - 3 = 16 < 18 with the margin
    for verification, want in (("secondary_sources_only", E.STATUS_UNCERTIFIED), ("verified", E.STATUS_PASS)):
        cdir, _ = _stub_cfg(tmp_path, verification=verification)
        spec = E.load_spec("test:jet_stub", cdir, alt_limits=True)
        mp = Path(tmp_path) / "k.yaml"
        mp.write_text(yaml.safe_dump({"margins": {"k": {"applied_margin": 3.0, "unit": "mm", "direction": "over"}}}))
        for models in (_stub_models(margin_fn=lambda pool, w, T: 3.0), _stub_models(margin_keys=["k"])):
            marg = E.load_margins(mp)
            def st(w):
                return {x.key: x for x in E.evaluate_blend(spec, POOL3, w, models, marg)}["smoke_point_mm"]
            assert st(w_met).status == want, (verification, st(w_met).reason)
            assert st(w_marg).status == E.STATUS_FAIL and st(w_miss).status == E.STATUS_FAIL


def test_requirement_without_margin_is_infeasible(tmp_path):
    """Review fix B: a requirement whose model is not by_construction and has no margin is never treated as margin 0: the optimiser's
    alternative branch is dead (_Sides.alt None) and the engine does not pass the alternative; a by_construction requirement (naphthalenes)
    legitimately has margin 0 in both."""
    cdir, mp = _stub_cfg(tmp_path)
    marg = E.load_margins(mp)
    spec = E.load_spec("test:jet_stub", cdir, alt_limits=True)
    ok = E._Sides(spec, POOL3, _stub_models(margin_fn=_in_domain_3mm), marg)
    assert ok.alt["smoke_point_mm"] is not None
    bad_models = _stub_models(margin_fn=_in_domain_3mm, naph_by_construction=False)
    bad = E._Sides(spec, POOL3, bad_models, marg)
    assert bad.alt["smoke_point_mm"] is None
    w = np.array([0.2, 0.4, 0.4])
    sm = {x.key: x for x in E.evaluate_blend(spec, POOL3, w, bad_models, marg)}["smoke_point_mm"]
    assert sm.status not in E.OK_STATUSES
    good = {x.key: x for x in E.evaluate_blend(spec, POOL3, w, _stub_models(margin_fn=_in_domain_3mm), marg)}["smoke_point_mm"]
    assert good.status in E.OK_STATUSES                                         # no verification flag in the stub spec: the alternative passes
    assert ok.vector(w)[1] > 0 > bad.vector(w)[1]                                # smoke slack: alt branch alive (+) vs dead (primary miss)
    res = _stub_design(tmp_path, bad_models)
    assert not res.passed


def test_none_margin_is_infeasible(tmp_path):
    res = _stub_design(tmp_path, _stub_models(margin_fn=lambda pool, w, T: None))
    sm = {c.key: c for c in res.constraints}["smoke_point_mm"]
    assert sm.status not in E.OK_STATUSES and sm.applied_margin is None   # fail (nominal 25 mm miss) or uncertified, never pass
    assert not res.passed
    assert res.optimiser["best_min_slack"] <= E.BAD_SLACK + 1e-9      # infeasible, not "margin 0"
    zero = _stub_design(tmp_path, _stub_models(margin_fn=lambda pool, w, T: 0.0))
    assert zero.optimiser["best_min_slack"] > E.BAD_SLACK + 1.0      # a zero margin is NOT the same as no margin


def test_margin_fn_slack_matches_sides(tmp_path):
    """_Sides (optimiser) and evaluate_blend (judge) give one slack per recipe, alt on and off, to 1e-9."""
    cdir, mp = _stub_cfg(tmp_path)
    margins = E.load_margins(mp)
    models = _stub_models(margin_fn=_in_domain_3mm)
    for alt in (False, True):
        spec = E.load_spec("test:jet_stub", cdir, alt_limits=alt)
        sides = E._Sides(spec, POOL3, models, margins)
        for w in ([0.05, 0.5, 0.45], [0.119, 0.4, 0.481], [0.3, 0.3, 0.4], [0.7, 0.2, 0.1]):
            w = np.array(w)
            vec = sides.vector(w)
            cons = E.evaluate_blend(spec, POOL3, w, models, margins)
            by = {c.key: c for c in cons}
            # sides order follows spec constraint order: aromatics (1 side), smoke (1 side)
            assert _same(vec[0], by["aromatics_vol_pct"].slack)
            if w[0] <= 0.6:
                assert _same(vec[1], by["smoke_point_mm"].slack), (alt, w, vec, by["smoke_point_mm"].slack)
            else:
                assert vec[1] == E.BAD_SLACK                          # margin_fn None -> infeasible side


def test_class_mix_basis():
    """The recipe-side class fractions use the GCxGC basis: cycloalkane = mono+di+tri-cycloparaffin, one mass total.
    A real fuel's own class fractions, written as a recipe, pass the gate (for every mapped fuel)."""
    f = SM.class_fractions()
    comp = pd.read_csv(SM.COMP)
    fid = int(f.index[0])
    sub = comp[comp.fuel_id == fid]
    tot = sub.mass_pct.sum()
    cyc = sub[sub.hc_class.isin(["monocycloparaffin", "dicycloparaffin", "tricycloparaffin"])].mass_pct.sum() / tot
    assert abs(f.loc[fid, "cyclo"] - cyc) < 1e-12
    classes = {"c": "cycloalkane", "n": "n_paraffin", "i": "iso_paraffin", "a": "aromatic"}
    pool = ["c", "n", "i", "a"]
    for fid, r in f.iterrows():
        w = np.array([r.cyclo, r.nparaffin, r.isoparaffin, max(0.0, 1.0 - r.cyclo - r.nparaffin - r.isoparaffin)])
        fr = SM.recipe_class_fractions(pool, w * 3.7, classes)           # un-normalised weights are normalised
        assert np.allclose(fr, (r.cyclo, r.nparaffin, r.isoparaffin), atol=1e-9)
        ok, why = SM.class_mix_ok_for_recipe(pool, w, classes)
        assert ok and why is None, fid
    ok, why = SM.class_mix_ok_for_recipe(pool, np.array([0.75, 0.05, 0.1, 0.1]), classes)   # 75 % cycloalkane
    assert not ok and why == SM.REASON
    # the library's decalins and butylcyclohexane are cycloalkane (GCxGC dicyclo/monocyclo)
    from fuelblend.design import models_v4 as M
    L = M.library()
    assert {L.cls(n) for n in ("decalin", "decalin_trans", "butylcyclohexane")} == {"cycloalkane"}


class _FakeSmokeModels:
    """smoke margin lookup with an in-memory margins dict (the real yaml stays untouched)."""


def _smoke_fn_stub(pool, w, T):
    return 22.0


def test_smoke_margin_key_replaces_lab_key():
    classes = {"c": "cycloalkane", "n": "n_paraffin", "i": "iso_paraffin", "a": "aromatic"}
    pool = ["c", "n", "i", "a"]
    ok_w = np.array([0.2, 0.2, 0.4, 0.2])
    bad_w = np.array([0.8, 0.05, 0.05, 0.1])
    lab = {LR.SMOKE_MARGIN_KEY: {"applied_margin": 7.8, "unit": "mm", "direction": "over"}}
    both = dict(lab)
    both[LR.SMOKE_V5_KEY] = {"applied_margin": 4.11, "unit": "mm", "direction": "over"}
    in_dom = lambda: True
    f_lab = LR.smoke_recipe_margin(pool, ok_w, classes, lab, in_domain=True)
    f_new = LR.smoke_recipe_margin(pool, ok_w, classes, both, in_domain=True)
    assert f_lab == (7.8, LR.SMOKE_MARGIN_KEY, None)               # key absent: today's 7.8 mm, no class-mix gate
    assert f_new == (4.11, LR.SMOKE_V5_KEY, None)                  # key present REPLACES (not max over both)
    assert LR.smoke_recipe_margin(pool, bad_w, classes, both, in_domain=True)[0] is None
    assert LR.smoke_recipe_margin(pool, bad_w, classes, both, in_domain=True)[2] == SM.REASON
    assert LR.smoke_recipe_margin(pool, bad_w, classes, lab, in_domain=True)[0] == 7.8     # no gate without the new key
    assert LR.smoke_recipe_margin(pool, ok_w, classes, both, in_domain=False)[0] is None   # out of regression domain
    # present but unusable entry: no fallback to the lab key
    bad = dict(lab)
    bad[LR.SMOKE_V5_KEY] = {"value": "inf"}
    assert LR.smoke_recipe_margin(pool, ok_w, classes, bad, in_domain=True)[0] is None


@pytest.mark.parametrize("alt", [False, True])
@pytest.mark.parametrize("folder,spec,pool,weights", _cases())
def test_margin_fn_matches_rejudge(folder, spec, pool, weights, alt):
    models = _models_for(pool)
    w = np.array([weights.get(n, 0.0) for n in models["pool"]], float)
    w /= w.sum()
    old = LR.evaluate_recipe(spec, models, weights, margins_path=MARGINS, alt_limits=alt)
    new = LR.evaluate_recipe(spec, models, weights, margins_path=MARGINS, alt_limits=alt, recipe_margins=True)
    assert new.optimiser.get("margin_fn_consistency") == []
    # the raw engine with margin_fn (no post-hoc rejudge at all)
    marg = E.load_margins(MARGINS)
    sp = E.load_spec(spec, alt_limits=alt)
    raw = {c.key: c for c in E.evaluate_blend(sp, models["pool"], w, LR.with_recipe_margins(models, marg), marg)}
    o, n = _by_key(old), _by_key(new)
    for key in o:
        for tag, other in (("annotated", n[key]), ("raw", raw[key])):
            if o[key].prop in ("viscosity", "smoke_point"):
                assert other.status == o[key].status, (key, tag, other.reason, o[key].reason)
                assert _same(other.slack, o[key].slack), (key, tag)
                assert _same(other.applied_margin, o[key].applied_margin), (key, tag)
                assert _same(other.prediction, o[key].prediction), (key, tag)
            else:
                assert other.status == o[key].status and (tag == "raw" or other.reason == o[key].reason)
                assert _same(other.slack, o[key].slack) and _same(other.applied_margin, o[key].applied_margin)


def test_margin_fn_matches_rejudge_smoke_in_domain():
    """In-domain smoke (aromatics 8-23 vol %): lab-key margin through margin_fn == the post-hoc rejudge, alt on and off."""
    pool = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
    models = _models_for(pool)
    weights = {"n-decane": 0.2, "n-undecane": 0.2, "n-dodecane": 0.1, "butylcyclohexane": 0.3, "decalin_trans": 0.0,
               "1,2,4-trimethylbenzene": 0.2}
    marg = E.load_margins(MARGINS)
    for spec in ("astm:jet_a", "bis:jet_a1"):
        for alt in (False, True):
            old = LR.evaluate_recipe(spec, models, weights, margins_path=MARGINS, alt_limits=alt)
            new = LR.evaluate_recipe(spec, models, weights, margins_path=MARGINS, alt_limits=alt, recipe_margins=True)
            assert new.optimiser.get("margin_fn_consistency") == []
            o, n = _by_key(old)["smoke_point_mm"], _by_key(new)["smoke_point_mm"]
            assert o.applied_margin is not None and abs(o.applied_margin - 7.8009965628100595) < 1e-12
            assert _same(o.slack, n.slack) and _same(o.applied_margin, n.applied_margin)
            assert o.status == n.status, (spec, alt, o.status, n.status)     # review fix A: independent flags-off path vs margin_fn path, both ways


def test_margin_fn_absent_unchanged():
    """R-SC-4: with recipe_margins off (default) the results equal the frozen pre-change results."""
    frozen = json.loads(BASELINE.read_text())
    assert dataclasses.fields(E.PropertyModel)                      # imports fine
    import inspect
    assert inspect.signature(LR.design_routed).parameters["recipe_margins"].default is False
    assert inspect.signature(LR.evaluate_recipe).parameters["recipe_margins"].default is False
    for case in _cases():
        folder, spec, pool, weights = case.values
        key = case.id
        models = _models_for(pool)
        for alt in (False, True):
            res = LR.evaluate_recipe(spec, models, weights, margins_path=MARGINS, alt_limits=alt)
            want = frozen["evaluate_recipe"][f"{key}|alt={alt}"]
            got = {c.key: dict(status=c.status, pred=c.prediction, slack=c.slack, margin=c.applied_margin, reason=c.reason,
                               model=c.model, margin_keys=list(c.margin_keys)) for c in res.constraints}
            assert set(got) == set(want)
            for k in want:
                g, w_ = got[k], want[k]
                if alt and k == "smoke_point_mm":
                    # review fix A: the alt-branch wording changed ("met but not certified" instead of "not passed"); the status must not
                    # be worse than before (a recipe meeting the alt limit on an unverified spec is uncertified, never fail)
                    assert g["model"] == w_["model"], (key, k)
                    assert E._status_rank(g["status"]) >= E._status_rank(w_["status"]), (key, k, g["status"], w_["status"])
                else:
                    assert g["status"] == w_["status"] and g["reason"] == w_["reason"] and g["model"] == w_["model"], (key, k)
                assert g["margin_keys"] == w_["margin_keys"], (key, k)
                for f in ("pred", "slack", "margin"):
                    if alt and f == "slack" and k == "smoke_point_mm":
                        # Task 1.5 rule 4: with the alt flag ON the reported smoke slack is max(primary, alternative branch)
                        assert g[f] >= w_[f] - 1e-12, (key, k, f)
                        continue
                    assert _same(g[f], w_[f], 1e-12), (key, k, f)
    # robustness table, flag-free path
    pool = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
    models = _models_for(pool)
    w = [0.15, 0.13, 0.03, 0.69, 0.0, 0.0]
    rows = E.robustness_table("astm:jet_a", pool, w, models, margins_path=MARGINS)
    want = frozen["robustness"]
    assert [(r["constraint"], r["model"], r["status"]) for r in rows] == [tuple(x[:3]) for x in want]
    for r, x in zip(rows, want):
        assert _same(r["prediction"], x[3], 1e-12) and _same(r["applied_margin"], x[4], 1e-12)
    # stub design: weights and best slack unchanged (alt off, no margin_fn)
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        models = _stub_models(margin_keys=["k"])
        cdir, mp = _stub_cfg(tmp, arom_min=10.0, with_alt=False)
        mp.write_text(yaml.safe_dump({"margins": {"k": {"applied_margin": 2.0, "unit": "mm", "direction": "over"}}}))
        r = E.design("test:jet_stub", POOL3, models, n_starts=5, seed=42, margins_path=mp, configs_dir=cdir)
        st = frozen["stub"]
        assert np.allclose([r.weights[n] for n in POOL3], st["weights"], atol=1e-12)
        assert abs(r.optimiser["best_min_slack"] - st["best"]) < 1e-12


def test_predict_route_smoke_key_replaces_lab_key(monkeypatch):
    """User decision R-SC-2b, 'Yes, replace (Recommended)': predict_route._smoke with the real yaml (v5 key present) gives 4.107069860688682 mm
    for an in-domain recipe inside the class-mix range and an uncertified result (reason class_mix_outside_real_fuel_range) for an in-domain
    recipe outside it; with the key absent (in-memory config) the 7.8 mm lab key is used exactly as before, with no class-mix gate."""
    from fuelblend.api import predict_route as pr
    real = pr._v5_margins()
    assert pr.SMOKE_V5_KEY in real and real[pr.SMOKE_V5_KEY]["applied_margin"] == 4.107069860688682
    names = pd.read_csv(ROOT / "data" / "processed" / "v5" / "pure_props_v5.csv").drop_duplicates("compound").set_index("compound").smiles
    comp = {names[n]: x for n, x in {"n-decane": 0.2, "n-dodecane": 0.1, "isooctane": 0.2, "butylcyclohexane": 0.3,
                                     "1,2,4-trimethylbenzene": 0.2}.items()}
    comp_bad = {names[n]: x for n, x in {"n-decane": 0.05, "butylcyclohexane": 0.75, "1,2,4-trimethylbenzene": 0.2}.items()}   # no iso-alkane
    r1 = pr.predict_fuel(comp, T_K=293.15, fuel_type="jet")["smoke_point_mm"]
    assert r1["margin"] == 4.107069860688682 and r1["in_domain"] is True
    r2 = pr.predict_fuel(comp_bad, T_K=293.15, fuel_type="jet")["smoke_point_mm"]
    assert r2.get("margin") is None and r2.get("value") is not None and r2["uncertified"] is True
    assert any("class_mix_outside_real_fuel_range" in x for x in r2["domain_reasons"])
    lab_only = {k: v for k, v in real.items() if k != pr.SMOKE_V5_KEY}
    monkeypatch.setattr(pr, "_v5_margins", lambda: lab_only)
    r0 = pr.predict_fuel(comp, T_K=293.15, fuel_type="jet")["smoke_point_mm"]
    assert r0["margin"] == pytest.approx(7.8009965628100595)                     # key absent: today's behaviour
    r3 = pr.predict_fuel(comp_bad, T_K=293.15, fuel_type="jet")["smoke_point_mm"]
    assert r3["margin"] == pytest.approx(7.8009965628100595)                     # no class-mix gate without the new key


def test_engine_smoke_margin_uses_v5_key_in_memory():
    """with_recipe_margins on the real route models: v5 key present (in-memory dict) replaces 7.8 mm for an in-domain recipe inside the
    class-mix range; a butylcyclohexane-heavy recipe gets no margin (reason class_mix_outside_real_fuel_range)."""
    pool = ["n-decane", "n-dodecane", "isooctane", "butylcyclohexane", "1,2,4-trimethylbenzene"]
    models = _models_for(pool)
    marg = E.load_margins(MARGINS)
    withkey = {**marg, LR.SMOKE_V5_KEY: {"applied_margin": 4.11, "unit": "mm", "direction": "over"}}
    spec = E.load_spec("astm:jet_a")
    ok = np.array([0.2, 0.1, 0.2, 0.3, 0.2])
    c = {x.key: x for x in E.evaluate_blend(spec, models["pool"], ok, LR.with_recipe_margins(models, withkey), withkey)}["smoke_point_mm"]
    assert c.applied_margin == pytest.approx(4.11)
    c0 = {x.key: x for x in E.evaluate_blend(spec, models["pool"], ok, LR.with_recipe_margins(models, marg), marg)}["smoke_point_mm"]
    assert c0.applied_margin == pytest.approx(7.8009965628100595)
    bad = np.array([0.05, 0.0, 0.0, 0.75, 0.2])
    cb = {x.key: x for x in E.evaluate_blend(spec, models["pool"], bad, LR.with_recipe_margins(models, withkey), withkey)}["smoke_point_mm"]
    assert cb.applied_margin is None and SM.REASON in cb.reason and cb.status != E.STATUS_PASS


def test_recipe_margins_viscosity_label_is_ml_when_in_domain():
    """Task 1.6 check: with recipe_margins the in-domain viscosity row carries the ML label (value and margin are the ML ones,
    identical to the recipe_margins-off result); a Refutas label must never sit beside the ML margin."""
    pool = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
    models = LR.build_route_models(pool)
    w = {n: 1.0 / len(pool) for n in pool}
    rows = {}
    for rm in (False, True):
        r = LR.evaluate_recipe("astm:jet_a", models, w, margins_path=REAL_MARGINS, recipe_margins=rm)
        rows[rm] = {c.key: c for c in r.constraints}["kin_visc_m20C_cSt"]
    assert rows[True].model == LR.VISC_ML_NAME == rows[False].model
    assert rows[True].prediction == pytest.approx(rows[False].prediction, rel=1e-12)
    assert rows[True].applied_margin == pytest.approx(rows[False].applied_margin, rel=1e-12)


def _rejudge_stub(tmp_path, monkeypatch, engine_reason, engine_status=None):
    """Smoke row of a recipe that misses the 25 mm limit, meets the verified 18 mm alternative with the 3 mm lab margin, with an engine reason."""
    cdir, _ = _stub_cfg(tmp_path, verification="verified")
    monkeypatch.setattr(E, "CONFIGS", cdir)
    mp = Path(tmp_path) / "k2.yaml"
    mp.write_text(yaml.safe_dump({"margins": {LR.SMOKE_MARGIN_KEY: {"applied_margin": 3.0, "unit": "mm", "direction": "over"}}}))
    marg = E.load_margins(mp)
    w = np.array([0.2, 0.4, 0.4])                      # smoke 22 mm: 22 - 3 = 19 >= 18
    spec = E.load_spec("test:jet_stub", cdir, alt_limits=True)
    models = _stub_models()
    cons = E.evaluate_blend(spec, POOL3, w, models, marg)
    sm = {c.key: c for c in cons}["smoke_point_mm"]
    sm.reason = (sm.reason + "; " if sm.reason else "") + engine_reason
    if engine_status is not None:
        sm.status = engine_status
    res = E.DesignResult(weights=dict(zip(POOL3, map(float, w))), predictions={}, margins={}, passed=False, extrapolating=False,
                         binding_constraints=[], constraints=cons, spec_name="test:jet_stub", pool=list(POOL3), optimiser={})
    LR._rejudge_smoke(res, marg, models, alt_limits=True)
    return {c.key: c for c in res.constraints}["smoke_point_mm"]


def test_rejudge_smoke_alt_met_verified_keeps_other_engine_blockers(tmp_path, monkeypatch):
    """Review IMPORTANT-1: alt branch met on a VERIFIED spec passes only when the engine reported no other blocker."""
    clean = _rejudge_stub(tmp_path, monkeypatch, "")
    assert clean.status == E.STATUS_PASS, clean.reason
    for why in ("extrapolating: 0.50 wt fraction in classes absent from this model's validation data",
                "non-verified constant(s): x.rho15=disputed",
                "margin direction does not protect the side of the limit"):
        blocked = _rejudge_stub(tmp_path, monkeypatch, why)
        assert blocked.status == E.STATUS_UNCERTIFIED and blocked.status not in E.OK_STATUSES, (why, blocked.reason)


def _plateau_models():
    """Smoke margin 3 mm only for w[0] <= 0.15, else None (the centroid start sits on the None plateau)."""
    fn = lambda pool, w, T: 3.0 if float(w[0]) <= 0.15 else None            # noqa: E731
    m = _stub_models(margin_fn=fn)
    m["props"]["smoke_point"] = dataclasses.replace(
        m["props"]["smoke_point"], margin_stand_in=3.0, margin_penalty_fn=lambda pool, w, T: max(float(w[0]) - 0.15, 0.0))
    return m


def test_graded_objective_moves_off_the_plateau_judge_unchanged(tmp_path):
    """Review IMPORTANT-2: with a None margin the old objective is the constant BAD_SLACK (SLSQP stays at its start); the graded objective
    (optimiser only) has a gradient and moves; the judge still gives uncertified for a None margin, never margin 0."""
    cdir, mp = _stub_cfg(tmp_path)
    spec = E.load_spec("test:jet_stub", cdir, alt_limits=True)
    marg = E.load_margins(mp)
    models = _plateau_models()
    c = np.full(3, 1 / 3)
    old = E._Sides(spec, POOL3, models, marg)                    # graded off: bit-identical old behaviour
    assert (old.vector(c) <= E.BAD_SLACK).any() and old.vector(c).min() == E.BAD_SLACK
    w_old, s_old, _ = E._optimise(old, 3, 1, 42)
    assert np.allclose(w_old, c) and s_old == E.BAD_SLACK                      # stuck at its start on the plateau
    new = E._Sides(spec, POOL3, models, marg, graded=True)
    assert new.vector(c).min() > E.BAD_SLACK
    w_new, s_new, _ = E._optimise(new, 3, 1, 42)
    assert not np.allclose(w_new, c) and s_new > E.BAD_SLACK and w_new[0] < c[0]            # moved toward the certifiable region
    # judge unchanged: the plateau recipe has no margin -> uncertified (or fail), never pass, margin None
    sm = {x.key: x for x in E.evaluate_blend(spec, POOL3, c, models, marg)}["smoke_point_mm"]
    assert sm.applied_margin is None and sm.status not in E.OK_STATUSES


def test_graded_objective_default_off_in_design(tmp_path):
    cdir, mp = _stub_cfg(tmp_path)
    r = E.design("test:jet_stub", POOL3, _plateau_models(), n_starts=3, seed=42, margins_path=mp, configs_dir=cdir, alt_limits=True)
    assert r.optimiser["graded_objective"] is False
    r2 = E.design("test:jet_stub", POOL3, _plateau_models(), n_starts=3, seed=42, margins_path=mp, configs_dir=cdir, alt_limits=True,
                  graded_objective=True)
    assert r2.optimiser["graded_objective"] is True


def test_graded_objective_is_continuous_at_the_margin_edge(tmp_path):
    """Review 3: no jump where the margin turns None (the old 0.01 floor); the penalty is applied inside and outside."""
    cdir, mp = _stub_cfg(tmp_path)
    spec = E.load_spec("test:jet_stub", cdir, alt_limits=True)
    marg = E.load_margins(mp)
    models = _plateau_models()
    s = E._Sides(spec, POOL3, models, marg, graded=True)
    def at(x):
        w = np.array([x, (1 - x) / 2, (1 - x) / 2])
        return s.vector(w)
    a, b = at(0.15), at(0.1500001)
    assert np.max(np.abs(a - b)) < 1e-4                          # continuous across the edge
    assert at(0.20).min() < at(0.16).min()                       # and the penalty keeps growing outside


def test_smoke_penalty_fn_inner_buffer_class_mix():
    """Review 3: the class-mix distance is measured to the box shrunk 0.005 at non-zero edges (iso-paraffin lower edge 0.007)."""
    from fuelblend.design import library_route as LR
    pool = ["a", "b", "c"]
    classes = {"a": "cycloalkane", "b": "n_paraffin", "c": "iso_paraffin"}
    pen = LR._smoke_penalty_fn(lambda p, w, T: 22.0, classes, lambda p, w, T: 12.0)
    assert pen(pool, np.array([0.60, 0.39, 0.01]), None) > 0.0       # inside the box but inside the buffer
    assert pen(pool, np.array([0.60, 0.378, 0.022]), None) == 0.0    # well inside
    assert pen(pool, np.array([0.61, 0.39, 0.0]), None) > pen(pool, np.array([0.60, 0.39, 0.01]), None)   # outside > buffer zone


def test_alt_point_view_relabels_alt_met_smoke(tmp_path):
    """Review 3: a smoke point estimate below the primary 25 mm but meeting the alt branch (>= 18, naphthalenes <= 3) is reported as
    alt met (design and API), not as a failed limit; a smoke below 18 stays failed."""
    import pandas as pd
    cdir, mp = _stub_cfg(tmp_path)
    marg = E.load_margins(mp)
    w = np.array([0.1, 0.45, 0.45])
    models = _stub_models()
    res = E.DesignResult(weights=dict(zip(POOL3, map(float, w))), predictions={}, margins={}, passed=False, extrapolating=False,
                         binding_constraints=[], constraints=[], spec_name="test:jet_stub", pool=list(POOL3), optimiser={})
    def fc(dp, ap):
        return pd.DataFrame([dict(constraint="smoke_point_mm", design_pred=dp, forward_pred=ap, limit=">= 25", lo=25.0, hi=float("nan"),
                                  design_meets=False, api_meets=False)])
    out, keys = LR.apply_alt_point_view(fc(24.9, 24.8), res, models, marg, "test:jet_stub", configs_dir=cdir)
    assert keys == ["smoke_point_mm"] and bool(out.loc[0, "design_meets"]) and bool(out.loc[0, "api_meets"])
    out, keys = LR.apply_alt_point_view(fc(24.9, 17.0), res, models, marg, "test:jet_stub", configs_dir=cdir)
    assert keys == [] and not bool(out.loc[0, "api_meets"])
    out, keys = LR.apply_alt_point_view(fc(17.5, 24.8), res, models, marg, "test:jet_stub", configs_dir=cdir)
    assert keys == [] and not bool(out.loc[0, "design_meets"])


def test_graded_penalty_reaches_the_objective_when_its_own_slack_is_not_binding(tmp_path):
    """Review 3: the smoke row's own slack (alt branch) is far from the minimum, so a penalty subtracted from that row alone never
    changes the max-min optimum.  The penalty must lower the whole objective (optimiser only)."""
    cdir, mp = _stub_cfg(tmp_path)
    spec = E.load_spec("test:jet_stub", cdir, alt_limits=True)
    marg = E.load_margins(mp)
    def models(pen):
        m = _stub_models(margin_fn=lambda pool, w, T: 3.0)             # margin never None
        m["props"]["smoke_point"] = dataclasses.replace(m["props"]["smoke_point"], margin_stand_in=3.0, margin_penalty_fn=pen)
        return m
    w_no, _, _ = E._optimise(E._Sides(spec, POOL3, models(lambda pool, w, T: 0.0), marg, graded=True), 3, 5, 42)
    w_pen, _, _ = E._optimise(E._Sides(spec, POOL3, models(lambda pool, w, T: 20.0 * max(float(w[0]) - 0.10, 0.0)), marg, graded=True), 3, 5, 42)
    assert w_pen[0] < w_no[0] - 0.005, (w_pen, w_no)


def test_rejudge_smoke_alt_met_scope_or_unknown_reason_never_passes(tmp_path, monkeypatch):
    """Review 4 (fail-closed): a STATUS_SCOPE row, or ANY reason not on the benign list, must not be upgraded to PASS."""
    scope = "margin validated on ['jet'] fuels only; this is a diesel spec, so the pass is not a certification"
    for status, why in ((E.STATUS_SCOPE, scope), (None, scope), (None, "some future engine blocker text")):
        r = _rejudge_stub(tmp_path, monkeypatch, why, engine_status=status)
        assert r.status == E.STATUS_UNCERTIFIED and r.status not in E.OK_STATUSES, (status, why, r.reason)
    assert _rejudge_stub(tmp_path, monkeypatch, "").status == E.STATUS_PASS            # the clean case still passes


def test_point_estimate_view_alt_does_not_mask_unrelated_violation():
    """Review 4: an alt-met smoke row is not a violation, but an unrelated violation still counts (nothing is subtracted/clamped)."""
    import types
    mk = lambda key, p, lo=None, hi=None: types.SimpleNamespace(key=key, prediction=p, lo=lo, hi=hi)   # noqa: E731
    cons = [mk("smoke_point_mm", 24.9, lo=25.0), mk("density_15C_kg_m3", 760.0, lo=775.0, hi=840.0), mk("flash_point_C", 50.0, lo=38.0)]
    assert LR.point_estimate_view_alt(cons, ["smoke_point_mm"]) == (1, 0)          # only density violates
    assert LR.point_estimate_view_alt(cons[:1] + cons[2:], ["smoke_point_mm"]) == (0, 0)
    assert LR.point_estimate_view_alt(cons, []) == (2, 0)
