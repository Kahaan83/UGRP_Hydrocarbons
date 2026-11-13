"""v5 Task 6: RouteLibrary design on SMILES/InChIKey pools and the design-vs-forward check (R-API-4)."""
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import inchi

from fuelblend.design import engine_v4 as E
from fuelblend.design import library_route as LR
from fuelblend.design import models_v4 as M

ROOT = Path(__file__).resolve().parents[1]
MARGINS = ROOT / "configs" / "uncertainty_v5.yaml"
V5POOL = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
SPECS = {"jet_astm_d1655__v5pool": "astm:jet_a", "jet_bis_is1571__v5pool": "bis:jet_a1"}


def _smiles(names):
    d = pd.read_csv(LR.PURE_PROPS_V5).drop_duplicates("compound").set_index("compound").smiles
    return [d[n] for n in names]


@lru_cache(maxsize=None)
def _design_unswapped(pool_key, spec, n_starts=30):
    L = LR.RouteLibrary(list(pool_key))
    models = M.build_models(L)
    return E.design(spec, L.pool, models, n_starts=n_starts, seed=42, margins_path=MARGINS), L.pool


def test_name_pool_unchanged():
    """RouteLibrary on table names gives the library designs of D4 (build_models, not the freezing swap, not the flash overlay).
    Against the plain Library: 1e-9.  Against designs_v4.csv: 5.1e-5, because the CSV weights are rounded to 4 decimals."""
    d4 = pd.read_csv(ROOT / "results" / "v5" / "D4" / "designs_v4.csv").set_index("design_id")
    base = M.build_models(M.library(str(LR.PURE_PROPS_V5)))
    for did, spec in SPECS.items():
        res, pool = _design_unswapped(tuple(V5POOL), spec)
        ref = E.design(spec, V5POOL, base, n_starts=30, seed=42, margins_path=MARGINS)
        for n in V5POOL:
            assert abs(res.weights[n] - ref.weights[n]) < 1e-9
        csv = json.loads(d4.loc[did, "weights"])                 # the committed CSV is rounded to 4 decimals
        for n in V5POOL:
            assert abs(res.weights[n] - csv[n]) < 5.1e-5


def test_smiles_pool_equals_name_pool():
    spec = "astm:jet_a"
    by_name = LR.build_route_models(V5POOL)
    ik = [inchi.MolToInchiKey(Chem.MolFromSmiles(s)) for s in _smiles(V5POOL)]
    mixed = [_smiles(V5POOL)[0], ik[1], ik[2], V5POOL[3], ik[4], _smiles(V5POOL)[5]]
    by_id = LR.build_route_models(mixed)
    assert by_id["pool"] == V5POOL
    assert not by_id["library"].nonverified_used              # nothing recorded by the build probes (freezing records on use)
    a = E.design(spec, by_name["pool"], by_name, n_starts=8, seed=42, margins_path=MARGINS)
    b = E.design(spec, by_id["pool"], by_id, n_starts=8, seed=42, margins_path=MARGINS)
    for n in V5POOL:
        assert abs(a.weights[n] - b.weights[n]) < 1e-9
    assert {x.split("=")[0].split(".")[1] for x in by_id["library"].nonverified_used} <= {"Tm", "dHfus"}   # freezing only


def _first_absent_nalkane():
    L = LR.RouteLibrary(V5POOL)
    tab = set(L._table_ik)
    for n in range(5, 19):
        s = "C" * n
        ik = inchi.MolToInchiKey(Chem.MolFromSmiles(s))
        if ik in tab:
            continue
        try:
            LR.RouteLibrary(V5POOL + [s])
            models = LR.build_route_models(V5POOL + [s])
        except E.ModelUnavailable:
            continue
        return s, models
    raise AssertionError("no usable absent n-alkane")


def test_nontable_compound_uncertifies():
    s, models = _first_absent_nalkane()
    assert s == "C" * 9                                           # n-nonane: first n-alkane absent from the table
    res = E.design("astm:jet_a", models["pool"], models, n_starts=8, seed=42, margins_path=MARGINS)
    cons = {c.key: c for c in res.constraints}
    assert len(models["pool"]) == 7
    for key in ("density_15C_kg_m3", "flash_point_C", "freezing_point_C", "aromatics_vol_pct" if "aromatics_vol_pct" in cons else "density_15C_kg_m3"):
        assert cons[key].status == "uncertified", key
        assert "non-verified constant" in cons[key].reason, (key, cons[key].reason)
        assert s in cons[key].reason
    assert res.passed is False
    assert any(c.status != "uncertified" for c in res.constraints)  # the verified-only constraints keep their own status


def test_nontable_constant_statuses():
    s, models = _first_absent_nalkane()
    L = models["library"]
    assert L.status(s, "MW") == "verified" and L.status(s, "rho15") == "route_grid"
    assert L.status(s, "Tm") == "coutinho_table" and L.status(s, "dHfus") == "coutinho_table"
    assert L.status(s, "flash") == "route_anchored"
    for p in ("RON", "MON", "DCN", "Psat37_8"):
        assert L.status(s, p) == "absent"
        with pytest.raises(E.ModelUnavailable):
            L.v(s, p)


def test_member_without_own_data_refused():
    c40 = "C" * 40
    with pytest.raises(E.ModelUnavailable) as ei:
        LR.build_route_models(V5POOL + [c40])
    assert c40 in str(ei.value)


@pytest.mark.parametrize("did", list(SPECS))
def test_forward_check_d4(did):
    d4 = pd.read_csv(ROOT / "results" / "v5" / "D4" / "designs_v4.csv").set_index("design_id")
    w = json.loads(d4.loc[did, "weights"])
    models = LR.build_route_models(V5POOL)
    res = LR.evaluate_recipe(SPECS[did], models, w, margins_path=MARGINS)
    df = LR.design_forward_check(res)
    assert list(df.columns) == LR.CHECK_COLUMNS
    judged = df[df.tol.notna()]
    assert set(judged.constraint) >= {"flash_point_C[tol]", "nhoc_MJ_kg[compH_rdkit_MW]"} and len(judged) >= 3
    assert (judged["pass"] == True).all(), df.to_string()     # noqa: E712
    # every constrained property has a limit row judged on BOTH the design and the API value
    lim = df[df.limit != ""]
    assert set(lim.constraint) == {c.key for c in res.constraints}, df.to_string()
    for key in ("density_15C_kg_m3", "freezing_point_C", "smoke_point_mm"):
        r = df[df.constraint == key].iloc[0]
        assert isinstance(r.cause, str) and r.cause.strip()
        assert r.design_meets in (True, False) and r.api_meets in (True, False)


def test_category_b_designs_meet_limits_on_api_value():
    """Review C1: a design that is called 'meets limits' must meet every limit at the API (predict_fuel) point estimate too."""
    from fuelblend.design import library_route as lr
    models = LR.build_route_models(V5POOL)
    for spec in ("astm:jet_a", "bis:jet_a1"):
        res = LR.design_routed(spec, models, n_starts=8, seed=42, margins_path=MARGINS)
        fc = LR.design_forward_check(res)
        lim = fc[fc.limit != ""]
        both = lr.limits_met_by_both(fc)
        assert both == bool(lim.design_meets.all() and lim.api_meets.all())
        if both:
            assert (lim.forward_pred.map(np.isfinite)).all()


def test_design_freezing_equals_api_freezing():
    """The design freezing model is the API function (ideal SvL, mole fraction, designated verified-else-Coutinho tmap): same value for a recipe."""
    from fuelblend.api import predict_route as pr
    models = LR.build_route_models(V5POOL)
    w = {"n-undecane": 0.25, "n-dodecane": 0.15, "butylcyclohexane": 0.6}
    res = LR.evaluate_recipe("astm:jet_a", models, w, margins_path=MARGINS)
    c = next(c for c in res.constraints if c.prop == "freezing")
    comp = {s: res.weights[n] for n, s in zip(V5POOL, _smiles(V5POOL)) if res.weights[n] > 0}
    f = pr.predict_fuel(comp, T_K=288.15, fuel_type="jet")["freezing_point_C"]["value"]
    assert abs(c.prediction - f) < 1e-9
    fz = models["props"]["freezing"]                      # v5 delta_freeze_margin: verified constants carry the lab-liquidus margin key
    assert fz.margin_keys == ("freezing_point_C.lab_liquidus.svl_verified",) and fz.margin_missing_reason == ""


def test_forward_check_none_safe(monkeypatch):
    """A forward value of None gives api_meets False with a cause, not an exception."""
    from fuelblend.api import predict_route as pr
    models = LR.build_route_models(V5POOL)
    res = LR.evaluate_recipe("astm:jet_a", models, {"n-undecane": 0.3, "butylcyclohexane": 0.7}, margins_path=MARGINS)
    real = pr.predict_fuel

    def fake(*a, **k):
        out = real(*a, **k)
        out["viscosity_kin_cSt"] = dict(out["viscosity_kin_cSt"], value=None, status="T_out_of_grid", domain_reasons=["T"])
        out["smoke_point_mm"] = dict(out["smoke_point_mm"], value=None, status="no_evaluation")
        out["flash_point_C"] = dict(out["flash_point_C"], value=None, status="low_coverage", variants={})
        return out
    monkeypatch.setattr(pr, "predict_fuel", fake)
    df = LR.design_forward_check(res)
    for key in ("kin_visc_m20C_cSt", "smoke_point_mm", "flash_point_C"):
        r = df[df.constraint == key].iloc[0]
        assert r.api_meets is False or r.api_meets == False   # noqa: E712
        assert "None" in r.cause
    assert not LR.limits_met_by_both(df)


def test_flash_on_overlay_and_note():
    from fuelblend.api import predict_route as pr
    models = LR.build_route_models(V5POOL)
    res = LR.evaluate_recipe("astm:jet_a", models, {"n-decane": 0.3, "n-dodecane": 0.3, "butylcyclohexane": 0.4}, margins_path=MARGINS)
    comp = {s: res.weights[n] for n, s in zip(V5POOL, _smiles(V5POOL)) if res.weights[n] > 0}
    fl = pr.predict_fuel(comp, T_K=288.15, fuel_type="jet")["flash_point_C"]
    c = next(c for c in res.constraints if c.prop == "flash")
    assert abs(c.prediction - fl["value"]) < 0.01
    if c.applied_margin is not None and np.isfinite(c.applied_margin):
        assert LR.LAB_MARGIN_NOTE in c.reason


def test_default_n_starts():
    assert LR.default_n_starts(range(21)) == 8 and LR.default_n_starts(range(6)) == 30


def test_probe_failure_uncertifies(monkeypatch):
    models = LR.build_route_models(V5POOL)
    L = models["library"]
    pm = models["props"]["density"]

    def boom(pool, w, T_C):
        raise RuntimeError("boom")
    bad = LR._flag(L, __import__("dataclasses").replace(pm, fn=boom))
    assert "probe failed" in bad.margin_missing_reason and not bad.margin_keys_by_side
    models["props"]["density"] = bad
    res = LR.evaluate_recipe("astm:jet_a", models, {"n-decane": 0.5, "butylcyclohexane": 0.5}, margins_path=MARGINS)
    c = next(c for c in res.constraints if c.prop == "density")
    assert c.status == "uncertified" and "probe failed" in c.reason


def test_visc_T_mismatch_raises():
    models = LR.build_route_models(V5POOL)
    res = LR.evaluate_recipe("astm:jet_a", models, {"n-decane": 0.5, "butylcyclohexane": 0.5}, margins_path=MARGINS)
    with pytest.raises(ValueError):
        LR.design_forward_check(res, T_C_visc=-40.0)


def test_smoke_margin_design_side_domain_gated():
    """Smoke margin on the design side, domain- and class-mix-gated.  User decision R-SC-2b, 'Yes, replace (Recommended)': inside the
    regression domain AND the real-fuel class-mix range the margin is 4.107069860688682 mm (configs/uncertainty_v5.yaml, replacing the
    7.8 mm lab margin); in the domain but outside the class-mix range it is uncertified (reason class_mix_outside_real_fuel_range, margin
    None); out of the domain (0 vol % aromatics) it stays uncertified with the OutOfDomainWarning reason, as before."""
    import math
    from fuelblend.design import library_route as LR
    new_margin = 4.107069860688682
    models = LR.build_route_models(["n-decane", "n-dodecane", "isooctane", "butylcyclohexane", "1,2,4-trimethylbenzene"])
    mp = ROOT / "configs/uncertainty_v5.yaml"
    out = LR.evaluate_recipe("astm:jet_a", models, {"n-decane": 0.5, "n-dodecane": 0.5}, margins_path=mp)
    s = next(c for c in out.constraints if c.key == "smoke_point_mm")
    assert s.status == "uncertified" and s.applied_margin is None and "OutOfDomainWarning" in s.reason
    ind = LR.evaluate_recipe("astm:jet_a", models, {"n-decane": 0.2, "n-dodecane": 0.1, "isooctane": 0.2, "butylcyclohexane": 0.3,
                                                    "1,2,4-trimethylbenzene": 0.2}, margins_path=mp)
    t = next(c for c in ind.constraints if c.key == "smoke_point_mm")
    assert "OutOfDomainWarning" not in t.reason and t.applied_margin == new_margin
    assert t.slack == pytest.approx((t.prediction - new_margin - 25.0) / 25.0)
    assert t.status == ("fail" if t.slack < 0 else t.status)
    nomix = LR.evaluate_recipe("astm:jet_a", models, {"n-decane": 0.45, "n-dodecane": 0.4, "1,2,4-trimethylbenzene": 0.15}, margins_path=mp)
    u = next(c for c in nomix.constraints if c.key == "smoke_point_mm")        # in the aromatics domain, but no iso-paraffin
    assert "OutOfDomainWarning" not in u.reason and u.applied_margin is None and u.status != "pass"


def test_design_viscosity_is_api_designated_with_ml_margin():
    """Design viscosity = predict_fuel designated value (ML pair residual in its domain, cc5f5ef) and is judged with the
    R-MLV-6 ln margin of the T bin (0.287 below 263 K); out of the ML domain Refutas and its yaml margin stay."""
    import math
    from fuelblend.api import predict_route as pr
    models = LR.build_route_models(V5POOL)
    w = {"n-undecane": 0.25, "n-dodecane": 0.15, "butylcyclohexane": 0.6}
    res = LR.evaluate_recipe("astm:jet_a", models, w, margins_path=MARGINS)
    c = next(c for c in res.constraints if c.prop == "viscosity")
    comp = {s: res.weights[n] for n, s in zip(V5POOL, _smiles(V5POOL)) if res.weights[n] > 0}
    api = pr.predict_fuel(comp, T_K=253.15, fuel_type="jet")["viscosity_kin_cSt"]
    assert api["method"] == "refutas+ml_pair_gn"
    assert abs(c.prediction / api["value"] - 1) < 1e-6
    assert c.model == "refutas+ml_pair_gn" and c.applied_margin == pytest.approx(api["margin"])
    assert c.applied_margin == pytest.approx(0.28747778594508044)
    assert c.slack == pytest.approx((8.0 - c.prediction * math.exp(c.applied_margin)) / 8.0)
    assert c.status == "pass"
    st = LR.visc_ml_status(["CCCCCCCCCC"], [1.0], 150.0)     # outside the 233-393 K training range: Refutas kept
    assert st["in_domain"] is False and st["margin"] is None and st["g"] == 0.0
