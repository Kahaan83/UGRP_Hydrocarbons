"""v5 decision 1 (Task 6): the designated freezing route is UNIFAC Schroeder-van Laar in the API and the design judge."""
from pathlib import Path

import numpy as np
import pytest
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi

from fuelblend.api import predict_route
from fuelblend.api.predict_route import predict_fuel
from fuelblend.baselines import flash_unifac as fu
from fuelblend.design import engine_v4 as E
from fuelblend.design import library_route
from fuelblend.eval import composition_route as CR
from fuelblend.eval import freeze_multicomp as FM

ROOT = Path(__file__).resolve().parents[1]
UNIFAC_KEY = "freezing_point_C.lab_liquidus.svl_unifac"
IDEAL_KEY = "freezing_point_C.lab_liquidus.svl_verified"
_M = yaml.safe_load((ROOT / "configs" / "uncertainty_v5.yaml").read_text(encoding="utf-8"))["margins"]
MARGIN_UNIFAC = float(_M[UNIFAC_KEY]["applied_margin"])
MARGIN_IDEAL = float(_M[IDEAL_KEY]["applied_margin"])

DECANE, DODECANE, TDEC = "CCCCCCCCCC", "CCCCCCCCCCCC", "C1CC[C@H]2CCCC[C@@H]2C1"
BCH, TMB = "CCCCC1CCCCC1", "Cc1ccc(C)c(C)c1"
BLEND = {DECANE: 0.5, DODECANE: 0.3, TDEC: 0.2}          # n-decane / n-dodecane / trans-decalin, all with verified constants


def _ik(s):
    return inchi.MolToInchiKey(Chem.MolFromSmiles(s))


def _x(smiles, w):
    w = np.asarray(w, float)
    mw = np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smiles])
    x = (w / mw)
    return x / x.sum()


@pytest.fixture
def clear_caches():
    predict_route._tmap.cache_clear()
    predict_route._tmap_freezing.cache_clear()
    yield
    predict_route._tmap.cache_clear()
    predict_route._tmap_freezing.cache_clear()


def test_designated_model_is_unifac():
    assert CR.DESIGNATED_FREEZING_MODEL == "svl_unifac"
    assert CR.FREEZING_MODELS == {"ideal_svl": IDEAL_KEY, "svl_unifac": UNIFAC_KEY}


def test_margin_key_single_source():
    assert predict_route.freezing_margin_key() == library_route.freezing_margin_key() == UNIFAC_KEY
    assert not hasattr(predict_route, "FREEZING_MARGIN_KEY") and not hasattr(library_route, "FREEZING_MARGIN_KEY")


def test_api_freezing_unifac_with_ideal_variant():
    tm = CR.build_tmap_verified()
    assert all(_ik(s) in tm for s in BLEND)
    fp = predict_fuel(BLEND)["freezing_point_C"]
    assert fp["method"] == "svl_unifac_mole_fraction" and fp["uncertified"] is False
    assert abs(fp["safe_bound"] - fp["value"] - MARGIN_UNIFAC) < 1e-9
    assert fp["margin"] == pytest.approx(MARGIN_UNIFAC) and fp["margin_direction"] == "under"
    iv = fp["variants"]["ideal_svl"]
    assert iv["value"] is not None and abs(iv["safe_bound"] - iv["value"] - MARGIN_IDEAL) < 1e-9
    assert "coutinho_constants" in fp["variants"]
    smi = list(BLEND)
    x = _x(smi, list(BLEND.values()))
    iks = [_ik(s) for s in smi]
    cc = [(tm[k]["tm_K"], tm[k]["dh_fus_J_mol"]) for k in iks]
    ref = FM.svl_activity_C(x, cc, fu.make_unifac_gamma_fn(smi)[0])      # the delta_freeze_adopt model with thermo's own gamma
    assert fp["value"] == pytest.approx(ref, abs=1e-6)
    assert iv["value"] == pytest.approx(CR.ideal_sle_C(x, iks, tm), abs=1e-6)


def test_unifac_fallback_is_uncertified(monkeypatch):                    # Review Focus 1
    real = fu.unifac_subgroups
    can = Chem.MolToSmiles(Chem.MolFromSmiles(TDEC))
    monkeypatch.setattr(fu, "unifac_subgroups",
                        lambda s: None if Chem.MolToSmiles(Chem.MolFromSmiles(s)) == can else real(s))
    fp = predict_fuel(BLEND)["freezing_point_C"]
    assert fp["uncertified"] is True and "unifac_fallback" in " ".join(fp["warnings"])
    assert np.isnan(fp["hi"]) and "safe_bound" not in fp and fp["margin"] is None   # never scored under the UNIFAC margin
    assert fp["value"] is not None and fp["in_domain"] is False


def test_unifac_fallback_design_uncertified(monkeypatch):                 # Review Focus 1, design judge
    real = fu.unifac_subgroups
    can = Chem.MolToSmiles(Chem.MolFromSmiles(BCH))
    monkeypatch.setattr(fu, "unifac_subgroups",
                        lambda s: None if Chem.MolToSmiles(Chem.MolFromSmiles(s)) == can else real(s))
    models = library_route.build_route_models(["n-undecane", "n-dodecane", "butylcyclohexane"])
    fz = models["props"]["freezing"]
    assert "UNIFAC=fallback" in fz.margin_missing_reason
    margins = E.load_margins(ROOT / "configs" / "uncertainty_v5.yaml")
    c = E.evaluate_constraint("freezing_point_C", {"max": 0.0}, models["pool"], np.array([0.4, 0.2, 0.4]), models, margins, kind="jet")
    assert c.status == "uncertified"


def test_zero_weight_member_matches_removed():                            # Review Focus 5
    tm = CR.build_tmap_verified()
    smi = [DECANE, DODECANE, TDEC]
    iks = [_ik(s) for s in smi]
    x2 = _x(smi[:2], [0.6, 0.4])
    v1, i1 = CR.freezing_C(x=list(x2) + [0.0], iks=iks, smiles=smi, tmap=tm)
    v2, i2 = CR.freezing_C(x=list(x2), iks=iks[:2], smiles=smi[:2], tmap=tm)
    assert np.isfinite(v1) and abs(v1 - v2) < 1e-6
    assert i1 == {"model": "svl_unifac", "unifac_fallback": []} == i2


def test_zero_weight_member_design_matches_removed():                     # Review Focus 5, design freezing function
    m3 = library_route.build_route_models(["n-decane", "n-dodecane", "butylcyclohexane"])
    m2 = library_route.build_route_models(["n-decane", "n-dodecane"])
    v3 = m3["props"]["freezing"].fn(m3["pool"], np.array([0.6, 0.4, 0.0]), None)
    v2 = m2["props"]["freezing"].fn(m2["pool"], np.array([0.6, 0.4]), None)
    assert abs(v3 - v2) < 1e-6


def test_design_freezing_equals_api_freezing_unifac():
    """UNIFAC twin of test_v5_design_forward::test_design_freezing_equals_api_freezing: the design judge and predict_fuel give the same
    value; the recipe has three zero-weight pool members (pool-order gamma with x = 0 vs the API's positive-only composition)."""
    pool = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
    models = library_route.build_route_models(pool)
    w = {"n-undecane": 0.25, "n-dodecane": 0.15, "butylcyclohexane": 0.6}
    res = library_route.evaluate_recipe("astm:jet_a", models, w, margins_path=ROOT / "configs" / "uncertainty_v5.yaml")
    c = next(c for c in res.constraints if c.prop == "freezing")
    smi = dict(zip(models["pool"], [models["library"].smiles[n] for n in models["pool"]]))
    comp = {smi[n]: res.weights[n] for n in models["pool"] if res.weights[n] > 0}
    f = predict_fuel(comp, T_K=288.15, fuel_type="jet")["freezing_point_C"]
    assert f["method"] == "svl_unifac_mole_fraction" and abs(c.prediction - f["value"]) < 1e-9
    fz = models["props"]["freezing"]
    assert fz.margin_keys == (UNIFAC_KEY,) and fz.margin_missing_reason == ""
    assert c.applied_margin == pytest.approx(MARGIN_UNIFAC)


def test_vec_gamma_matches_thermo():
    rng = np.random.default_rng(7)
    pool = [DECANE, DODECANE, TDEC, BCH, TMB, "Cc1ccccc1", "CC(C)CC(C)(C)C", "c1ccc2ccccc2c1"]
    for _ in range(15):
        k = int(rng.integers(2, len(pool) + 1))
        smi = list(rng.choice(pool, k, replace=False))
        x = rng.random(k)
        x[rng.random(k) < 0.25] = 0.0
        if (x > 0).sum() == 0:
            continue
        x = x / x.sum()
        g_ref, fb_ref = fu.make_unifac_gamma_fn(smi)
        g_vec, fb_vec = fu.make_unifac_gamma_fn_vec(smi)
        assert fb_ref == fb_vec
        for T in (150.0, 230.0, 300.0, 410.0):
            np.testing.assert_allclose(g_vec(T, x), g_ref(T, x), rtol=1e-12, atol=0)


def test_freezing_C_matches_reference_model():
    tm = CR.build_tmap_verified()
    for smi, w in (([DECANE, DODECANE, TDEC], [0.5, 0.3, 0.2]), ([DODECANE, BCH, TMB], [0.3, 0.4, 0.3]),
                   ([DECANE, TMB], [0.9, 0.1])):
        x = _x(smi, w)
        iks = [_ik(s) for s in smi]
        cc = [(tm[k]["tm_K"], tm[k]["dh_fus_J_mol"]) for k in iks]
        ref = FM.svl_activity_C(x, cc, fu.make_unifac_gamma_fn(smi)[0])
        v, _ = CR.freezing_C(x, iks, smi, tm)
        assert v == pytest.approx(ref, abs=1e-8)
        vi, info = CR.freezing_C(x, iks, smi, tm, model="ideal_svl")
        assert vi == CR.ideal_sle_C(x, iks, tm) and info["model"] == "ideal_svl"


def test_switch_after_cache_clear(monkeypatch, clear_caches):             # Review Focus 3
    monkeypatch.setattr(CR, "DESIGNATED_FREEZING_MODEL", "ideal_svl")
    predict_route._tmap_freezing.cache_clear()
    fp = predict_fuel(BLEND)["freezing_point_C"]
    assert fp["method"] == "ideal_svl_mole_fraction" and "ideal_svl" not in fp["variants"]
    assert abs(fp["safe_bound"] - fp["value"] - MARGIN_IDEAL) < 1e-9
    assert predict_route.freezing_margin_key() == library_route.freezing_margin_key() == IDEAL_KEY
    monkeypatch.setattr(CR, "DESIGNATED_FREEZING_MODEL", "svl_unifac")
    predict_route._tmap_freezing.cache_clear()
    assert predict_fuel(BLEND)["freezing_point_C"]["method"] == "svl_unifac_mole_fraction"
