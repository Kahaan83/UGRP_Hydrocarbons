"""Tests for predict_blend.py (T2 wrapper). Run: pytest tests/test_predict_blend.py -v

Molecules chosen for test coverage:
  n-decane  SMILES=CCCCCCCCCC   inchikey=DIOQZVSQGTUSAI-UHFFFAOYSA-N  (has pure rows, Antoine, LFL)
  n-dodecane SMILES=CCCCCCCCCCCC inchikey=SNRUBQQJIBEYMU-UHFFFAOYSA-N  (has pure rows, Antoine, LFL)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

# n-decane and n-dodecane: both in compounds.csv with Antoine and LFL
DEC_SMILES = "CCCCCCCCCC"
DOD_SMILES = "CCCCCCCCCCCC"
DEC_IK = "DIOQZVSQGTUSAI-UHFFFAOYSA-N"
DOD_IK = "SNRUBQQJIBEYMU-UHFFFAOYSA-N"


def _predict(*args, **kwargs):
    """Import-time deferred so the module need not be importable before it exists."""
    from fuelblend.api.predict_blend import predict
    return predict(*args, **kwargs)


# ── Test 1: fractions auto-normalize ─────────────────────────────────────────

def test_auto_normalize():
    """Off-by-tiny fractions (e.g. 0.5001 + 0.4999) normalize without error."""
    composition = [(DEC_SMILES, 0.5001), (DOD_SMILES, 0.4999)]
    result = _predict(composition, T_K=298.15)
    # check all three properties present and numeric
    for prop in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"):
        assert prop in result, f"missing {prop}"
        v = result[prop]["value"]
        assert v is not None and isinstance(v, float), f"{prop}.value must be float, got {v!r}"
    # no global warnings about normalization failure
    assert isinstance(result["warnings"], list)


# ── Test 2: unknown compound → warning, no crash ─────────────────────────────

def test_unknown_compound_warns():
    """An unrecognised SMILES triggers a warning and still returns a result dict."""
    composition = [(DEC_SMILES, 0.5), ("CCCOOOO[Zn]CCCC", 0.5)]  # not a HC
    result = _predict(composition, T_K=298.15)
    # global warnings should contain something about unknown/non-hydrocarbon
    all_warnings = result.get("warnings", []) + sum(
        (result.get(p, {}).get("warnings", []) for p in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C")),
        [],
    )
    found = any(
        any(kw in str(w).lower() for kw in ("unknown", "non_hc", "non-hc", "nonhc", "not found", "no smiles"))
        for w in all_warnings
    )
    assert found, f"Expected an unknown-compound warning. Got: {all_warnings}"


# ── Test 3: T outside domain → warning on every property ─────────────────────

def test_T_out_of_domain_warns():
    """T=600 K is outside [230, 470 K] → each property should carry a warning."""
    composition = [(DEC_SMILES, 0.5), (DOD_SMILES, 0.5)]
    result = _predict(composition, T_K=600.0)
    for prop in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"):
        warnings = result.get(prop, {}).get("warnings", [])
        found = any("T_out" in str(w) or "temperature" in str(w).lower() or "230" in str(w) or "470" in str(w)
                    for w in warnings)
        assert found, f"{prop}: expected T-range warning, got {warnings}"


# ── Test 4: flash is on the anchored-Antoine route; the nominal window warning replaces antoine_out_of_range ─────────

def test_flash_antoine_out_of_range_warns():
    """v5: flash is ideal Le Chatelier on the anchored Antoine registry, independent of T_K.  There is no fallback model:
    the method is always the route method, and a flash outside the nominal 0-100 C Antoine window carries the warning
    'antoine_nominal_window_exceeded' (a flash inside it does not).  T_K outside 230-470 K still carries the T-domain warning."""
    composition = [(DEC_SMILES, 0.5), (DOD_SMILES, 0.5)]
    result = _predict(composition, T_K=298.15)
    flash = result.get("flash_point_C", {})
    assert flash["method"] == "ideal_le_chatelier_anchored" and isinstance(flash["value"], float)
    out_of_window = not (0.0 <= flash["value"] <= 100.0)
    assert any("antoine_nominal_window_exceeded" in str(w) for w in flash["warnings"]) == out_of_window
    assert not any("antoine_out_of_range" in str(w) for w in flash["warnings"])
    # a heavy blend flashes above 100 C: the warning must appear
    heavy = _predict([("CCCCCCCCCCCCCCCCCC", 1.0)], T_K=298.15)["flash_point_C"]     # n-octadecane
    assert heavy["value"] is not None and not (0.0 <= heavy["value"] <= 100.0)
    assert any("antoine_nominal_window_exceeded" in str(w) for w in heavy["warnings"])
    result2 = _predict(composition, T_K=800.0)
    assert result2["flash_point_C"]["value"] == flash["value"]                         # flash does not depend on T_K
    all_w = result2.get("warnings", []) + result2.get("flash_point_C", {}).get("warnings", [])
    assert any("T_out" in str(w) or "230" in str(w) or "470" in str(w) for w in all_w), (
        f"Expected T-domain warning at T=800K, got {all_w}")


# ── Test 5: determinism ───────────────────────────────────────────────────────

def test_determinism():
    """Two calls with identical inputs return identical outputs."""
    composition = [(DEC_SMILES, 0.6), (DOD_SMILES, 0.4)]
    r1 = _predict(composition, T_K=310.0)
    r2 = _predict(composition, T_K=310.0)
    for prop in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"):
        assert r1[prop]["value"] == r2[prop]["value"], \
            f"{prop}.value changed between calls: {r1[prop]['value']} vs {r2[prop]['value']}"
        assert r1[prop].get("std") == r2[prop].get("std"), f"{prop}.std not deterministic"


# ── v4 Task 18: no lazy retrain; intervals from uncertainty_v4.yaml applied_margin only ──────────────

import math

import yaml


def _pb():
    import fuelblend.api.predict_blend as pb
    return pb


def test_missing_bundle_raises_not_trains(monkeypatch, tmp_path):
    """A missing mixture-residual bundle must raise FileNotFoundError naming the build command, and must not train."""
    pb = _pb()

    def _boom(*a, **k):
        raise AssertionError("a training function was invoked at load time")

    from sklearn.ensemble import HistGradientBoostingRegressor
    monkeypatch.setattr(HistGradientBoostingRegressor, "fit", _boom)
    for name in ("_ensure_mix_models", "_train_mix_models"):
        assert not hasattr(pb, name), f"{name} (lazy training) must be removed"
    monkeypatch.setattr(pb, "_MIX_MODEL_DIR", tmp_path / "empty_mix")
    monkeypatch.setattr(pb, "_mix_models", None)
    with pytest.raises(FileNotFoundError) as ei:
        pb._load_mix_models()
    msg = str(ei.value)
    assert "train_mixture_residual.py" in msg and ".venv/Scripts/python.exe" in msg
    assert not (tmp_path / "empty_mix").exists(), "nothing may be written (no training / no mkdir)"
    # D3b: no residual correction ships; the kwarg only adds a warning and the route value is unchanged
    base = _predict([(DEC_SMILES, 0.5), (DOD_SMILES, 0.5)], T_K=293.15)
    flagged = _predict([(DEC_SMILES, 0.5), (DOD_SMILES, 0.5)], T_K=293.15, use_mixture_residual=True)
    for prop in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"):
        assert flagged[prop]["value"] == base[prop]["value"]
        assert "mixture_residual_not_shipped_D3b" in flagged[prop]["warnings"]
    # pure-component bundle: same contract
    monkeypatch.setattr(pb, "_PURE_MODEL_DIR", tmp_path / "empty_pure")
    monkeypatch.setattr(pb, "_pure_models", None)
    with pytest.raises(FileNotFoundError, match="train_pure_models.py"):
        pb._load_pure_models()


def _fake_margins(flash_lab=4.0, flash_d93=0.0, visc=0.3, dens="inf"):
    def ent(v, unit, direction):
        e = {"value": v, "unit": unit, "direction": direction}
        if v != "inf":
            e["applied_margin"] = v
        return e
    return {
        "flash_point_C.lab_mixtures.ideal": ent(flash_lab, "degC", "over"),
        "flash_point_C.real_fuels_D93_reconstructed.ideal": ent(flash_d93, "degC", "over"),
        "viscosity_kin_cSt.crc647.recalibrated_refutas.all_T": ent(visc, "ln-ratio", "two_sided"),
        "density_kg_m3.crc647_ideal_volume.all_T": ent(dens, "kg/m3", "two_sided"),
    }


def test_interval_from_v4_yaml(monkeypatch):
    """interval = value +/- applied_margin; flash uses max(lab, D93); 'inf' margin -> (nan, nan), not in domain."""
    pb = _pb()
    comp = [(DEC_SMILES, 0.5), (DOD_SMILES, 0.5)]
    monkeypatch.setattr(pb, "_margins_cache", _fake_margins(flash_lab=4.0, flash_d93=0.0, visc=0.3, dens="inf"))
    r = _predict(comp, T_K=293.15, fuel_type="jet", composition_kind="reconstructed")   # 20 C: inside the viscosity window

    f = r["flash_point_C"]
    assert f["method"] == "ideal_le_chatelier_anchored"
    assert f["margin"] == 4.0                                 # max(4.0, 0.0): the 0 D93 floor is never used alone
    assert f["lo"] == pytest.approx(f["value"] - 4.0) and f["hi"] == pytest.approx(f["value"] + 4.0)
    assert f["in_domain"] is True and f["safe_bound"] == f["lo"]
    assert f["tuple"] == (f["value"], f["lo"], f["hi"], "ideal_le_chatelier_anchored", True)

    v = r["viscosity_dyn_mPas"]
    assert v["method"] == "refutas_recalibrated"
    assert v["lo"] == pytest.approx(v["value"] * math.exp(-0.3)) and v["hi"] == pytest.approx(v["value"] * math.exp(0.3))
    assert v["in_domain"] is True

    d = r["density_kg_m3"]                                    # 'inf' margin: no interval, never a number
    assert math.isnan(d["lo"]) and math.isnan(d["hi"]) and d["in_domain"] is False and d["margin"] is None

    # the real yaml: every applied_margin read is a finite number or the entry is unusable
    real = yaml.safe_load((ROOT / "configs" / "uncertainty_v4.yaml").read_text(encoding="utf-8"))["margins"]
    m, _, _, why = pb._applied_margin(["flash_point_C.lab_mixtures.ideal",
                                       "flash_point_C.real_fuels_D93_reconstructed.ideal"], real)
    assert m == pytest.approx(max(real["flash_point_C.lab_mixtures.ideal"]["applied_margin"], 0.0)) and why == ""
    m2, _, _, why2 = pb._applied_margin(["flash_point_C.real_fuels_D56.ideal"], real)   # value: inf
    assert m2 is None and "no valid margin" in why2


def test_not_certified_cases(monkeypatch):
    """Non-jet / undeclared fuel -> in_domain False; no-margin properties -> (nan, nan, in_domain False)."""
    pb = _pb()
    monkeypatch.setattr(pb, "_margins_cache", _fake_margins())
    comp = [(DEC_SMILES, 0.5), (DOD_SMILES, 0.5)]
    for ft in ("diesel", "gasoline", None):
        r = _predict(comp, T_K=293.15, fuel_type=ft)
        assert r["flash_point_C"]["in_domain"] is False
        assert r["flash_point_C"]["domain_reasons"]
    t = pb.predict_tuples(comp, T_K=293.15, fuel_type="jet")
    for prop in ("freezing_point_C", "smoke_point_mm", "ron", "mon", "dcn", "dvpe_kpa", "pour_point_C"):
        val, lo, hi, name, ok = t[prop]
        assert val is None and math.isnan(lo) and math.isnan(hi) and ok is False
    # temperature outside the viscosity window: raw Refutas, no certified interval
    r = _predict(comp, T_K=310.0, fuel_type="jet", composition_kind="reconstructed")
    assert r["viscosity_dyn_mPas"]["method"] != "refutas_recalibrated"
    assert r["viscosity_dyn_mPas"]["in_domain"] is False


def test_recal_beta_matches_models_v4_library():
    """The light recalibration fit equals the engine's Library fit (single source of truth for the correction)."""
    pb = _pb()
    from fuelblend.design.models_v4 import library
    b = library().recal_beta
    assert pb._recal_beta() == pytest.approx((float(b[0]), float(b[1])))


def test_pure_component_viscosity_follows_composition_kind():
    """known: raw Refutas for any n_comp (no single-component special case).  reconstructed: recalibrated path, uniformly."""
    pb = _pb()
    raw = pb._refutas([DEC_IK], [1.0], 293.15)[0]
    r = pb._predict_viscosity([DEC_IK], [1.0], [1.0], 293.15)               # default composition_kind="known"
    assert r["method"] == "refutas_raw" and r["value"] == pytest.approx(raw)
    mix = pb._predict_viscosity([DEC_IK, DOD_IK], [0.5, 0.5], [0.5, 0.5], 293.15)
    assert mix["method"] == "refutas_raw"
    rec = pb._predict_viscosity([DEC_IK, DOD_IK], [0.5, 0.5], [0.5, 0.5], 293.15, composition_kind="reconstructed")
    assert rec["method"] == "refutas_recalibrated"
