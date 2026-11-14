"""v4 review fix 4: predict(composition_kind='known') uses raw Refutas and lab margins; 'reconstructed' keeps the CRC path."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.api import predict_blend as pb  # noqa: E402
from fuelblend.api.predict_blend import predict  # noqa: E402

DEC, DOD = "CCCCCCCCCC", "CCCCCCCCCCCC"
ML = "refutas+ml_pair_gn"   # v5 (cc5f5ef, user decision): designated viscosity in the ML domain; Refutas kept as variants['refutas_raw']


def _refutas(r: dict) -> dict:
    """The Refutas result: the designated entry itself, or variants['refutas_raw'] where the ML value is designated."""
    return r["variants"]["refutas_raw"] if r["method"] == ML else r


def test_binary_is_continuous_with_pure():
    # 283.15 K: the pure-component rho/eta fits of decane AND dodecane both cover it, so the raw Refutas path is exercised.
    # (At 253.15 K the density fits start at 273.15 K / dodecane melts at 263.6 K, so the brief's decane/dodecane pair
    # silently drops to the gn_mass / gn_fallback paths; octane is used for that colder check.)
    pure = predict([(DEC, 1.0)], T_K=283.15, composition_kind="known")["viscosity_dyn_mPas"]
    near = predict([(DEC, 0.999), (DOD, 0.001)], T_K=283.15, composition_kind="known")["viscosity_dyn_mPas"]
    assert _refutas(pure)["method"] == _refutas(near)["method"] == "refutas_raw"
    assert pure["method"] == near["method"]
    assert abs(near["value"] / pure["value"] - 1) < 0.01
    assert abs(_refutas(near)["value"] / _refutas(pure)["value"] - 1) < 0.01
    # Cold check (253.15 K): v5 reads the frozen 230-470 K grid, so the cold path is raw Refutas too and stays continuous
    # (the v4 gn_mass / gn_fallback paths, a 12 % jump, are gone).  The cold Refutas lab bins have no margin: Refutas is
    # uncertified there; the designated ML value carries its own R-MLV-6 margin (T < 263 K bin).
    pc = predict([(DEC, 1.0)], T_K=253.15, composition_kind="known")["viscosity_dyn_mPas"]
    nc = predict([(DEC, 0.999), ("CCCCCCCC", 0.001)], T_K=253.15, composition_kind="known")["viscosity_dyn_mPas"]
    assert _refutas(pc)["method"] == _refutas(nc)["method"] == "refutas_raw"
    assert _refutas(pc)["in_domain"] is False and _refutas(nc)["in_domain"] is False
    assert abs(nc["value"] / pc["value"] - 1) < 0.01
    assert abs(_refutas(nc)["value"] / _refutas(pc)["value"] - 1) < 0.01


def test_known_composition_never_uses_recalibration():
    r = predict([(DEC, 0.5), (DOD, 0.5)], T_K=283.15, composition_kind="known")
    assert r["viscosity_dyn_mPas"]["method"] in ("refutas_raw", ML)
    assert _refutas(r["viscosity_dyn_mPas"])["method"] == "refutas_raw"
    r1 = predict([(DEC, 1.0)], T_K=293.15, composition_kind="known")
    assert _refutas(r1["viscosity_dyn_mPas"])["method"] == "refutas_raw"   # no single-component special case


def test_known_is_default():
    r = predict([(DEC, 0.5), (DOD, 0.5)], T_K=293.15)
    assert _refutas(r["viscosity_dyn_mPas"])["method"] == "refutas_raw"


def test_reconstructed_keeps_recalibrated_path_and_jet_requirement():
    r = predict([(DEC, 0.5), (DOD, 0.5)], T_K=293.15, fuel_type="jet", composition_kind="reconstructed")
    assert r["viscosity_dyn_mPas"]["method"] == "refutas_recalibrated" and r["viscosity_dyn_mPas"]["in_domain"] is True
    r2 = predict([(DEC, 0.5), (DOD, 0.5)], T_K=293.15, composition_kind="reconstructed")
    assert r2["viscosity_dyn_mPas"]["in_domain"] is False             # reconstructed still needs fuel_type='jet'


def test_bad_composition_kind_rejected():
    with pytest.raises(ValueError):
        predict([(DEC, 1.0)], composition_kind="guess")


def _ent(v, unit, d):
    e = {"value": v, "unit": unit, "direction": d}
    if v != "inf":
        e["applied_margin"] = v
    return e


def test_known_uses_lab_margins_by_temperature_bin(monkeypatch):
    margins = {
        "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K": _ent("inf", "ln-ratio", "two_sided"),
        "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_263_283K": _ent("inf", "ln-ratio", "two_sided"),
        "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_ge_283K": _ent(0.2, "ln-ratio", "two_sided"),
        "density_kg_m3.lab_mixtures.ideal_volume.over": _ent(4.0, "kg/m3", "over"),
        "density_kg_m3.lab_mixtures.ideal_volume.under": _ent(2.0, "kg/m3", "under"),
        "flash_point_C.lab_mixtures.ideal": _ent(4.0, "degC", "over"),
        "flash_point_C.real_fuels_D93_reconstructed.ideal": _ent(0.0, "degC", "over"),
    }
    monkeypatch.setattr(pb, "_margins_cache", margins)
    comp = [(DEC, 0.5), (DOD, 0.5)]
    warm = predict(comp, T_K=303.15, fuel_type=None, composition_kind="known")
    v = _refutas(warm["viscosity_dyn_mPas"])                       # Refutas keeps the yaml lab margin (monkeypatched 0.2)
    assert v["margin"] == pytest.approx(0.2) and v["lo"] == pytest.approx(v["value"] * math.exp(-0.2))
    assert v["in_domain"] is True                                   # lab margins are not jet specific
    d = warm["density_kg_m3"]
    assert d["lo"] == pytest.approx(d["value"] - 4.0) and d["hi"] == pytest.approx(d["value"] + 2.0)
    assert d["in_domain"] is True
    for T in (253.15, 273.15):                                      # both cold bins are inf: nothing borrowed
        c = _refutas(predict(comp, T_K=T, composition_kind="known")["viscosity_dyn_mPas"])
        assert c["margin"] is None and math.isnan(c["lo"]) and c["in_domain"] is False


def test_known_non_hydrocarbon_has_no_lab_margin(monkeypatch):
    monkeypatch.setattr(pb, "_margins_cache", {
        "density_kg_m3.lab_mixtures.ideal_volume.over": _ent(4.0, "kg/m3", "over"),
        "density_kg_m3.lab_mixtures.ideal_volume.under": _ent(2.0, "kg/m3", "under")})
    r = predict([("CCO", 0.5), (DEC, 0.5)], T_K=293.15, composition_kind="known")
    d = r["density_kg_m3"]
    assert d["in_domain"] is False and d["margin"] is None and math.isnan(d["lo"])
