"""ML pair-residual viscosity in predict_fuel: designated in domain (user decision 2026-11-13), Refutas kept as a variant."""
import math

import pytest

from fuelblend.api import predict_route as pr

DEC, DOD, TOL = "CCCCCCCCCC", "CCCCCCCCCCCC", "Cc1ccccc1"
BLEND_A = {DEC: 0.5, TOL: 0.5}
BLEND_B = {DEC: 0.4, DOD: 0.3, TOL: 0.3}


def _fresh(blend, T, kind="known"):
    return pr.predict_fuel(blend, T_K=T, composition_kind=kind)


@pytest.mark.parametrize("blend", [BLEND_A, BLEND_B])
def test_variant_present_and_designated_unchanged(blend):
    r = _fresh(blend, 313.15)
    v = r["viscosity_kin_cSt"]
    ml = v["variants"]["ml_pair_residual"]
    assert ml["status"] == "ok" and ml["in_domain"] is True and ml["margin_unit"] == "ln"
    assert ml["margin"] == pytest.approx(0.13426357, abs=1e-6)
    # designated is the ML value; the Refutas result is kept as variants['refutas_raw']
    from fuelblend.eval.composition_route import refutas_nu_safe
    ctx = pr.get_ctx()
    comps = pr.resolve(blend)
    pm = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    nu, _ = refutas_nu_safe(ctx.prov, pm.idx, pm.w, 313.15)
    rr = v["variants"]["refutas_raw"]
    assert rr["method"] == "refutas_raw" and rr["value"] == pytest.approx(float(nu), rel=1e-12)
    assert ml["value"] == pytest.approx(rr["value"] * math.exp(ml["ln_correction"]))
    assert v["method"] == "refutas+ml_pair_gn" and v["value"] == ml["value"] and v["in_domain"] is True
    assert v["lo"] == pytest.approx(v["value"] * math.exp(-ml["margin"])) and v["hi"] == pytest.approx(v["value"] * math.exp(ml["margin"]))
    d = r["viscosity_dyn_mPas"]
    assert d["method"] == "refutas+ml_pair_gn"
    assert d["value"] == pytest.approx(d["variants"]["refutas_raw"]["value"] * math.exp(ml["ln_correction"]))


def test_reconstructed_keeps_refutas_designated():
    v = _fresh(BLEND_A, 313.15, "reconstructed")["viscosity_kin_cSt"]
    assert v["method"] != "refutas+ml_pair_gn"
    assert any("ml_not_designated" in w for w in v["warnings"])


def test_cold_margin_bins():
    for T, m in ((253.15, 0.28747778), (273.15, 0.60843819)):
        ml = _fresh(BLEND_A, T)["viscosity_kin_cSt"]["variants"]["ml_pair_residual"]
        if ml["status"] == "ok":
            assert ml["margin"] == pytest.approx(m, abs=1e-6)


def test_out_of_range_returns_none_with_warning():
    odd = "CC(C)CCCC(C)c1ccccc1"          # no own grid data
    ml = _fresh({DEC: 0.5, odd: 0.5}, 313.15)["viscosity_kin_cSt"]["variants"].get("ml_pair_residual")
    if ml is not None:
        assert ml["value"] is None and ml["margin"] is None and ml["in_domain"] is False
        assert any("ml_pair_residual_out_of_range" in w for w in ml["warnings"])
    ml2 = _fresh(BLEND_A, 330.0 + 100)["viscosity_kin_cSt"]["variants"].get("ml_pair_residual")  # beyond the training range
    if ml2 is not None:
        assert ml2["margin"] is None


def test_reconstructed_has_no_margin():
    ml = _fresh(BLEND_A, 313.15, "reconstructed")["viscosity_kin_cSt"]["variants"]["ml_pair_residual"]
    if ml["value"] is not None:
        assert ml["in_domain"] is False and ml["margin"] is None
        assert any("reconstructed" in w for w in ml["warnings"])


def test_no_ysi_ml_variant():
    r = _fresh(BLEND_A, 313.15)
    assert "ysi_route_ml_filled" not in r["smoke_point_mm"]["variants"]


def test_density_ml_variant_binary_has_margin_designated_unchanged():
    d = _fresh(BLEND_A, 298.15)["density_kg_m3"]
    assert d["method"] == "dens_idvol"
    v = d["variants"]["ml_pair_excess_volume"]
    if v["status"] == "ok":
        assert v["in_domain"] is True and v["margin"] == pytest.approx(4.4248287, abs=1e-6)
        assert abs(v["value"] - d["value"]) < 20.0


def test_density_ml_ternary_no_margin():
    v = _fresh(BLEND_B, 298.15)["density_kg_m3"]["variants"]["ml_pair_excess_volume"]
    assert v["margin"] is None and v["in_domain"] is False
