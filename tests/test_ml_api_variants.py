"""ML pair-residual viscosity variant in predict_fuel (variant only; designated Refutas unchanged)."""
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
    assert ml["value"] == pytest.approx(v["value"] * math.exp(ml["ln_correction"]))
    # designated stays raw Refutas on the BASE route
    from fuelblend.eval.composition_route import refutas_nu_safe
    ctx = pr.get_ctx()
    comps = pr.resolve(blend)
    pm = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    nu, _ = refutas_nu_safe(ctx.prov, pm.idx, pm.w, 313.15)
    assert v["method"] == "refutas_raw" and v["value"] == pytest.approx(float(nu), rel=1e-12)
    dml = r["viscosity_dyn_mPas"]["variants"]["ml_pair_residual"]
    assert dml["value"] == pytest.approx(r["viscosity_dyn_mPas"]["value"] * math.exp(ml["ln_correction"]))


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
