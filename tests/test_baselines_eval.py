"""Unit tests for scripts/baselines_eval.py.  Run: python -m pytest tests/test_baselines_eval.py

ALL numbers in this file are INVENTED FIXTURES (hand-picked pure viscosities, densities and Antoine
constants). They exist to prove the code implements the formulas in AGENTS.md section 2 exactly; they are
NOT data and no result in reports/BASELINES.md comes from this file.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "scripts"))
import baselines_eval as B  # noqa: E402

# --------------------------------------------------------------------------- fixtures
T0 = 300.0
ETA = {"A": {300.0: 1.0, 350.0: 0.7}, "B": {300.0: 2.0, 350.0: 1.2}, "C": {300.0: 0.02, 350.0: 0.015}}
RHO = {"A": {300.0: 800.0, 350.0: 780.0}, "B": {300.0: 900.0, 350.0: 880.0}, "C": {300.0: 800.0, 350.0: 780.0}}
MW = {"A": 100.0, "B": 200.0, "C": 100.0}


def pure_rows():
    rows = []
    for k, tab in ETA.items():
        for T, v in tab.items():
            rows.append(dict(inchikey=k, property="viscosity_dyn_mPas", T_K=T, value=v))
    for k, tab in RHO.items():
        for T, v in tab.items():
            rows.append(dict(inchikey=k, property="density_kg_m3", T_K=T, value=v))
    return pd.DataFrame(rows)


def compound_registry():
    return pd.DataFrame([
        dict(inchikey="A", mw=100.0, antoine_A=6.0, antoine_B=1600.0, antoine_C=-50.0, lfl_vol_pct=2.0),
        dict(inchikey="B", mw=200.0, antoine_A=6.0, antoine_B=1600.0, antoine_C=-50.0, lfl_vol_pct=1.0),
        dict(inchikey="C", mw=100.0, antoine_A=np.nan, antoine_B=np.nan, antoine_C=np.nan, lfl_vol_pct=1.5),
        dict(inchikey="D", mw=100.0, antoine_A=6.0, antoine_B=1600.0, antoine_C=-50.0, lfl_vol_pct=np.nan),
    ]).set_index("inchikey")


# ------------------------------------------------------------------- pure curves
def test_pure_curve_reproduces_its_own_nodes():
    fits = B.build_fits(pure_rows())
    assert abs(B.eta_at(fits, "A", 300.0) - 1.0) < 1e-10
    assert abs(B.eta_at(fits, "B", 350.0) - 1.2) < 1e-10
    assert abs(B.rho_at(fits, "A", 300.0) - 800.0) < 1e-9
    # outside the measured range (+2 K pad) -> no value, no extrapolation
    assert B.eta_at(fits, "A", 400.0) is None
    assert B.rho_at(fits, "B", 250.0) is None


def test_log_linear_interpolation_is_monotone_between_nodes():
    fits = B.build_fits(pure_rows())
    mid = B.eta_at(fits, "A", 320.0)
    assert 0.7 < mid < 1.0  # invented fixture: viscosity falls with T


# --------------------------------------------------------------- mole conversion
def test_mole_from_mass_hand_numbers():
    x = B.mole_from_mass(["A", "B"], [0.5, 0.5], MW)
    assert abs(x[0] - 2.0 / 3.0) < 1e-12 and abs(x[1] - 1.0 / 3.0) < 1e-12
    x = B.mole_from_mass(["A", "B"], [0.3, 0.7], MW)
    assert abs(x[0] - 0.3 / 100 / (0.3 / 100 + 0.7 / 200)) < 1e-12
    assert B.mole_from_mass(["A", "Z"], [0.5, 0.5], MW) is None  # unknown molar mass


# ------------------------------------------------------------------ GN baselines
def test_gn_mass_recovers_exact_log_mass_formula():
    """Two invented binary mixtures: gn_mass must equal exp(sum w_i ln eta_i) exactly."""
    fits = B.build_fits(pure_rows())
    for w in ([0.3, 0.7], [0.6, 0.4]):
        got, why = B.gn_from_frac(["A", "B"], w, fits, T0)
        want = math.exp(w[0] * math.log(ETA["A"][T0]) + w[1] * math.log(ETA["B"][T0]))
        assert why == "" and abs(got - want) < 1e-12, (got, want)


def test_gn_mole_uses_mole_fractions():
    fits = B.build_fits(pure_rows())
    w = [0.3, 0.7]
    x = B.mole_from_mass(["A", "B"], w, MW)
    got, _ = B.gn_from_frac(["A", "B"], x, fits, T0)
    want = math.exp(x[0] * math.log(1.0) + x[1] * math.log(2.0))
    assert abs(got - want) < 1e-12


def test_gn_skips_records_without_component_curve():
    fits = B.build_fits(pure_rows())
    got, why = B.gn_from_frac(["A", "B"], [0.5, 0.5], fits, 400.0)
    assert got is None and why == "no_pure_eta_at_T"


# ------------------------------------------------------------------ Refutas + density
def test_refutas_hand_computed():
    fits = B.build_fits(pure_rows())
    w = [0.3, 0.7]
    eta_dyn, nu_mix, rho_mix, why = B.refutas(["A", "B"], w, fits, T0)
    nu_a = ETA["A"][T0] / (RHO["A"][T0] / 1000.0)
    nu_b = ETA["B"][T0] / (RHO["B"][T0] / 1000.0)
    vbn = w[0] * (14.534 * math.log(math.log(nu_a + 0.8)) + 10.975) + \
          w[1] * (14.534 * math.log(math.log(nu_b + 0.8)) + 10.975)
    want_nu = math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8
    want_rho = 1.0 / (w[0] / RHO["A"][T0] + w[1] / RHO["B"][T0])
    assert why == ""
    assert abs(nu_mix - want_nu) < 1e-10
    assert abs(rho_mix - want_rho) < 1e-9
    assert abs(eta_dyn - want_nu * want_rho / 1000.0) < 1e-10
    # the blend sits between the pure kinematic viscosities
    assert min(nu_a, nu_b) < nu_mix < max(nu_a, nu_b)


def test_refutas_rejects_nu_below_0p2_and_missing_density():
    fits = B.build_fits(pure_rows())
    _, _, _, why = B.refutas(["A", "C"], [0.5, 0.5], fits, T0)   # C: 0.02 mPa.s / 0.8 = 0.025 cSt
    assert why == "nu_lt_0p2"
    no_rho = B.build_fits(pure_rows()[pure_rows().property == "viscosity_dyn_mPas"])
    _, _, _, why = B.refutas(["A", "B"], [0.5, 0.5], no_rho, T0)
    assert why == "no_pure_rho_at_T"


def test_dens_idvol_is_harmonic_mass_mixing():
    fits = B.build_fits(pure_rows())
    got, why = B.dens_idvol(["A", "B"], [0.3, 0.7], fits, T0)
    want = 1.0 / (0.3 / 800.0 + 0.7 / 900.0)
    assert why == "" and abs(got - want) < 1e-9
    # pure record is the identity: v_i = 1
    got1, _ = B.dens_idvol(["A"], [1.0], fits, T0)
    assert abs(got1 - 800.0) < 1e-9


# ------------------------------------------------------------------ flash point
def _flash_expected(x, A, BB, C, lfl):
    """Independent inversion of H=1 for the special case of shared Antoine constants."""
    psat = B.P_TOTAL / sum(xi / l for xi, l in zip(x, lfl))
    return BB / (A - math.log10(psat)) - C - 273.15


def test_psat_kpa_convention():
    assert abs(B.psat_kpa(6.0, 1000.0, -50.0, 300.0) - 100.0) < 1e-9   # log10(P/kPa) = A - B/(T+C)


def test_chat_flash_binary_matches_hand_computation():
    """x = [0.5, 0.5], LFL 2 % and 1 %, same Antoine -> H=1 at a hand-computable T."""
    comp = compound_registry()
    x = [0.5, 0.5]
    pred, why = B.flash_prediction(["A", "B"], x, comp)
    want = _flash_expected(x, 6.0, 1600.0, -50.0, [0.02, 0.01])
    assert why == "" and abs(pred - want) < 1e-3, (pred, want)
    pred_lt, why_lt = B.flash_prediction(["A", "B"], x, comp, lfl_temp=True)
    assert why_lt == "" and pred_lt < pred        # Zabetakis LFL(T) lowers the flash point


def test_chat_flash_pure_identity():
    comp = compound_registry()
    pred, why = B.flash_prediction(["A"], [1.0], comp)
    want = _flash_expected([1.0], 6.0, 1600.0, -50.0, [0.02])
    assert why == "" and abs(pred - want) < 1e-3


def test_chat_flash_reports_missing_inputs():
    comp = compound_registry()
    assert B.flash_prediction(["C"], [1.0], comp)[1] == "no_antoine"
    assert B.flash_prediction(["D"], [1.0], comp)[1] == "unknown_LFL"
    assert B.flash_prediction(["A", "A"], [0.5, 0.5], comp)[1] == ""   # duplicate key is still solvable


# --------------------------------------------------------------- record pass + ids
def _model_table():
    def row(prop, sys_, keys, x, w, val, kind, T=T0, P=np.nan):
        return dict(property=prop, system_id=sys_, T_K=T, P_kPa=P, kind=kind, n_comp=len(keys),
                    value=val, conflict=False, tier="A", sources="fixture",
                    comp_inchikeys=json.dumps(keys), x_mole=json.dumps(x), w_mass=json.dumps(w),
                    split="S1")
    rows = [row("viscosity_dyn_mPas", "pureA", ["A"], [1.0], [1.0], 1.0, "pure"),
            row("viscosity_dyn_mPas", "pureB", ["B"], [1.0], [1.0], 2.0, "pure"),
            row("viscosity_dyn_mPas", "mixAB", ["A", "B"], [2 / 3, 1 / 3], [0.3, 0.7], 1.8, "binary"),
            row("flash_point_C", "mixABf", ["A", "B"], [0.5, 0.5], [0.5, 0.5], 40.0, "binary", T=np.nan),
            row("density_kg_m3", "fuel:9", [], [], [], 800.0, "real_fuel", T=np.nan)]
    return pd.DataFrame(rows)


def test_predict_records_wiring_and_unique_ids():
    mt = _model_table()
    comp = compound_registry().assign(mw=[100.0, 200.0, 100.0, 100.0])
    fits = B.build_fits(pure_rows())
    pred, skip = B.predict_records(mt, comp, fits)
    assert pred.baseline_row_id.is_unique and pred.eval_key.is_unique
    mix = pred[pred.system_id == "mixAB"].iloc[0]
    want = math.exp(0.3 * math.log(1.0) + 0.7 * math.log(2.0))
    assert abs(mix.gn_mass - want) < 1e-12
    x_ab = B.mole_from_mass(["A", "B"], [0.3, 0.7], MW)   # recomputed from w_mass, not the stored x_mole
    assert abs(mix.gn_mole - math.exp(x_ab[0] * math.log(1.0) + x_ab[1] * math.log(2.0))) < 1e-12
    fl = pred[pred.system_id == "mixABf"].iloc[0]
    x_abf = B.mole_from_mass(["A", "B"], [0.5, 0.5], MW)
    assert abs(fl.chat_nominal - _flash_expected(x_abf, 6.0, 1600.0, -50.0, [0.02, 0.01])) < 1e-3
    assert fl.chat_lflT < fl.chat_nominal
    real = pred[pred.kind == "real_fuel"].iloc[0]
    assert pd.isna(real.dens_idvol) and skip[("density_kg_m3", "dens_idvol", "no_molecule_list")] == 1


# ------------------------------------------------------------------- split exclusion
def test_apply_splits_drops_locked_system_ids():
    mt = pd.DataFrame({"property": ["density_kg_m3"] * 3, "system_id": ["a", "b", "c"],
                       "kind": ["binary"] * 3, "value": [1.0, 2.0, 3.0]})
    splits = {"S1": [0, 1], "S2": [2], "S3": []}
    locked = {"multi_system_ids_locked": ["b"], "real_fuel_ids_locked": []}
    out, info = B.apply_splits(mt, splits, locked)
    assert "b" not in set(out.system_id)                       # never evaluated
    assert set(out.system_id) == {"a", "c"}
    assert out.set_index("system_id").loc["c", "split"] == "S2"
    assert info["locked_rows_removed_from_splits"] == 1        # the manifest had leaked it

    locked2 = {"multi_system_ids_locked": ["zzz"], "real_fuel_ids_locked": []}
    out2, info2 = B.apply_splits(mt, splits, locked2)
    assert info2["locked_rows_removed_from_splits"] == 0 and len(out2) == 3


def test_split_labels_absent_from_real_manifest(mt_and_splits=None):
    """Against the real manifests: no locked id survives in a split (read-only check)."""
    mt, comp, splits, locked = B.load_inputs()
    out, info = B.apply_splits(mt, splits, locked)
    lk = set(locked["multi_system_ids_locked"])
    assert not (set(out.loc[out.split.isin(B.SPLIT_NAMES), "system_id"]) & lk)
    assert set(out.split.unique()) <= set(B.SPLIT_NAMES) | {""}
