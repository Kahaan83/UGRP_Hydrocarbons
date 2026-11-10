"""v4 Task 12: density / net heat of combustion / freezing point on the v4 corpus (scripts/v4/eval_bulk_props_v4.py)."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("eval_bulk_props_v4", ROOT / "scripts" / "v4" / "eval_bulk_props_v4.py")
ev = importlib.util.module_from_spec(_spec)
sys.modules["eval_bulk_props_v4"] = ev
_spec.loader.exec_module(ev)


@pytest.fixture(scope="module")
def run():
    fuels, role, comp, props, _ = ev.load_inputs()
    dens = ev.eval_density(comp, props, role, fuels, ev.PureDensity(comp))
    hoc = ev.eval_hoc(comp, props, role, fuels)
    fr = ev.eval_freezing(comp, props, role, fuels)
    return dict(fuels=fuels, role=role, comp=comp, dens=dens, hoc=hoc, fr=fr)


def test_v4_inputs_are_clean(run):
    for key in ("dens", "hoc", "fr"):
        df = run[key]
        assert len(df) > 0
        ev.assert_inputs_clean(df)                       # raw_file targets, no quarantined provenance, no dropped fuel
        assert not df.composition_provenance.str.contains("quarantin").any()
    inv = set(pd.read_csv(ev.INV).fuel_id)
    assert set(run["fuels"]) <= inv
    assert all(f.startswith("CRC647") for f in run["fuels"])
    assert not (set(run["fuels"]) & set(pd.read_csv(ev.DROPPED).fuel_id))


def _raw_value_exists(raw, fuel_id, prop, value, T_C=None, scale=1.0, col="value_std", tol=0.011):
    code = int(fuel_id.split("_")[1])
    r = raw[(raw.code == code) & (raw.property == prop)]
    if T_C is not None:
        r = r[r.T_C.isna() | ((r.T_C - T_C).abs() < 0.06)]   # D4052 rows carry no T in the raw file (60 F)
    return bool(((r[col] * scale - value).abs() < tol).any())


def test_evaluated_targets_equal_the_raw_csv_and_are_not_placeholders(run):
    raw = pd.read_csv(ROOT / "data" / "raw" / "pour_point" / "crc647_fuel_properties.csv", dtype={"fuel_code": str})
    raw["code"] = raw.fuel_code.astype(int)
    for r in run["dens"].itertuples():
        assert _raw_value_exists(raw, r.fuel_id, "density_kg_m3", r.measured_kg_m3, T_C=r.T_C), (r.fuel_id, r.T_C)
        assert r.measured_kg_m3 not in (795.0, 804.0)          # v3 placeholder constants (CLM-25)
    for r in run["fr"].itertuples():
        assert _raw_value_exists(raw, r.fuel_id, "freezing_point_C", r.measured_C), r.fuel_id
    for r in run["hoc"].itertuples():
        assert _raw_value_exists(raw, r.fuel_id, "net_heat_content_btu_lb", r.measured_ncv_MJ_kg,
                                 scale=0.002326, tol=0.01), r.fuel_id
    assert run["dens"].measured_kg_m3.nunique() > 30            # not a handful of repeated constants


def test_dev_only_offset_is_not_fit_on_locked(run):
    fr = run["fr"]
    dv = fr[fr.slice == "dev"]
    assert fr.dev_fit_offset_C.iloc[0] == pytest.approx((dv.wat_ideal_library_C - dv.measured_C).mean())
    lk = fr[fr.slice == "locked"]
    # changing locked targets must not move the offset: recompute from dev rows only
    assert (lk.wat_ideal_library_C - lk.measured_C).mean() != pytest.approx(fr.dev_fit_offset_C.iloc[0], abs=1e-6)


def test_dev_locked_follow_v4_splits_and_never_mix(run):
    sp = json.load(open(ev.SPLITS))
    for key in ("dens", "hoc", "fr"):
        df = run[key]
        assert set(df.slice) == {"dev", "locked"}
        for r in df.itertuples():
            assert (sp["fuel_role"][r.fuel_id] == "locked") == (r.slice == "locked")
    assert run["hoc"].slice.value_counts().to_dict() == {"dev": 45, "locked": 12}
    assert run["fr"].slice.value_counts().to_dict() == {"dev": 45, "locked": 12}
    assert "POSF-10325" not in set(run["hoc"].fuel_id)     # n=1 headline fuel: no statistics


def test_summary_reports_dev_and_locked_separately(run):
    s = ev.summarize(run["hoc"], {"m": ("err_v2_const_H", "measured_ncv_MJ_kg")})
    n = s.set_index("slice").N
    assert n["dev"] == 45 and n["locked"] == 12 and n["all_dev+locked"] == 57


def test_bias_stats_toy_case():
    meas = np.array([100.0, 200.0, 400.0, 50.0, 10.0])
    pred = meas + np.array([2.0, -4.0, 8.0, -1.0, 1.0])
    s = ev.bias_stats(pred, meas)
    assert s["N"] == 5
    assert s["MAE"] == pytest.approx(3.2)
    assert s["RMSE"] == pytest.approx(np.sqrt((4 + 16 + 64 + 1 + 1) / 5))
    assert s["bias"] == pytest.approx(1.2)
    assert s["MAPE_pct"] == pytest.approx(np.mean([2 / 100, 4 / 200, 8 / 400, 1 / 50, 1 / 10]) * 100)
    assert s["MAE_lo95"] <= s["MAE"] <= s["MAE_hi95"]
    assert ev.bias_stats([1.0, np.nan], [1.0, 2.0])["N"] == 1      # non-finite pairs dropped, not counted


def test_ideal_volume_density_toy():
    # equal-mass 1:1 of 700 and 900 kg/m3 -> harmonic mean
    w, rho = np.array([0.5, 0.5]), np.array([700.0, 900.0])
    assert 1.0 / np.sum(w / rho) == pytest.approx(2 * 700 * 900 / 1600)


def test_ideal_sle_mole_fraction_toy():
    tmap = {"A": dict(tm_K=263.6, dh_fus_J_mol=37400.0, n_carbon=12, cls="n_paraffin"),
            "B": dict(tm_K=200.0, dh_fus_J_mol=10000.0, n_carbon=8, cls="iso_paraffin")}
    # pure crystalliser: ln x = 0 -> liquidus equals its melting point
    assert ev.ideal_sle_mole_C([1.0], ["A"], [170.3], tmap) == pytest.approx(263.6 - 273.15)
    # dilution to x=0.3 (equal MW): 1/T = 1/Tm - (R/dH) ln x  (Schroeder-van Laar)
    t = ev.ideal_sle_mole_C([0.3, 0.7], ["A", "B"], [100.0, 100.0], tmap)
    assert t == pytest.approx(1.0 / (1 / 263.6 - ev.R_GAS / 37400.0 * np.log(0.3)) - 273.15)


def test_isotope_normalisation_and_class():
    from rdkit import Chem
    m, deut = ev.norm_mol("[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])[2H]")
    assert deut and ev.hc_class(m) == "n_paraffin"
    assert Chem.MolToSmiles(m) == "CCC"


def test_library_coutinho_is_the_one_used():
    from fuelblend.thermo import coutinho_wax
    assert ev.predict_freezing_point_coutinho is coutinho_wax.predict_freezing_point_coutinho
    assert ev.get_pure_thermal_properties is coutinho_wax.get_pure_thermal_properties


def test_bad_inputs_are_rejected():
    df = pd.DataFrame(dict(fuel_id=["CRC647_2"], target_provenance=["formula_derived"],
                           composition_provenance=["reconstructed"], split_role=["fold0"]))
    with pytest.raises(AssertionError):
        ev.assert_inputs_clean(df)
    df = df.assign(target_provenance="raw_file", composition_provenance="quarantined")
    with pytest.raises(AssertionError):
        ev.assert_inputs_clean(df)
    df = df.assign(composition_provenance="reconstructed", fuel_id="POSF-10325")
    with pytest.raises(AssertionError):
        ev.assert_inputs_clean(df)
