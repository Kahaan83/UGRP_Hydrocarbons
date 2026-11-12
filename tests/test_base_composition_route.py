"""v5 BASE Task 2: shared composition route (src/fuelblend/eval/composition_route.py)."""
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
V5 = ROOT / "data" / "real_fuels" / "v5"
S3 = ROOT / "results" / "v5" / "S3"
OUTDIR = ROOT / "results" / "v5" / "BASE"

from fuelblend.eval.composition_route import (  # noqa: E402
    COVERAGE_MIN, MappedFuel, RouteCtx, build_tmap, kind_bins, map_bins, predict_properties, refutas_nu_safe)

ALL_TARGETS_9T = pd.DataFrame(
    [dict(property="flash_point_C", T_K=np.nan), dict(property="freezing_point_C", T_K=np.nan)]
    + [dict(property=p, T_K=T) for p in ("viscosity_kin_cSt", "density_kg_m3")
       for T in (233.15, 243.15, 253.15, 263.15, 273.15, 283.15, 288.15, 298.15, 313.15)])


@pytest.fixture(scope="module")
def parts():
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    return reg, prov, reps


@pytest.fixture(scope="module")
def reps(parts):
    return parts[2]


@pytest.fixture(scope="module")
def route_ctx(parts):
    reg, prov, reps = parts
    idx = sorted({i for v in reps.by.values() for i in v})
    return RouteCtx(reg, prov, build_tmap(reg, idx))


@pytest.fixture(scope="module")
def propane_idx(parts):
    return parts[0].smiles.index("CCC")


@pytest.fixture(scope="module")
def decane_idx(parts):
    return parts[0].smiles.index("CCCCCCCCCC")


@pytest.fixture(scope="module")
def av3322_fuel_01(parts):
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    fid = sorted(comp.fuel_id.unique())[0]
    return map_bins(kind_bins("av3322", comp[comp.fuel_id == fid]), parts[2])


class _NoDiaromaticReps:
    """Stub with a representative for everything except aromatic_di (index 0)."""
    def get(self, kind, n_carbon, max_gap=3):
        return None if kind == "aromatic_di" else 0 if kind == "n_paraffin" else 1


def test_zero_mass_rep_does_not_poison_refutas(route_ctx, propane_idx, decane_idx):
    nu, rho = refutas_nu_safe(route_ctx.prov, np.array([propane_idx, decane_idx]), np.array([0.0, 1.0]), 313.15)
    assert np.isfinite(nu) and np.isfinite(rho)


def test_status_ok_implies_finite(route_ctx, av3322_fuel_01):
    out = predict_properties(route_ctx, av3322_fuel_01, ALL_TARGETS_9T)
    assert (out.status == "ok").any()
    assert np.isfinite(out.loc[out.status == "ok", "pred"]).all()


def test_map_bins_is_property_independent():
    bins = pd.DataFrame([dict(hc_class="n_paraffin", n_carbon=10, mass_pct=60.0),
                         dict(hc_class="iso_paraffin", n_carbon=11, mass_pct=37.0),
                         dict(hc_class="aromatic_di", n_carbon=11, mass_pct=3.0)])
    m = map_bins(bins, _NoDiaromaticReps())
    assert m.unmapped_by_class["aromatic_di"] == pytest.approx(0.03)
    assert m.unmapped_fraction == pytest.approx(0.03)
    assert m.w.sum() == pytest.approx(1.0)
    assert not m.low_coverage and COVERAGE_MIN == 0.90


def test_low_coverage_flag_when_over_ten_percent_unmapped():
    bins = pd.DataFrame([dict(hc_class="n_paraffin", n_carbon=10, mass_pct=85.0),
                         dict(hc_class="no_rep:tricycloalkane", n_carbon=12, mass_pct=15.0)])
    m = map_bins(bins, _NoDiaromaticReps())
    assert m.unmapped_fraction == pytest.approx(0.15) and m.low_coverage


def _score_like_s3(parts, route_ctx):
    reg, prov, reps = parts
    pred = pd.read_csv(S3 / "corrected_2026-11-12_predictions.csv")
    fuels = pd.read_csv(S3 / "corrected_2026-11-12_fuels.csv")
    comps = {"AV-33-22": pd.read_csv(V5 / "crc_av3322_composition.csv"),
             "FuelLib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "AJFTD": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    key = {"AV-33-22": "av3322", "FuelLib": "fuellib", "AJFTD": "ajftd"}
    rows = []
    for f in fuels[fuels.M_ok].itertuples():
        src = key[f.source]
        fid = int(f.fuel_id) if f.source == "AJFTD" else f.fuel_id
        c = comps[f.source]
        mf = map_bins(kind_bins(src, c[c.fuel_id == fid]), reps)
        t = pred[(pred.entity == f.entity) & pred.status_M.eq("ok")][["entity", "property", "T_K"]].drop_duplicates()
        if t.empty:
            continue
        t = t.assign(T_K=t.T_K.astype(float))
        out = predict_properties(route_ctx, mf, t)
        out["entity"] = f.entity
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


def test_reproduces_s3_corrected_M_predictions(parts, route_ctx):
    ref = pd.read_csv(S3 / "corrected_2026-11-12_predictions.csv")
    ref = ref[ref.status_M == "ok"].copy()
    new = _score_like_s3(parts, route_ctx).rename(columns={"pred": "pred_new", "status": "status_new"})
    ref["Tk"] = ref.T_K.round(3).fillna(-1.0)
    new["Tk"] = new.T_K.astype(float).round(3).fillna(-1.0)
    j = ref.merge(new[["entity", "property", "Tk", "pred_new", "status_new"]], on=["entity", "property", "Tk"], how="inner")
    assert len(j) > 100
    j["match"] = np.isclose(j.pred_M, j.pred_new, rtol=1e-5, atol=1e-6) & (j.status_new == "ok")
    bad = j[~j.match]
    ajftd_bad = bad[bad.source == "AJFTD"]
    assert ajftd_bad.empty, ajftd_bad[["entity", "property", "T_K", "pred_M", "pred_new"]].head()
    if len(bad):
        d = bad[["entity", "property", "T_K", "pred_M", "pred_new"]].rename(columns={"pred_M": "pred_s3"})
        d["cause"] = ("S2 maps for AV-33-22/FuelLib mark bins unmapped when the representative lacks Antoine or LFL "
                      "(flash-specific); BASE maps by SubReps.get only")
        OUTDIR.mkdir(parents=True, exist_ok=True)
        d.to_csv(OUTDIR / "route_vs_s3_differences.csv", index=False)
        assert d.cause.str.len().gt(0).all()


def test_lfl_override_changes_only_flash(route_ctx, av3322_fuel_01, parts):
    reg = parts[0]
    base = predict_properties(route_ctx, av3322_fuel_01, ALL_TARGETS_9T)
    ik = reg.ik[int(av3322_fuel_01.idx[np.argmax(av3322_fuel_01.w)])]
    old = route_ctx.flash_reg.loc[route_ctx.flash_reg.inchikey == ik, "lfl_vol_pct"].iloc[0]
    over = predict_properties(replace(route_ctx, lfl_override={ik: float(old) * 2.0}), av3322_fuel_01, ALL_TARGETS_9T)
    nf = (base.property != "flash_point_C").to_numpy()
    assert np.array_equal(base.pred.to_numpy()[nf], over.pred.to_numpy()[nf], equal_nan=True)
    assert base.pred.to_numpy()[~nf][0] != over.pred.to_numpy()[~nf][0]
