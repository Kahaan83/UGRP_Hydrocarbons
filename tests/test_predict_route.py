"""Task 1 (v5 api route): route context and per-property component mapping."""
import os

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import inchi

from fuelblend.api import predict_route as pr
from fuelblend.baselines.flash_unifac import _lookup
from fuelblend.eval.composition_route import kind_bins, map_bins

DECANE = "CCCCCCCCCC"
DODECANE = "CCCCCCCCCCCC"
TOLUENE = "Cc1ccccc1"


def _ik(smi):
    return inchi.MolToInchiKey(Chem.MolFromSmiles(smi))


@pytest.fixture(scope="module")
def ctx():
    return pr.get_ctx()


def test_resolve_smiles_and_inchikey_agree(ctx):
    a = pr.resolve({DECANE: 1.0})
    b = pr.resolve({_ik(DECANE): 1.0})
    assert len(a) == len(b) == 1 and a[0] == b[0]
    assert a[0].w == pytest.approx(1.0) and a[0].is_hc and a[0].n_carbon == 10
    assert a[0].subkind == "n_paraffin" and a[0].reg_idx is not None


def test_resolve_merges_duplicates(ctx):
    c = pr.resolve([(DECANE, 1.0), (_ik(DECANE), 1.0), ("C(CCCCCCCC)C", 2.0), (TOLUENE, 0.0), (DODECANE, 4.0)])
    assert len(c) == 2
    d = {x.n_carbon: x.w for x in c}
    assert d[10] == pytest.approx(0.5) and d[12] == pytest.approx(0.5)


def test_resolve_bad_smiles_raises():
    with pytest.raises(ValueError):
        pr.resolve({"not_a_smiles(((": 1.0})
    with pytest.raises(ValueError):
        pr.resolve({"AAAAAAAAAAAAAA-BBBBBBBBBB-N": 1.0})


def test_map_for_own_data_no_substitution(ctx):
    comps = pr.resolve({DECANE: 0.4, DODECANE: 0.4, TOLUENE: 0.2})
    p = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert p.substituted_mass == 0 and p.unmapped_mass == 0 and p.own_mass == pytest.approx(1.0)
    assert p.w.sum() == pytest.approx(1.0)


def test_map_for_substitutes_missing(ctx):
    odd = "CC(C)CCCC(C)c1ccccc1"          # hydrocarbon absent from the registry, so no grid
    comps = pr.resolve({DECANE: 0.5, odd: 0.5})
    odd_c = [c for c in comps if c.ik == _ik(odd)][0]
    assert odd_c.reg_idx is None and odd_c.is_hc
    p = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert len(p.substitutions) == 1
    a, b, gap = p.substitutions[0]
    assert a == odd_c.ik and b != a and gap <= 3
    assert p.substituted_mass == pytest.approx(0.5) and p.own_mass == pytest.approx(0.5)


def test_non_hc_unmapped(ctx):
    comps = pr.resolve({DECANE: 0.8, "CCO": 0.2})
    p = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert p.unmapped_mass == pytest.approx(0.2) and p.w.sum() == pytest.approx(1.0)


def test_map_for_route_equivalence(ctx):
    root = pr.ROOT
    fu = pd.read_csv(root / "results/v5/BASE/t2_fuels.csv")
    fu = fu[~fu.low_coverage]
    V5 = root / "data/real_fuels/v5"
    srcs = {"av3322": "crc_av3322_composition.csv", "fuellib": "fuellib_gcxgc_compositions.csv", "ajftd": "ajftd_gcxgc_composition.csv"}
    comps = {}
    n = 0
    for r in fu.itertuples():
        if r.source not in srcs:
            continue
        if r.source not in comps:
            comps[r.source] = pd.read_csv(V5 / srcs[r.source])
        fid = int(r.fuel_id) if r.source == "ajftd" else r.fuel_id
        c = comps[r.source]
        mf = map_bins(kind_bins(r.source, c[c.fuel_id == fid]), ctx.reps)
        cs = pr.resolve({ctx.reg.ik[int(i)]: float(w) for i, w in zip(mf.idx, mf.w)})
        p = pr.map_for(cs, pr.predicates(ctx)["HAS_FLASH"], ctx.reps)
        assert np.array_equal(np.sort(p.idx), np.sort(mf.idx))
        order = np.argsort(p.idx)
        assert np.allclose(p.w[order], mf.w[np.argsort(mf.idx)])
        n += 1
        if n == 5:
            break
    assert n == 5


def test_provider_never_writes():
    path = pr.GRID_NPZ
    m0 = os.stat(path).st_mtime_ns
    pr.get_ctx.cache_clear()
    c = pr.get_ctx()
    c.prov.lookup(np.array([[0]]), 300.0)
    c.prov.ensure([0, 1, 2])
    assert os.stat(path).st_mtime_ns == m0


def test_provider_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        pr.ReadOnlyProvider(pr.get_ctx().reg, tmp_path / "nope.npz")


def test_non_hc_with_grid_is_unmapped(ctx):
    etbe = "CCOC(C)(C)C"
    c = pr.resolve({DECANE: 0.8, etbe: 0.2})
    e = [x for x in c if not x.is_hc][0]
    if e.reg_idx is None or not ctx.prov.valid[e.reg_idx]:
        pytest.skip("ETBE has no grid in this registry")
    p = pr.map_for(c, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert p.unmapped_mass == pytest.approx(0.2) and p.own_mass == pytest.approx(0.8)


def test_rep_without_own_data_is_unmapped(ctx):
    P = pr.predicates(ctx)
    cand = [i for i in range(len(ctx.reg.ik)) if P["HAS_GRID"](i) and P["HAS_FLASH"](i) and not P["HAS_SLE"](i)
            and pr.subkind(ctx.reg.smiles[i]) is not None]
    if not cand:
        pytest.skip("no grid+flash compound lacking SLE")
    c = pr.resolve({ctx.reg.ik[cand[0]]: 1.0})
    p = pr.map_for(c, P["HAS_SLE"], ctx.reps)
    assert p.substituted_mass == 0 and p.unmapped_mass == pytest.approx(1.0) and p.idx.size == 0


def test_resolve_mole_fractions(ctx):
    c = pr.resolve({DECANE: 0.5, DODECANE: 0.5}, mole_fractions=True)
    d = {x.n_carbon: x.w for x in c}
    m10, m12 = ctx.mw(c[0].smiles if c[0].n_carbon == 10 else c[1].smiles), ctx.mw(c[0].smiles if c[0].n_carbon == 12 else c[1].smiles)
    assert d[10] == pytest.approx(m10 / (m10 + m12))


# ---------------------------------------------------------------- Task 2: predict_fuel and predict() on the route
import inspect
import math
import time

from fuelblend.api import predict_blend as pb
from fuelblend.baselines.flash_unifac import flash_point_unifac_detail
from fuelblend.eval.composition_route import refutas_nu_safe

MESITYLENE = "Cc1cc(C)cc(C)c1"
BLEND3 = {DECANE: 0.5, DODECANE: 0.3, TOLUENE: 0.2}
ISOOCTANE = "CC(C)CC(C)(C)C"
BUTYLCYCLOHEXANE = "CCCCC1CCCCC1"
# inside the real-fuel class-mix range (cyclo 0-0.709, n-paraffin 0-0.508, iso-paraffin 0.007-0.995) and inside the smoke regression domain
BLEND_IN_MIX = {DECANE: 0.2, DODECANE: 0.1, ISOOCTANE: 0.2, BUTYLCYCLOHEXANE: 0.3, TOLUENE: 0.2}
SMOKE_V5_MARGIN = 4.107069860688682          # configs/uncertainty_v5.yaml smoke_point_mm.ajftd_v5_in_domain.bulk_regression
CLASS_MIX_REASON = "class_mix_outside_real_fuel_range"
PROPS4 = ("viscosity_kin_cSt", "viscosity_dyn_mPas", "density_kg_m3", "flash_point_C")


def test_three_component_physics(ctx):
    r = pr.predict_fuel(BLEND3, T_K=313.15)
    assert set(r) >= {"T_K", "warnings", "components"}
    f = r["flash_point_C"]
    assert f["method"] == "ideal_le_chatelier_anchored" and f["status"] == "ok"
    assert f["coverage"]["substituted_mass"] == 0 and f["coverage"]["unmapped_mass"] == 0
    comps = pr.resolve(BLEND3)
    pm = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    nu, rho = refutas_nu_safe(ctx.prov, pm.idx, pm.w, 313.15)
    def _refutas(p):  # designated may be the ML value (user decision 2026-11-13); Refutas is then variants['refutas_raw']
        d = r[p]
        return d["variants"]["refutas_raw"]["value"] if d["method"] == "refutas+ml_pair_gn" else d["value"]
    assert _refutas("viscosity_kin_cSt") == pytest.approx(nu, abs=1e-9, rel=0)
    assert r["density_kg_m3"]["value"] == pytest.approx(rho, abs=1e-9, rel=0)
    assert _refutas("viscosity_dyn_mPas") == pytest.approx(nu * rho / 1000.0, abs=1e-9, rel=0)
    for p in PROPS4:
        d = r[p]
        assert set(d) >= {"value", "unit", "method", "status", "coverage", "lo", "hi", "margin", "in_domain",
                          "indicative_interval", "variants", "warnings"}
        assert d["status"] == "ok" and d["indicative_interval"] is None
        assert set(d["coverage"]) == {"own_mass", "substituted_mass", "unmapped_mass", "substitutions"}


def test_flash_overlay_designated(ctx):
    r = pr.predict_fuel(BLEND3, T_K=313.15)
    comps = pr.resolve(BLEND3)
    pm = pr.map_for(comps, pr.predicates(ctx)["HAS_FLASH"], ctx.reps)
    smi = [ctx.reg.smiles[i] for i in pm.idx]
    x = pm.w / np.array([ctx.mw(s) for s in smi])
    x = x / x.sum()
    ov = flash_point_unifac_detail(smi, x, registry=pr._overlay_ctx().flash_reg, use_unifac=False)["flash_C"]
    plain = flash_point_unifac_detail(smi, x, registry=ctx.flash_reg, use_unifac=False)["flash_C"]
    f = r["flash_point_C"]
    assert f["value"] == pytest.approx(ov, abs=1e-9, rel=0)
    assert f["variants"]["v4_designated"] == pytest.approx(plain, abs=1e-9, rel=0)


def test_mesitylene_warning():
    r = pr.predict_fuel({DECANE: 0.5, DODECANE: 0.4, MESITYLENE: 0.1}, T_K=298.15, fuel_type="jet")
    f = r["flash_point_C"]
    assert "antoine_suspect:mesitylene" in f["warnings"]
    assert np.isfinite(f["lo"]) and "lab_margin_set_by_one_system" in f["warnings"]
    r2 = pr.predict_fuel(BLEND3, T_K=298.15, fuel_type="jet")
    assert "antoine_suspect:mesitylene" not in r2["flash_point_C"]["warnings"]
    assert "lab_margin_set_by_one_system" in r2["flash_point_C"]["warnings"]


def test_flash_independent_of_T():
    vals = [pr.predict_fuel(BLEND3, T_K=T)["flash_point_C"]["value"] for T in (233.15, 288.15, 313.15)]
    assert vals[0] == vals[1] == vals[2] and math.isfinite(vals[0])


def test_sixty_component_blend_runs_physics(ctx):
    P = pr.predicates(ctx)
    cand = [i for i in range(len(ctx.reg.ik)) if (P["HAS_GRID"](i) or P["HAS_FLASH"](i)) and pr.subkind(ctx.reg.smiles[i])
            and 8 <= pr._mol_info(ctx.reg.smiles[i])[2] <= 16 and pr._mol_info(ctx.reg.smiles[i])[3]]
    assert len(cand) >= 60
    pick = np.random.default_rng(0).choice(cand, 60, replace=False)
    comp = {ctx.reg.ik[int(i)]: 1.0 for i in pick}
    t0 = time.time()
    r = pr.predict_fuel(comp, T_K=298.15)
    assert time.time() - t0 < 20
    for p in PROPS4:
        d = r[p]
        assert d["method"] != "none" and d["status"] == "ok" and d["value"] is not None and math.isfinite(d["value"])


def _unregistered_hcs(ctx, n):
    P = pr.predicates(ctx)
    out = []
    for nc in range(9, 16):
        for a in range(1, nc // 2):
            for b in range(a + 1, nc - 1):
                s = "".join("C(C)" if k in (a, b) else "C" for k in range(nc))
                m = Chem.MolFromSmiles(s)
                if m is None or inchi.MolToInchiKey(m) in ctx.reg.ik2i:
                    continue
                c = pr.resolve({s: 1.0})[0]
                if (pr.map_for([c], P["HAS_GRID"], ctx.reps).substituted_mass == 1.0
                        and pr.map_for([c], P["HAS_FLASH"], ctx.reps).substituted_mass == 1.0):
                    out.append(s)
                if len(out) == n:
                    return out
    return out


def test_all_substituted_blend(ctx):
    smis = _unregistered_hcs(ctx, 60)
    assert len(smis) == 60
    r = pr.predict_fuel({s: 1.0 for s in smis}, T_K=298.15, fuel_type="jet")
    for p in PROPS4:
        d = r[p]
        assert d["coverage"]["substituted_mass"] == pytest.approx(1.0) and d["in_domain"] is False
        assert d["value"] is not None and math.isfinite(d["value"]) and d["status"] == "ok"
        assert math.isnan(d["lo"]) and math.isnan(d["hi"]) and d["margin"] is None


def test_duplicates_and_zero_weights():
    messy = pr.predict_fuel([(DECANE, 1.0), (_ik(DECANE), 2.0), ("C(CCCCCCCC)C", 0.0), (TOLUENE, 0.0), (DODECANE, 4.0)],
                            T_K=298.15)
    ref = pr.predict_fuel({DECANE: 0.5, DODECANE: 0.5}, T_K=298.15)
    assert len(messy["components"]) == 2
    for p in PROPS4:
        assert math.isfinite(messy[p]["value"])
    ref = pr.predict_fuel({DECANE: 3.0 / 7.0, DODECANE: 4.0 / 7.0}, T_K=298.15)
    for p in PROPS4:
        assert messy[p]["value"] == pytest.approx(ref[p]["value"], rel=1e-12)


def test_T_outside_grid():
    for T in (200.0, 500.0):
        r = pr.predict_fuel(BLEND3, T_K=T)
        for p in ("viscosity_kin_cSt", "viscosity_dyn_mPas", "density_kg_m3"):
            assert r[p]["value"] is None and r[p]["status"] == "T_out_of_grid" and r[p]["in_domain"] is False
        assert r["flash_point_C"]["status"] == "ok" and math.isfinite(r["flash_point_C"]["value"])


def test_non_hc_component():
    r5 = pr.predict_fuel({DECANE: 0.55, DODECANE: 0.4, "CCO": 0.05}, T_K=298.15, fuel_type="jet")
    assert any("non_hc" in w for w in r5["warnings"])
    for p in ("viscosity_kin_cSt", "density_kg_m3", "flash_point_C"):
        d = r5[p]
        assert d["status"] == "ok" and d["coverage"]["unmapped_mass"] == pytest.approx(0.05)
        assert d["in_domain"] is False and math.isnan(d["lo"]) and d["margin"] is None
    r20 = pr.predict_fuel({DECANE: 0.4, DODECANE: 0.4, "CCO": 0.2}, T_K=298.15)
    for p in ("viscosity_kin_cSt", "density_kg_m3", "flash_point_C"):
        assert r20[p]["status"] == "low_coverage" and r20[p]["value"] is None


def test_bad_smiles_raises():
    with pytest.raises(ValueError):
        pr.predict_fuel({"not_a_smiles(((": 1.0})


def test_no_gbm_mix_anywhere():
    assert "gbm_mix" not in inspect.getsource(pb)
    assert not hasattr(pb, "_flash_lechatelier")


def test_predict_wrapper_keys():
    r = pb.predict([(DECANE, 0.5), (DODECANE, 0.5)], T_K=298.15)
    assert {"viscosity_dyn_mPas", "density_kg_m3", "flash_point_C", "uncertified", "T_K", "warnings"} <= set(r)
    for p in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"):
        assert isinstance(r[p]["value"], float) and "tuple" in r[p] and "model_name" in r[p]
    assert r["flash_point_C"]["method"] == "ideal_le_chatelier_anchored"
    f = pr.predict_fuel({DECANE: 0.5, DODECANE: 0.5}, T_K=298.15)["flash_point_C"]["value"]
    assert r["flash_point_C"]["value"] == f


# ====================================================================================================================
# Task 3: freezing, heat, smoke, aromatics; K1 and F1x variants
# ====================================================================================================================
from fuelblend.api import predict_route as _pr3
from fuelblend.eval import freezing_f1 as ff
from fuelblend.eval.composition_route import ideal_sle_C, build_tmap, build_tmap_designated


def _t2_fuels(n, need_heat=False):
    """First n non-low-coverage T2 fuels as ({inchikey: w}, MappedFuel, entity)."""
    root = pr.ROOT
    fu = pd.read_csv(root / "results/v5/BASE/t2_fuels.csv")
    fu = fu[~fu.low_coverage]
    if need_heat:
        rows = pd.read_csv(root / "results/v5/BASE/rows_t2_bulk.csv")
        h = rows[(rows.property == "heat_of_combustion_net_MJ_kg") & (rows.model == "ga_comp_H") & (rows.status == "ok")]
        fu = fu[fu.entity.isin(set(h.entity))]
    V5 = root / "data/real_fuels/v5"
    srcs = {"av3322": "crc_av3322_composition.csv", "fuellib": "fuellib_gcxgc_compositions.csv", "ajftd": "ajftd_gcxgc_composition.csv"}
    comps, out = {}, []
    c0 = pr.get_ctx()
    for r in fu.itertuples():
        if r.source not in srcs:
            continue
        if r.source not in comps:
            comps[r.source] = pd.read_csv(V5 / srcs[r.source])
        fid = int(r.fuel_id) if r.source == "ajftd" else r.fuel_id
        c = comps[r.source]
        mf = map_bins(kind_bins(r.source, c[c.fuel_id == fid]), c0.reps)
        out.append(({c0.reg.ik[int(i)]: float(w) for i, w in zip(mf.idx, mf.w)}, mf, r.entity))
        if len(out) == n:
            break
    return out


def test_freezing_svl_matches_route(ctx):
    reg = ctx.reg
    tmap = build_tmap_designated(build_tmap(reg, [i for i in range(len(reg.ik)) if reg.usable[i]]))   # verified where present, else Coutinho
    fuels = _t2_fuels(5)
    assert len(fuels) == 5
    n_finite = 0
    for comp, mf, _ in fuels:
        iks = [reg.ik[int(i)] for i in mf.idx]
        x = mf.w / np.array([ctx.mw(reg.smiles[int(i)]) for i in mf.idx])
        x = x / x.sum()
        ref = ideal_sle_C(x, iks, tmap)
        r = pr.predict_fuel(comp, T_K=298.15, fuel_type="jet")["freezing_point_C"]
        if np.isfinite(ref):
            n_finite += 1
            assert r["value"] is not None and abs(r["value"] - ref) < 1e-9
        else:
            assert r["value"] is None
    assert n_finite >= 3


def test_f1x_coefficients_reproduce():
    a, b = ff.fit_f1x()
    assert np.isfinite(a) and np.isfinite(b)
    s5 = pd.read_csv(pr.ROOT / "results/v5/S5/s5_predictions.csv")
    h = s5[(s5.role == "headline") & s5.F0x.notna() & s5.freezing_C.notna()].reset_index(drop=True)
    assert len(h) == 49
    for i in (0, 10, 30):
        tr = h[h.group != h.group[i]]
        bb, aa = np.polyfit(tr.F0x, tr.freezing_C, 1)
        assert abs((aa + bb * h.F0x[i]) - h.F1x[i]) < 1e-6
    d = pd.read_csv(pr.ROOT / "results/v5/BASE/delta_freezing/delta.csv")
    row = d[(d.set == "pooled") & (d.model == "F1x_coutinho_complete")].iloc[0]
    assert ff.F1X_LOGO_MAE == float(row.mae) and 3.0 < ff.F1X_LOGO_MAE < 3.1


def test_f1x_variant_jet_with_nalkanes_only():
    blend = {DECANE: 0.5, DODECANE: 0.3, TOLUENE: 0.2}
    r = pr.predict_fuel(blend, fuel_type="jet")["freezing_point_C"]
    v = r["variants"]["f1x_recalibrated"]
    a, b = ff.fit_f1x()
    assert abs(v["value"] - (a + b * v["T_first_C"])) < 1e-9
    assert v["logo_mae_C"] == ff.F1X_LOGO_MAE
    assert "f1x_recalibrated" not in pr.predict_fuel(blend, fuel_type=None)["freezing_point_C"]["variants"]
    assert "f1x_recalibrated" not in pr.predict_fuel({TOLUENE: 1.0}, fuel_type="jet")["freezing_point_C"]["variants"]


def test_heat_matches_bulk_rows(ctx):
    rows = pd.read_csv(pr.ROOT / "results/v5/BASE/rows_t2_bulk.csv")
    h = rows[(rows.property == "heat_of_combustion_net_MJ_kg") & (rows.model == "ga_comp_H") & (rows.status == "ok")]
    fuels = _t2_fuels(5, need_heat=True)
    assert len(fuels) == 5
    for comp, mf, ent in fuels:
        ref = h[h.entity == ent].pred.iloc[0]
        r = pr.predict_fuel(comp, fuel_type="jet")["heat_net_MJ_kg"]
        assert r["value"] is not None and abs(r["value"] - ref) < 5e-7, (ent, r["value"], ref)   # csv stores 8 significant digits (1e-9 is unreachable from the file)
        assert r["method"] == "ga_comp_H"
    r = pr.predict_fuel({DECANE: 0.5, DODECANE: 0.5}, fuel_type="jet")["heat_net_MJ_kg"]
    assert np.isfinite(r["lo"]) and r["in_domain"] is True and r["margin"] > 0
    assert pr.predict_fuel({DECANE: 1.0}, fuel_type=None)["heat_net_MJ_kg"]["value"] is None


def test_smoke_designated_bulk_regression(ctx):
    P = pr.predicates(ctx)
    reg = ctx.reg
    no_ysi = next(reg.ik[i] for i in range(len(reg.ik)) if reg.usable[i] and not P["HAS_YSI"](i) and pr._mol_info(reg.smiles[i])[3]
                  and pr.subkind(reg.smiles[i]) is not None)
    r = pr.predict_fuel({TOLUENE: 0.5, no_ysi: 0.5})["smoke_point_mm"]
    # designated = bulk regression, independent of YSI coverage; YSI only a gated variant
    assert r["method"] == "ajftd_bulk_regression" and r["status"] == "ok" and np.isfinite(r["value"])
    assert r["variants"]["ysi_route"]["status"] == "no_evaluation" and r["variants"]["ysi_route"]["value"] is None
    assert r["coverage"]["ysi_own_mass"] < 0.9
    assert math.isnan(r["lo"]) and math.isnan(r["hi"]) and r["in_domain"] is False and r["uncertified"] is True
    assert any(w.startswith("out_of_domain") for w in r["warnings"])        # 50 % aromatics is far outside 8-23 vol %
    ok = pr.predict_fuel({TOLUENE: 0.1, DECANE: 0.9})["smoke_point_mm"]
    # user decision R-SC-2b ('Yes, replace (Recommended)'): the real-fuel margin applies only inside the class-mix range; decane + toluene has
    # no iso-paraffin (range 0.007-0.995), so it is in the aromatics domain but uncertified with the class-mix reason
    assert ok["method"] == "ajftd_bulk_regression" and ok["status"] == "ok"
    assert not any(w.startswith("out_of_domain") for w in ok["warnings"])
    assert ok["in_domain"] is False and ok["uncertified"] is True and any(CLASS_MIX_REASON in x for x in ok["domain_reasons"])
    assert ok["variants"]["ysi_route"]["status"] == "ok" and np.isfinite(ok["variants"]["ysi_route"]["value"])


def test_smoke_in_domain_no_warning():
    """User decision R-SC-2b, 'Yes, replace (Recommended)': in the aromatics domain AND inside the class-mix range the real-fuel margin
    4.107069860688682 mm replaces the 7.8 mm lab margin; in the aromatics domain but outside the class-mix range it is uncertified."""
    r = pr.predict_fuel(BLEND_IN_MIX, fuel_type="jet")
    s = r["smoke_point_mm"]
    a = r["aromatics_vol_pct"]["value"]
    assert 8.0 <= a <= 23.0 and s["method"] == "ajftd_bulk_regression" and s["status"] == "ok"
    assert s["inputs"]["aromatics_vol_pct"] == pytest.approx(a)
    assert not any(w.startswith("out_of_domain") for w in s["warnings"])
    assert s["margin"] == SMOKE_V5_MARGIN and s["in_domain"] is True
    r2 = pr.predict_fuel({DECANE: 0.85, TOLUENE: 0.15}, fuel_type="jet")["smoke_point_mm"]     # in aromatics domain, no iso-paraffin
    assert 8.0 <= pr.predict_fuel({DECANE: 0.85, TOLUENE: 0.15})["aromatics_vol_pct"]["value"] <= 23.0
    assert r2["in_domain"] is False and r2["uncertified"] is True and any(CLASS_MIX_REASON in x for x in r2["domain_reasons"])
    assert math.isnan(r2["lo"]) and math.isnan(r2["hi"])


def test_aromatics_by_construction():
    r = pr.predict_fuel({DECANE: 0.8, TOLUENE: 0.2})
    assert r["aromatics_wt_pct"]["value"] == pytest.approx(20.0, abs=1e-9)
    assert r["aromatics_wt_pct"]["method"] == "by_construction"
    v = r["aromatics_vol_pct"]
    assert v["method"] == "by_construction" and 15.0 < v["value"] < 20.0      # toluene is denser than decane -> its volume share is below its mass share


def test_k1_only_at_fitted_T():
    a = pr.predict_fuel(BLEND3, T_K=233.15)["viscosity_kin_cSt"]
    k = a["variants"]["k1_cold_offset_in_sample"]
    assert k["label"] == "in_sample_T2" and np.isfinite(k["value"]) and k["value"] != a["value"]
    base = a["variants"]["refutas_raw"]["value"] if a["method"] == "refutas+ml_pair_gn" else a["value"]
    assert abs(k["value"] - base * math.exp(k["offset_ln"])) < 1e-9   # K1 is an offset on raw Refutas
    assert "k1_cold_offset_in_sample" not in pr.predict_fuel(BLEND3, T_K=240.0)["viscosity_kin_cSt"]["variants"]


def test_uncertified_props_have_nan_interval():
    r = pr.predict_fuel(BLEND3, T_K=298.15, fuel_type="jet")
    assert math.isnan(r["freezing_point_C"]["lo"]) and math.isnan(r["freezing_point_C"]["hi"]) and r["freezing_point_C"]["in_domain"] is False
    assert r["freezing_point_C"]["uncertified"] is True
    w = pb.predict([(DECANE, 0.5), (DODECANE, 0.5)], fuel_type="jet")
    assert "freezing_point_C" not in pb.UNCERTIFIED_PROPERTIES and "smoke_point_mm" not in pb.UNCERTIFIED_PROPERTIES
    assert w["smoke_point_mm"]["in_domain"] is False and math.isnan(w["smoke_point_mm"]["lo"])
    assert w["uncertified"]["smoke_point_mm"]["in_domain"] is False
    # v5 delta_freeze_margin: decane + dodecane are both in tm_dhfus_verified.csv, so the lab-liquidus margin applies (known composition)
    fz = w["freezing_point_C"]
    assert fz["in_domain"] is True and fz["margin"] == pytest.approx(16.975450975521966) and fz["safe_bound"] == fz["hi"]


def test_k1_not_on_reconstructed():
    r = pr.predict_fuel(BLEND3, T_K=233.15, composition_kind="reconstructed")["viscosity_kin_cSt"]
    assert "k1_cold_offset_in_sample" not in r["variants"]


def test_aromatics_same_basis_with_non_hc():
    r = pr.predict_fuel({DECANE: 0.4, TOLUENE: 0.4, "CCO": 0.2})
    assert r["aromatics_wt_pct"]["value"] == pytest.approx(50.0, abs=1e-9)
    assert r["heat_net_MJ_kg"]["tuple"] is not None


def test_indicative_flash_q90():
    f = pd.read_csv(pr.ROOT / "results" / "v5" / "BASE" / "rows_t2_flash.csv")
    f = f[(f.status == "ok") & (~f.low_coverage.astype(bool)) & f.chosen.astype(bool)].drop_duplicates("group")
    q = float(np.quantile(np.abs(f.pred_designated - f.measured), 0.9))
    r = pr.predict_fuel(BLEND3, fuel_type="jet")["flash_point_C"]
    ii = r["indicative_interval"]
    assert ii["q90_abs"] == pytest.approx(q) and ii["n_groups"] == len(f)
    assert ii["kind"] == "empirical_T2_descriptive_not_held_out"
    assert ii["lo"] == pytest.approx(r["value"] - q) and ii["hi"] == pytest.approx(r["value"] + q)


def test_indicative_absent_without_jet():
    r = pr.predict_fuel(BLEND3)
    assert r["flash_point_C"]["value"] is not None and r["flash_point_C"]["indicative_interval"] is None
    assert all(v["indicative_interval"] is None for v in r.values() if isinstance(v, dict) and "indicative_interval" in v)
    r = pr.predict_fuel({DECANE: 0.8, "CCO": 0.2}, fuel_type="jet")          # unmapped 0.2 > 0.10
    assert r["flash_point_C"]["value"] is None and r["flash_point_C"]["indicative_interval"] is None


def test_lo_hi_untouched_by_indicative(monkeypatch):
    a = pr.predict_fuel(BLEND3, T_K=253.15, fuel_type="jet")
    monkeypatch.setattr(pr, "_indicative_interval", lambda *x: None)
    b = pr.predict_fuel(BLEND3, T_K=253.15, fuel_type="jet")
    assert b["flash_point_C"]["indicative_interval"] is None
    for k in ("flash_point_C", "viscosity_kin_cSt", "density_kg_m3", "heat_net_MJ_kg", "smoke_point_mm"):
        for f in ("lo", "hi", "in_domain"):
            x, y = a[k][f], b[k][f]
            assert (x == y) or (isinstance(x, float) and math.isnan(x) and math.isnan(y))
    assert a["viscosity_kin_cSt"]["indicative_interval"]["lo"] > 0


def _bulk(prop, T):
    b = pd.read_csv(pr.ROOT / "results" / "v5" / "BASE" / "rows_t2_bulk.csv")
    b = b[(b.property == prop) & (b.variant == "designated") & (b.status == "ok") & (~b.low_coverage.astype(bool)) & (b.T_K == T)]
    return b.sort_values("chosen", ascending=False, kind="stable").drop_duplicates("group")


def test_indicative_density_q90_from_csv():
    b = _bulk("density_kg_m3", 313.15)
    q = float(np.quantile(np.abs(b.pred - b.measured), 0.9))
    ii = pr._indicative_interval("density_kg_m3", 313.15, 800.0)
    assert ii["q90_abs"] == pytest.approx(q) and ii["n_groups"] == len(b) == b.group.nunique()
    assert ii["lo"] == pytest.approx(800.0 - q) and ii["hi"] == pytest.approx(800.0 + q)


def test_indicative_viscosity_multiplicative():
    b = _bulk("viscosity_kin_cSt", 253.15)
    q = float(np.quantile(np.abs(np.log(b.pred) - np.log(b.measured)), 0.9))
    ii = pr._indicative_interval("viscosity_kin_cSt", 253.15, 4.0)
    assert ii["q90_abs"] == pytest.approx(q)
    assert ii["lo"] == pytest.approx(4.0 * math.exp(-q)) and ii["hi"] == pytest.approx(4.0 * math.exp(q))
    assert ii["lo"] > 0


# ---------------------------------------------------------------------------------------------- whole-branch review fixes
def test_flash_window_exceeded_not_in_domain():
    """I1: a flash outside the nominal 0-100 C Antoine window keeps its value and warning but gets no interval."""
    for comp in ({"CCCCCCCCCCCCCC": 0.5, "CCCCCCCCCCCCCCCC": 0.5}, {"CCCCCC": 0.5, "CCCCCCC": 0.5}):
        f = pr.predict_fuel(comp, T_K=298.15, fuel_type="jet")["flash_point_C"]
        assert f["value"] is not None and not (0.0 <= f["value"] <= 100.0)
        assert any(w.startswith("antoine_nominal_window_exceeded") for w in f["warnings"])
        assert np.isnan(f["lo"]) and np.isnan(f["hi"]) and f["in_domain"] is False
        assert "antoine_nominal_window_exceeded" in f["domain_reasons"]
        assert np.isnan(f["tuple"][1]) and f["tuple"][4] is False
    g = pr.predict_fuel(BLEND3, T_K=298.15, fuel_type="jet")["flash_point_C"]          # inside the window: unchanged
    assert 0.0 <= g["value"] <= 100.0 and np.isfinite(g["lo"])


def test_mesitylene_warning_as_representative(monkeypatch):
    """M5: the warning also fires when mesitylene enters the flash mapping as a representative."""
    real = pr.map_for

    def fake(comps, has_own, reps):
        pm = real(comps, has_own, reps)
        pm.substitutions.append((comps[0].ik, pr.MESITYLENE_IK, 0))
        return pm
    monkeypatch.setattr(pr, "map_for", fake)
    f = pr.predict_fuel(BLEND3, T_K=298.15, fuel_type="jet")["flash_point_C"]
    assert "antoine_suspect:mesitylene" in f["warnings"]


def test_indicative_needs_five_groups_and_near_bin():
    """I2: no interval from a bin with fewer than 5 groups or one more than 10 K from T."""
    tab = pr._indicative_tables()
    rows = tab["viscosity_kin_cSt"]
    small = [r for r in rows if r[2] < pr.INDICATIVE_MIN_GROUPS and np.isfinite(r[0])]
    big = [r for r in rows if r[2] >= pr.INDICATIVE_MIN_GROUPS and np.isfinite(r[0])]
    assert small and big
    assert pr._indicative_interval("viscosity_kin_cSt", small[0][0], 4.0) is None
    ii = pr._indicative_interval("viscosity_kin_cSt", big[0][0], 4.0)
    assert ii["T_bin_K"] == big[0][0] and ii["n_groups"] >= 5
    far = max(r[0] for r in rows if np.isfinite(r[0])) + 50.0
    assert pr._indicative_interval("viscosity_kin_cSt", far, 4.0) is None
    r = pr.predict_fuel(BLEND3, T_K=298.15, fuel_type="jet")["viscosity_kin_cSt"]       # n_groups 1 at 298.15 K
    assert r["indicative_interval"] is None


def test_indicative_skipped_for_recalibrated_refutas():
    r = pr.predict_fuel(BLEND3, T_K=253.15, fuel_type="jet", composition_kind="reconstructed")["viscosity_kin_cSt"]
    assert r["method"] == "refutas_recalibrated" and r["indicative_interval"] is None


def test_grid_sha_mismatch_warns_not_raises(monkeypatch):
    """I4: a grid that differs from the pinned SHA-256 gives a top-level warning, no exception."""
    if not pr.GRID_NPZ.exists():
        pytest.skip("grid absent")
    assert pr.grid_sha_status()[0] is True
    assert not any(w.startswith("grid_sha_mismatch") for w in pr.predict_fuel(BLEND3)["warnings"])
    monkeypatch.setattr(pr, "GRID_SHA256", "0" * 64)
    pr.grid_sha_status.cache_clear()
    try:
        r = pr.predict_fuel(BLEND3)
        assert any(w.startswith("grid_sha_mismatch") for w in r["warnings"]) and r["flash_point_C"]["value"] is not None
    finally:
        monkeypatch.undo()
        pr.grid_sha_status.cache_clear()


def test_route_library_uses_read_only_provider():
    from fuelblend.design import library_route as LR
    L = LR.RouteLibrary(["n-decane", "n-dodecane"])
    assert isinstance(L.provider, pr.ReadOnlyProvider)


def test_predict_returns_viscosity_kin():
    from fuelblend.api import predict_blend as pb
    r = pb.predict([(DECANE, 0.5), (DODECANE, 0.5)], T_K=233.15)
    assert r["viscosity_kin_cSt"]["value"] > 0 and "k1_cold_offset_in_sample" in r["viscosity_kin_cSt"]["variants"]


def test_smoke_margin_in_domain_user_decision():
    """User decision R-SC-2b, 'Yes, replace (Recommended)': in-domain recipe inside the class-mix range gets the real-fuel margin
    4.107069860688682 mm (configs/uncertainty_v5.yaml smoke_point_mm.ajftd_v5_in_domain.bulk_regression), replacing the 7.8 mm lab margin;
    lo = value - margin, safe_bound = lo.  The lab key stays in the yaml (used when the new key is absent, see test_v5_margin_fn)."""
    s = pr.predict_fuel(BLEND_IN_MIX, T_K=288.15, fuel_type="jet")["smoke_point_mm"]
    assert 19.0 <= s["value"] <= 25.0 and s["in_domain"] is True
    assert s["margin"] == SMOKE_V5_MARGIN and s["margin_unit"] == "mm" and s["margin_direction"] == "over"
    assert s["lo"] == pytest.approx(s["value"] - SMOKE_V5_MARGIN, abs=1e-12)
    assert s["hi"] == pytest.approx(s["value"] + SMOKE_V5_MARGIN, abs=1e-12)
    assert s["safe_bound"] == s["lo"] and not s.get("uncertified")
    assert pr.SMOKE_MARGIN_NOTE_V5 in s["warnings"] and pr.SMOKE_MARGIN_NOTE not in s["warnings"]     # the note names the margin applied
    out = pr.predict_fuel(BLEND3, T_K=288.15, fuel_type="jet")["smoke_point_mm"]      # decane/dodecane/toluene: no iso-paraffin
    assert out["uncertified"] is True and any(CLASS_MIX_REASON in x for x in out["domain_reasons"]) and math.isnan(out["lo"])
    import yaml
    y = yaml.safe_load((pr.ROOT / "configs/uncertainty_v5.yaml").read_text(encoding="utf-8"))["margins"]
    assert y["smoke_point_mm.lab_mixtures.bulk_regression"]["applied_margin"] == 7.8009965628100595
    assert y["smoke_point_mm.ajftd_v5_in_domain.bulk_regression"]["applied_margin"] == SMOKE_V5_MARGIN
    assert y["flash_point_C.lab_mixtures.ideal"]["applied_margin"] == 2.206932808388615


def test_smoke_margin_not_applied_out_of_domain():
    """R-SM-2: aromatics outside 8-23 vol % (or result outside 19-25 mm) keeps NaN lo/hi, in_domain False, uncertified."""
    for blend in ({DECANE: 0.5, TOLUENE: 0.5}, {DECANE: 0.6, DODECANE: 0.4}):
        s = pr.predict_fuel(blend, T_K=288.15, fuel_type="jet")["smoke_point_mm"]
        assert math.isnan(s["lo"]) and math.isnan(s["hi"]) and s["in_domain"] is False and s["uncertified"] is True
        assert s["margin"] is None and any(w.startswith("out_of_domain") for w in s["warnings"])
        assert "smoke_margin_from_binaries_and_surrogates (user decision)" not in s["warnings"] and pr.SMOKE_MARGIN_NOTE_V5 not in s["warnings"]


def test_flash_margin_q2_user_decision():
    """R-Q2-1: in-domain known-composition flash margin is the Q2 value (35 systems, mesitylene excluded), lo = value - margin."""
    f = pr.predict_fuel(BLEND3, T_K=298.15, fuel_type="jet")["flash_point_C"]
    assert f["margin"] == 2.206932808388615
    assert f["lo"] == pytest.approx(f["value"] - 2.206932808388615, abs=1e-12)
    assert f["hi"] == pytest.approx(f["value"] + 2.206932808388615, abs=1e-12)
    assert f["in_domain"] is True
    assert "lab_margin_set_by_one_system" in f["warnings"]
    assert "flash_margin_excludes_mesitylene (user decision Q2)" in f["warnings"]
    import yaml
    y = yaml.safe_load((pr.ROOT / "configs/uncertainty_v5.yaml").read_text(encoding="utf-8"))["margins"]
    assert y["flash_point_C.lab_mixtures.ideal"]["applied_margin"] == 2.206932808388615
    assert y["flash_point_C.lab_mixtures.ideal"]["applied_margin_v4"] == 3.698768343688585
