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
    assert r["viscosity_kin_cSt"]["value"] == pytest.approx(nu, abs=1e-9, rel=0)
    assert r["density_kg_m3"]["value"] == pytest.approx(rho, abs=1e-9, rel=0)
    assert r["viscosity_dyn_mPas"]["value"] == pytest.approx(nu * rho / 1000.0, abs=1e-9, rel=0)
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
