"""v5 S3: group dedup, same property call for M and R, R-route input whitelist (scripts/v5/s3_composition_vs_physics.py)."""
import importlib.util
import inspect
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


s3 = _load("s3_composition_t", "scripts/v5/s3_composition_vs_physics.py")


def test_group_dedup_keeps_one_per_group():
    ent = pd.DataFrame([
        dict(entity="FuelLib:a", group="g1", source="FuelLib", n_pair_rows=5, pairable=True, group_locked=False),
        dict(entity="AJFTD:1", group="g1", source="AJFTD", n_pair_rows=4, pairable=True, group_locked=False),
        dict(entity="AJFTD:2", group="g1", source="AJFTD", n_pair_rows=5, pairable=True, group_locked=False),   # tie on rows: FuelLib wins by priority
        dict(entity="AV:1", group="g2", source="AV-33-22", n_pair_rows=3, pairable=True, group_locked=False),
        dict(entity="AJFTD:9", group="g3", source="AJFTD", n_pair_rows=9, pairable=False, group_locked=False),  # not pairable
        dict(entity="AJFTD:7", group="g4", source="AJFTD", n_pair_rows=6, pairable=True, group_locked=True),    # locked group
    ])
    sel = s3.select_one_per_group(ent)
    assert sel.group.is_unique
    assert set(sel.entity) == {"FuelLib:a", "AV:1"}
    ent.loc[ent.entity == "FuelLib:a", "n_pair_rows"] = 2                                                       # most rows now wins
    assert s3.select_one_per_group(ent).set_index("group").loc["g1", "entity"] == "AJFTD:2"


def test_union_find_groups():
    uf = s3.UF(["a", "b", "c", "d"])
    uf.union("a", "b"); uf.union("c", "b")
    assert uf.find("a") == uf.find("c") != uf.find("d")


def _toy_ctx():
    reg = SimpleNamespace(ik=["IK_C10", "IK_C12"], smiles=["CCCCCCCCCC", "CCCCCCCCCCCC"])
    refutas = lambda prov, idx, w, T: (2.0 + float(np.dot(w, np.asarray(idx) + 1)), 800.0, None, None)   # noqa: E731
    return s3.Ctx(reg, None, refutas, lambda x, iks, tmap: float(-50 + 10 * x[0]), {})


def test_m_and_r_use_same_property_call(monkeypatch):
    # structural: score_entity reaches the models only through predict_comp
    src = inspect.getsource(s3.score_entity)
    assert "predict_comp(" in src
    for banned in ("flash_point_unifac_detail", "refutas_nu", "ideal_sle_C"):
        assert banned not in src
    calls = []
    real = s3.predict_comp
    monkeypatch.setattr(s3, "predict_comp", lambda *a, **k: (calls.append(a[1:3]), real(*a, **k))[1])
    ctx = _toy_ctx()
    t = pd.DataFrame([dict(property="viscosity_kin_cSt", T_K=253.15, measured=3.0, target="visc_-20C"),
                      dict(property="density_kg_m3", T_K=288.15, measured=800.0, target="density_15C"),
                      dict(property="freezing_point_C", T_K=np.nan, measured=-50.0, target="freezing"),
                      dict(property="flash_point_C", T_K=np.nan, measured=45.0, target="flash")])
    comp = (np.array([0, 1]), np.array([0.5, 0.5]))
    out = s3.score_entity(ctx, comp, comp, t)
    assert len(calls) == 3                                           # M, R, R-leave-out: one function
    np.testing.assert_allclose(out.err_M, out.err_R)                 # identical composition -> identical error
    out2 = s3.score_entity(ctx, comp, (np.array([0]), np.array([1.0])), t)
    assert not np.allclose(out2.err_M, out2.err_R, equal_nan=True)   # a different composition changes it


def test_r_route_uses_only_whitelisted_bulk_rows():
    rows = []
    for prop, v in [("distillation_T0_C", 150.0), ("distillation_T10_C", 170.0), ("distillation_T50_C", 200.0), ("distillation_T90_C", 240.0),
                    ("distillation_T100_C", 265.0), ("aromatics_vol_pct", 18.0),
                    # forbidden for R: composition-derived or target-like quantities
                    ("distillation_simdis_T10_C", 1.0), ("pct_n_paraffin", 99.0), ("avg_mw_g_mol", 999.0), ("hydrogen_content_mass_pct", 14.0),
                    ("flash_point_C", 40.0), ("density_kg_m3", 800.0), ("avg_formula_C", 11.0)]:
        rows.append(dict(property=prop, value=v, method="D86" if "dist" in prop else "D1319"))
    d = s3.r_descriptors(pd.DataFrame(rows))
    assert set(d) <= set(s3.R_KEYS)
    assert d == {"d86_ibp": 150.0, "d86_t10": 170.0, "d86_t50": 200.0, "d86_t90": 240.0, "d86_fbp": 265.0, "aromatics_vol_pct": 18.0}
    assert not any(k in s3.R_PROPS for k in ("flash_point_C", "density_kg_m3", "freezing_point_C", "viscosity_kin_cSt"))
    src = inspect.getsource(s3.r_descriptors).lower()
    assert "composition" not in src and "gcxgc" not in src and "pct_" not in src


def test_bootstrap_is_seeded_and_signed():
    g = np.repeat(np.arange(10), 2)
    eM, eR = np.full(20, 1.0), np.full(20, 3.0)
    a = s3.paired_bootstrap(g, eM, eR)
    assert a == s3.paired_bootstrap(g, eM, eR)
    assert a[0] == 2.0 and a[1] == 2.0 and a[2] == 2.0
    assert s3.verdict_s31(5, 1.0, 0.5, 1.5, 0, 0) == "n too small"
