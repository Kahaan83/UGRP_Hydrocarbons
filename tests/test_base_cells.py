"""v5 BASE Part B: scoreboard cells.  One test function per cells_* script (Tasks 4-8)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "v5" / "BASE"

SCHEMA = ("property tier set model variant status n_groups n_rows bias mae rmse unsafe_rate trivial_model mae_trivial "
          "diff_trivial_minus_model ci_lo ci_hi verdict margin90 coverage90 n_excluded_censored n_excluded_nonfinite "
          "n_low_coverage evidence_path label").split()


def _schema_ok(c):
    assert not [k for k in SCHEMA if k not in c.columns]
    assert c.status.isin(["recomputed_v5", "aggregated_v5", "carried_v4", "no_evaluation"]).all()
    assert c.label.isin(["descriptive", "dev", "fold_b", "already_read", "exploratory"]).all()
    assert c.verdict.isin(["beats_trivial", "no_skill_shown", "worse_than_trivial", "n_too_small", "not_applicable"]).all()
    assert (c[c.n_groups < 8].verdict.isin(["n_too_small", "not_applicable"])).all()


def test_cells_t1_lab():
    c = pd.read_csv(BASE / "cells_t1_lab.csv")
    _schema_ok(c)
    v = c[(c.property == "viscosity_kin_cSt") & (c.set == "T_lt_263K") & (c.label == "dev")].iloc[0]
    assert v.n_groups == 20 and v.bias == pytest.approx(-0.1135, abs=1e-3) and v.margin90 == pytest.approx(0.4154, abs=1e-3)
    f = c[(c.property == "flash_point_C") & (c.variant == "designated")].iloc[0]
    assert f.n_groups == 31 and f.n_rows == 126 and f.bias == pytest.approx(-3.30, abs=0.01)
    s = c[(c.property == "freezing_point_C") & (c.variant == "designated")].iloc[0]
    assert np.isinf(s.margin90) and s.n_groups == 8
    d = c[(c.property == "density_kg_m3") & (c.set == "thermoml2020_new_systems")]
    # plan expected >1000 rows; observed: 1239 of 1599 non-pure rows are systems already in measurements_dedup, 63 rows have an invalid pure
    # property -> 297 scored rows in 4 systems (see the cell note)
    assert len(d) == 1 and d.n_rows.iloc[0] == 297 and d.n_groups.iloc[0] == 4 and d.label.eq("dev").all()


def test_cells_t2_flash():
    r = pd.read_csv(BASE / "rows_t2_flash.csv", dtype={"fuel_id": str})
    g1 = pd.read_csv(ROOT / "results/v5/G1/predictions.csv", dtype={"fuel_id": str})
    g1 = g1[(g1.property == "flash_point_C") & (g1.model == "ideal")]
    j = r[r.source == "ajftd_gcxgc"].merge(g1, on="fuel_id")
    assert len(j) == 37 and np.allclose(j.pred_designated, j.pred, atol=0.05)          # R-BASE-2
    c = pd.read_csv(BASE / "cells_t2_flash.csv")
    _schema_ok(c)
    assert set(c.variant) == {"designated", "lfl_overlay", "mesitylene_fixed", "lfl_overlay+mesitylene_fixed"}
    assert (c[c.n_groups < 8].verdict == "n_too_small").all()
    assert c.unsafe_rate.between(0, 1).all()
    assert (c.model == "api_as_shipped").any()


def test_cells_t2_bulk():
    vt = pd.read_csv(BASE / "viscosity_error_vs_T.csv")
    a = vt[vt.source == "av3322"].set_index("T_K")
    assert len(a) == 9
    tol = 0.02                                                        # 0.05 if route_vs_s3_differences.csv existed (it does not)
    assert a.loc[233.15, "bias_ln"] == pytest.approx(-0.637, abs=tol)      # C1 R-C1-5
    assert a.loc[313.15, "bias_ln"] == pytest.approx(0.066, abs=tol)
    c = pd.read_csv(BASE / "cells_t2_bulk.csv")
    _schema_ok(c)
    k1 = c[(c.property == "viscosity_kin_cSt") & (c.variant == "k1_cold_offset_loso") & (c.set == "ajftd_gcxgc") & (c.T_K == 253.15)]
    assert abs(k1.bias.iloc[0]) <= 0.10                                     # C1 R-C1-1
    fz = c[(c.property == "freezing_point_C") & (c.set == "pooled") & (c.label == "descriptive")]
    assert fz[fz.variant == "f1_recalibrated"].mae.iloc[0] == pytest.approx(3.46, abs=0.02)
    for p in ["heat_of_combustion_net_MJ_kg", "smoke_point_mm", "distillation_T10_C", "density_kg_m3"]:
        assert (c.property == p).any()
    r = pd.read_csv(BASE / "rows_t2_bulk.csv")
    assert np.isfinite(r.loc[r.status == "ok", "pred"]).all()


def test_cells_t3_bulk():
    c = pd.read_csv(BASE / "cells_t3_bulk.csv")
    _schema_ok(c)
    f = c[(c.property == "flash_point_C") & (c.model == "const") & (c.label == "dev") & (c.variant == "designated")].iloc[0]
    assert f.n_rows == 1728 and f.mae == pytest.approx(3.943, abs=1e-3)
    assert (c[c.label == "already_read"].status == "aggregated_v5").all()
    assert ((c.property == "flash_point_C") & (c.trivial_model == "logo_constant_within_method")).any()
    p = c[c.property == "pour_point_C"]
    assert set(p.model) >= {"const", "ridge", "freezing_offset"} and (p.n_rows == 57).all()


def test_cells_carried():
    c = pd.read_csv(BASE / "cells_carried.csv")
    _schema_ok(c)
    fl = c[(c.property == "flash_point_C") & (c.set == "crc647_r3") & (c.label == "dev") & (c.variant == "designated")].iloc[0]
    assert fl.n_groups == 45 and fl.bias == pytest.approx(0.44, abs=0.02) and fl.unsafe_rate == pytest.approx(0.133, abs=0.005)
    assert c[c.set == "crc647_r3"].label.isin(["dev", "fold_b"]).all()
    assert (c[c.set == "face_diesel"].verdict == "n_too_small").all()
    lf = c[(c.property == "flash_point_C") & (c.set == "crc647_r3") & (c.label == "dev") & (c.variant == "lfl_overlay")].iloc[0]
    assert lf.bias == pytest.approx(-0.726, abs=0.01)                       # L1 comparison.csv, after
    v = c[(c.property == "viscosity_kin_cSt") & (c.set == "crc647_r3") & (c.label == "dev") & (c.model == "Recalibrated Refutas") & c.T_K.isna()]
    assert len(v) == 1 and v.margin90.iloc[0] == pytest.approx(0.3465, abs=1e-3)
    assert len(c[(c.property == "viscosity_kin_cSt") & (c.set == "crc647_r3") & (c.label == "dev")].model.unique()) == 8
    d = c[c.property == "designs_certified"].iloc[0]
    assert d.n_rows == 4 and d.bias == 0
    assert (c[c.status == "no_evaluation"].reason.str.len() > 0).all()
