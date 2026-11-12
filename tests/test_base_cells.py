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
