"""v5 BASE Task 11: cold-viscosity margin delta."""
import importlib.util
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "v5" / "BASE" / "delta_cold_viscosity"
_spec = importlib.util.spec_from_file_location("delta_cv", ROOT / "scripts/v5/base/delta_cold_viscosity.py")
dcv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dcv)


def test_offset_never_uses_scored_group_or_source():
    rows = pd.read_csv(ROOT / "results/v5/BASE/rows_t2_bulk.csv")
    cf = dcv.crossfit(rows)
    assert len(cf) > 100
    assert (cf.train_source != cf.source).all()
    assert cf.apply(lambda r: r.group not in r.train_groups, axis=1).all()


def test_crossfit_matches_frozen_k1_rows():
    rows = pd.read_csv(ROOT / "results/v5/BASE/rows_t2_bulk.csv")
    cf = dcv.crossfit(rows)
    k = rows[(rows.variant == "k1_cold_offset_loso") & (rows.status == "ok") & (rows.property == "viscosity_kin_cSt")]
    m = cf.merge(k[["entity", "T_K", "pred"]], on=["entity", "T_K"], suffixes=("", "_frozen"))
    assert len(m) > 100
    assert np.allclose(m.pred_k1, m.pred_frozen, rtol=1e-6)


def test_margin_is_inf_below_nine_groups():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert dcv.margin(np.full(8, 0.1), np.arange(8)) == float("inf")
        assert np.isfinite(dcv.margin(np.linspace(0.05, 0.5, 9), np.arange(9)))


def test_outputs_exist_and_cold_table_complete():
    subprocess.run([sys.executable, str(ROOT / "scripts/v5/base/delta_cold_viscosity.py")], cwd=ROOT, check=True)
    t = pd.read_csv(D / "pure_cold_extrapolation.csv")
    assert {"smiles", "mean_mass_pct", "data_T_min_K", "source_at_253K", "source_at_233K", "err_vs_boehm_253K", "err_vs_boehm_233K"} <= set(t.columns)
    assert (t.mean_mass_pct > 2).all() and set(t.source_at_253K) <= {"fit", "ml_extrapolation"}
    assert (D / "proposed_margins.yaml").exists() and (D / "delta.csv").exists()
    assert not (ROOT / "configs" / "proposed_margins.yaml").exists()
