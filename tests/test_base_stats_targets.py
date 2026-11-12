"""v5 BASE Task 3: shared group statistics and the unified T2 target table."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.eval.group_stats import (cell_metrics, classify_method, flash_family, logo_constant,
                                        paired_group_bootstrap, verdict)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "v5" / "BASE"


def test_bootstrap_matches_s5_committed_row():
    s5 = pd.read_csv(ROOT / "results/v5/S5/s5_predictions.csv")
    h = s5[s5.role == "headline"].dropna(subset=["C", "F1", "freezing_C"])      # S5 scores the 49 groups with a finite F1
    assert len(h) == 49
    abs_c, abs_f1 = (h.C - h.freezing_C).abs(), (h.F1 - h.freezing_C).abs()
    d, lo, hi = paired_group_bootstrap(abs_c, abs_f1, h.group)
    assert d == pytest.approx(2.47, abs=0.01)
    assert lo > 0


def test_logo_constant_never_sees_own_group():
    y = np.array([1.0, 1.0, 9.0]); g = np.array(["a", "a", "b"])
    assert np.allclose(logo_constant(y, g), [9.0, 9.0, 1.0])


def test_logo_constant_within_strata():
    y = np.array([1.0, 3.0, 10.0, 20.0]); g = np.array(["a", "b", "c", "d"]); s = np.array(["x", "x", "y", "y"])
    assert np.allclose(logo_constant(y, g, s), [3.0, 1.0, 20.0, 10.0])


def test_verdict_small_n():
    assert verdict(7, 0.5, 1.0) == "n_too_small"
    assert verdict(8, 0.5, 1.0) == "beats_trivial"
    assert verdict(8, -1.0, -0.2) == "worse_than_trivial"
    assert verdict(8, -0.5, 0.4) == "no_skill_shown"


def test_cell_metrics_unsafe_rate():
    err = np.array([6.0, -2.0, 1.0, 7.0])
    m = cell_metrics(err, np.array(list("abcd")), unsafe=lambda e: e > 5)
    assert m["n_groups"] == 4 and m["n_rows"] == 4
    assert m["unsafe_rate"] == pytest.approx(0.5) and m["bias"] == pytest.approx(3.0)
    assert m["mae"] == pytest.approx(4.0)


def test_method_families_moved_unchanged():
    assert flash_family("D 93") == "D93_or_IP170" and flash_family("D56") == "D56"
    assert flash_family("2209") == "spec_table_220x" and flash_family("zzz") == "other"
    assert classify_method("D3338") == "CALCULATED" and classify_method("D 4809") == "MEASURED"


# ---------------------------------------------------------------- target table (built by scripts/v5/base/build_t2_targets.py)
def _t():
    return pd.read_csv(BASE / "t2_targets.csv")


def test_targets_exclude_censored_and_count_them():
    t = _t()
    assert t[t.property == "freezing_point_C"].qualifier_dropped.sum() >= 1
    cens = pd.read_csv(BASE / "t2_censored_excluded.csv")
    assert len(cens) >= 5 and cens.qualifier.notna().all()
    # a fuel whose only freezing values are censored has no freezing target row
    only = set(cens[cens.property == "freezing_point_C"].entity) - set(t[t.property == "freezing_point_C"].entity)
    assert only


def test_av3322_has_nine_temperatures():
    t = _t()
    a = t[(t.source == "av3322") & (t.property == "viscosity_kin_cSt")]
    assert a.groupby("entity").T_K.nunique().min() == 9


def test_one_median_target_per_key_with_method_family():
    t = _t()
    assert not t.duplicated(["entity", "property", "T_K"]).any()
    assert t[t.property.isin(["flash_point_C", "heat_of_combustion_net_MJ_kg"])].method_family.notna().all()


def test_one_chosen_fuel_per_group():
    f = pd.read_csv(BASE / "t2_fuels.csv")
    assert (f.groupby("group").chosen_for_group.sum() == 1).all()


def test_low_coverage_flag_uses_090():
    f = pd.read_csv(BASE / "t2_fuels.csv")
    assert (f.low_coverage == (f.unmapped_fraction > 0.10)).all()
