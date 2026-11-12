"""Fast tests for scripts/v5/s4b_heat_measured.py (no S4 data needed)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("s4b_mod", ROOT / "scripts" / "v5" / "s4b_heat_measured.py")
s4b = importlib.util.module_from_spec(_spec)
sys.modules["s4b_mod"] = s4b
_spec.loader.exec_module(s4b)


def test_method_class_assignment():
    assert s4b.classify_method("D3338") == "CALCULATED"
    assert s4b.classify_method("D 3338") == "CALCULATED"
    assert s4b.classify_method("D 4809") == "MEASURED"
    assert s4b.classify_method("D4809 (SwRI[32])") == "MEASURED"
    assert s4b.classify_method("D240") == "MEASURED"
    for code in ("400C", "400D", "231A", "231B", "400IPB", float("nan"), None, "XYZ"):
        assert s4b.classify_method(code) == "UNCLASSIFIED"


def test_labels_use_only_rows_of_the_class():
    lg = pd.DataFrame({"fuel_id": [1, 1, 2, 3], "property": [s4b.HEAT] * 4,
                       "method": ["D3338", "D4809", "D3338", "400C"], "value": [43.0, 43.5, 43.1, 43.2]})
    lab = s4b.build_heat_labels(lg, {}).set_index("fuel_id")
    assert lab.loc[1, "y_meas"] == 43.5 and lab.loc[1, "y_calc"] == 43.0
    assert np.isnan(lab.loc[2, "y_meas"]) and lab.loc[2, "y_calc"] == 43.1
    assert np.isnan(lab.loc[3, "y_meas"]) and np.isnan(lab.loc[3, "y_calc"])
    assert lab.loc[1, "fuel_class"] == "CALCULATED+MEASURED"


def test_within_constant_excludes_scored_group_and_matches_class():
    df = pd.DataFrame({"cell": ["a", "a", "b", "b", "c"], "cls": ["M", "C", "M", "C", "M"],
                       "y": [100.0, 1.0, 10.0, 2.0, 20.0]})
    p, nfb = s4b.within_const_logo(df, "y", "cls")
    assert nfb == 0
    # row 0 (cell a, class M): mean of M rows in b and c = 15, the scored value 100 must not enter
    assert p[0] == 15.0
    # row 1 (cell a, class C): only class C training row is b -> 2.0
    assert p[1] == 2.0
    # row 4 (cell c, class M): M rows in a and b = (100+10)/2
    assert p[4] == 55.0


def test_bootstrap_sign_positive_when_model_better():
    g = np.repeat(np.arange(10), 5)
    ref = np.full(50, 1.0)
    good = np.full(50, 0.5)
    d, lo, hi, ng = s4b.boot_pos(ref, good, g, B=500)
    assert d > 0 and lo > 0 and ng == 10
    d2, lo2, hi2, _ = s4b.boot_pos(good, ref, g, B=500)
    assert d2 < 0 and hi2 < 0


def test_fewer_than_8_groups_is_not_evaluable():
    n = 30
    P = pd.DataFrame({"fuel_id": np.tile(np.arange(n), 2), "cell": np.tile(np.arange(n) % 3, 2).astype(str),
                      "source_group": "x", "y": np.tile(np.ones(n), 2),
                      "pred": np.r_[np.full(n, 1.5), np.full(n, 1.1)],
                      "model": ["C_within"] * n + ["ridge"] * n})
    rows = s4b.score(P, "t")
    r = [x for x in rows if x["model"] == "ridge"][0]
    assert r["evaluable"] is False or not r["evaluable"]
    assert not r["skill"] and "not evaluable" in r["note"]
