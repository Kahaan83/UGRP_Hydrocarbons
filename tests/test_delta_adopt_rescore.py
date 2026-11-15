"""delta_adopt_rescore (results/v5/BASE/delta_adopt_rescore/README.md): the BASE re-score under the adopted defaults (decisions 2026-11-15).

Reads committed files only (frozen BASE, rerun/, drift.csv); nothing is recomputed and nothing is written.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "v5" / "BASE"
D = BASE / "delta_adopt_rescore"
RERUN = D / "rerun"
SB_KEYS = ["property", "tier", "set", "model", "variant", "label", "T_K", "cells_file", "status", "evidence_path"]

pytestmark = pytest.mark.skipif(not (D / "drift.csv").exists(), reason="delta_adopt_rescore results not written")


def _same(a, b) -> bool:
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        na = "" if a is None or (isinstance(a, float) and np.isnan(a)) else str(a)
        nb = "" if b is None or (isinstance(b, float) and np.isnan(b)) else str(b)
        return na == nb
    if fa == fb or (np.isnan(fa) and np.isnan(fb)):
        return True
    if np.isnan(fa) or np.isnan(fb):
        return False
    sc = max(abs(fa), abs(fb))
    return abs(fa - fb) <= 1e-6 * (sc if sc > 1 else 1.0)


def _load(p):
    t = pd.read_csv(p)
    for k in SB_KEYS:
        t[k] = t[k].astype(str)
    return t.set_index(SB_KEYS)


def test_scoreboard_drift_confined_to_drift_csv():
    """Every scoreboard cell that differs between the frozen BASE and the re-run is listed in drift.csv, and nothing else is."""
    o, n = _load(BASE / "scoreboard.csv"), _load(RERUN / "scoreboard.csv")
    assert set(o.index) == set(n.index)
    diff = set()
    for idx in o.index:
        for col in o.columns:
            if not _same(o.at[idx, col], n.at[idx, col]):
                diff.add(";".join(f"{k}={v}" for k, v in zip(SB_KEYS, idx)) + "|" + col)
    d = pd.read_csv(D / "drift.csv")
    listed = set(d[d.path == "results/v5/BASE/scoreboard.csv"].cell)
    assert diff == listed
    assert not d.cause.astype(str).str.startswith(("D1_", "D4_", "D5_", "D7_")).any()     # E1


def test_freezing_and_carried_cells_bit_identical():
    """E2: no freezing cell and no T1 / T3 / carried cell moves (no decision reaches them)."""
    o, n = _load(BASE / "scoreboard.csv"), _load(RERUN / "scoreboard.csv")
    prop = o.index.get_level_values("property")
    cf = o.index.get_level_values("cells_file")
    oo = o[(prop == "freezing_point_C") | ~cf.isin(["cells_t2_bulk.csv", "cells_t2_flash.csv"])]
    assert len(oo) > 100
    bad = [(i, c) for i in oo.index for c in oo.columns if not _same(oo.at[i, c], n.at[i, c])]
    assert not bad, bad[:5]


def test_no_error_rows_in_rerun():
    """Declared check: no predict() exception was turned into an error: row (the only path for decisions 1, 4, 5 into a BASE cell)."""
    b = pd.read_csv(RERUN / "rows_t2_bulk.csv")
    f = pd.read_csv(RERUN / "rows_t2_flash.csv")
    assert not b.api_method.astype(str).str.contains("error").any()
    assert not f.api_method.astype(str).str.contains("error").any()
