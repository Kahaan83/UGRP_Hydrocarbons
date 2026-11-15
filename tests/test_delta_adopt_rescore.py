"""delta_adopt_rescore (results/v5/BASE/delta_adopt_rescore/README.md): the BASE re-score under the adopted defaults (decisions 2026-11-15).

Adopted as the frozen v5 BASE on 2026-11-15 (user question U4, "Yea adopt"): the top-level results/v5/BASE files are the re-run files
(rerun/), and drift.csv documents every changed cell old (pre-adoption BASE, recoverable from git) -> new (adopted BASE).

Reads committed files only (adopted BASE, rerun/, drift.csv); nothing is recomputed and nothing is written.
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
# Every top-level file the --full re-run regenerates (rerun/ holds a copy of each; delta_files/ is not part of the BASE).
RERUN_FILES = sorted(p.name for p in RERUN.iterdir() if p.is_file())

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


def _same_prec(a, b) -> bool:
    """Equal up to CSV write precision.  The adopted files are the fixed point of the default `reproduce_base.py` (`--from-rows`) pass
    run by `--update-manifest` on the `--full` output (byte-identical to rerun/); re-scoring from the 8-significant-digit rows_t2_bulk.csv
    moves aggregates by at most 1.0e-5 absolute (density_error_vs_T.csv is written at 6 digits) / 1.2e-6 relative (U4, 2026-11-15)."""
    if _same(a, b):
        return True
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return False
    if np.isnan(fa) or np.isnan(fb):
        return False
    return abs(fa - fb) <= 1.1e-5 * max(1.0, abs(fa), abs(fb))


def _load(p):
    t = pd.read_csv(p)
    for k in SB_KEYS:
        t[k] = t[k].astype(str)
    return t.set_index(SB_KEYS)


def _cell_key(idx) -> str:
    return ";".join(f"{k}={v}" for k, v in zip(SB_KEYS, idx))


def test_adopted_base_equals_rerun():
    """U4: every file the re-run regenerated is, in the adopted BASE, the re-run file (cell-equal up to CSV precision, text-equal otherwise).
    Exactly 5 files carry precision-only differences (the --from-rows fixed point, see _same_prec); the other 15 are cell-equal at 1e-6."""
    assert len(RERUN_FILES) == 20
    prec_only = set()
    for name in RERUN_FILES:
        a, b = BASE / name, RERUN / name
        assert a.exists(), name
        if name.endswith(".csv"):
            ta, tb = pd.read_csv(a, low_memory=False), pd.read_csv(b, low_memory=False)
            assert list(ta.columns) == list(tb.columns) and ta.shape == tb.shape, name
            bad = [(i, c) for c in ta.columns for i in range(len(ta)) if not _same_prec(ta.at[i, c], tb.at[i, c])]
            assert not bad, (name, bad[:5])
            if any(not _same(ta.at[i, c], tb.at[i, c]) for c in ta.columns for i in range(len(ta))):
                prec_only.add(name)
        else:
            assert a.read_text(encoding="utf-8").replace("\r\n", "\n") == b.read_text(encoding="utf-8").replace("\r\n", "\n"), name
    assert prec_only <= {"cells_t2_bulk.csv", "density_error_vs_T.csv", "rows_t2_bulk.csv", "scoreboard.csv", "viscosity_error_vs_T.csv"}


def test_drift_csv_documents_old_to_new():
    """Every scoreboard cell listed in drift.csv carries the adopted value as `new` and a different pre-adoption value as `old`."""
    n = _load(BASE / "scoreboard.csv")
    d = pd.read_csv(D / "drift.csv", dtype={"old": str, "new": str})
    s = d[(d.path == "results/v5/BASE/scoreboard.csv") & (d.kind == "cell")]
    assert len(s) == 2021 and s.cell.is_unique
    have = {_cell_key(i): i for i in n.index}
    for r in s.itertuples():
        key, col = r.cell.rsplit("|", 1)
        assert key in have and col in n.columns, r.cell
        assert _same_prec(n.at[have[key], col], r.new), (r.cell, n.at[have[key], col], r.new)
        assert not _same(r.old, r.new), r.cell
    assert not d.cause.astype(str).str.startswith(("D1_", "D4_", "D5_", "D7_")).any()     # E1
    assert not d.cause.astype(str).str.startswith("unattributed").any()
    aj = s[s.cell.str.contains("property=flash_point_C;tier=T2_measured;set=ajftd_gcxgc;model=ideal_le_chatelier;variant=designated;")]
    b = aj[aj.cell.str.endswith("|bias")].iloc[0]
    assert float(b.old) == pytest.approx(-0.634, abs=1e-3) and float(b.new) == pytest.approx(-0.152, abs=1e-3)
    v = aj[aj.cell.str.endswith("|verdict")].iloc[0]
    assert (v.old, v.new) == ("no_skill_shown", "beats_trivial")


def test_freezing_and_carried_cells_not_in_drift():
    """E2: no drift.csv scoreboard cell is a freezing cell or a T1 / T3 / carried cell (no decision reaches them)."""
    n = _load(BASE / "scoreboard.csv")
    prop = n.index.get_level_values("property")
    cf = n.index.get_level_values("cells_file")
    untouched = {_cell_key(i) for i in n[(prop == "freezing_point_C") | ~cf.isin(["cells_t2_bulk.csv", "cells_t2_flash.csv"])].index}
    assert len(untouched) > 100
    d = pd.read_csv(D / "drift.csv")
    listed = {c.rsplit("|", 1)[0] for c in d[(d.path == "results/v5/BASE/scoreboard.csv") & (d.kind == "cell")].cell}
    assert not (listed & untouched)


def test_no_error_rows_in_adopted_base():
    """Declared check: no predict() exception was turned into an error: row (the only path for decisions 1, 4, 5 into a BASE cell)."""
    for root in (BASE, RERUN):
        b = pd.read_csv(root / "rows_t2_bulk.csv")
        f = pd.read_csv(root / "rows_t2_flash.csv")
        assert not b.api_method.astype(str).str.contains("error").any()
        assert not f.api_method.astype(str).str.contains("error").any()
