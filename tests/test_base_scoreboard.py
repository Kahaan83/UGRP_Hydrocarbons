"""v5 BASE Task 9: scoreboard assembly, required cells, data-gap list, drift check."""
import importlib.util
import shutil
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "v5" / "BASE"


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


bs = _load("build_scoreboard_under_test", "scripts/v5/base/build_scoreboard.py")
ra = _load("reproduce_all_under_test", "scripts/v4/reproduce_all.py")
REQUIRED_CELLS = bs.REQUIRED_CELLS


def test_every_required_cell_exists():                      # Review Focus 4, R-BASE-1
    s = pd.read_csv(BASE / "scoreboard.csv")
    have = set(zip(s.property.map(bs.norm_prop), s.tier))
    assert not [c for c in REQUIRED_CELLS if (bs.norm_prop(c[0]), c[1]) not in have]


def test_no_evaluation_rows_have_reasons():
    g = pd.read_csv(BASE / "data_gaps.csv")
    assert len(g) > 0
    assert g.reason.notna().all() and g.what_would_fill_it.notna().all()
    s = pd.read_csv(BASE / "scoreboard.csv")
    assert len(g) == int(((s.status == "no_evaluation") | (s.verdict == "n_too_small")).sum())


def test_small_n_never_gets_skill_verdict():                # Review Focus 3
    s = pd.read_csv(BASE / "scoreboard.csv")
    assert s[(s.n_groups < 8) & s.verdict.isin(["beats_trivial", "worse_than_trivial", "no_skill_shown"])].empty


def test_designated_columns_untouched_by_variants():        # Review Focus 5
    s = pd.read_csv(BASE / "scoreboard.csv")
    d = s[(s.property == "flash_point_C") & (s.tier == "T2_measured") & (s.set == "ajftd_gcxgc") & (s.model == "ideal_le_chatelier") & (s.variant == "designated")].iloc[0]
    # headline = one fuel per double-count group (31); G1's 37-row number is reproduced from rows_t2_flash.csv in R-BASE-2
    assert d.n_groups == 31 and d.bias == pytest.approx(-0.634, abs=0.01) and d.mae == pytest.approx(3.381, abs=0.01)
    v = s[(s.property == "flash_point_C") & (s.tier == "T2_measured") & (s.set == "ajftd_gcxgc") & (s.model == "ideal_le_chatelier") & (s.variant == "lfl_overlay")].iloc[0]
    assert v.bias != d.bias                                  # the variant is its own row; the designated row did not move
    r2 = bs.r_base_2()
    g1 = r2[r2.cell.str.contains("all 37")]
    assert g1.within.all() and len(g1) == 2                  # G1's -0.426 / 3.536 reproduced on its own 37 rows


def test_tampered_file_makes_reproduce_exit_1(tmp_path):    # same pattern as tests/test_reproduce_drift.py
    for f in ("scoreboard.csv", "data_gaps.csv"):
        shutil.copy(BASE / f, tmp_path / f)
    paths = [tmp_path / "scoreboard.csv", tmp_path / "data_gaps.csv"]
    man = tmp_path / "SUMS.txt"
    ra.write_manifest(man, root=tmp_path, paths=paths)
    assert ra.check_manifest(man, root=tmp_path, paths=paths) is False
    (tmp_path / "scoreboard.csv").write_text((tmp_path / "scoreboard.csv").read_text(encoding="utf-8") + "x", encoding="utf-8")
    assert ra.check_manifest(man, root=tmp_path, paths=paths) is True
