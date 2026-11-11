"""S4 split determinism, group-leakage and run-once guards (scripts/v5/eval_bulk_descriptor.py)."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("eval_bulk_descriptor", ROOT / "scripts" / "v5" / "eval_bulk_descriptor.py")
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)


@pytest.fixture(scope="module")
def built():
    t, dup = E.build_table()
    return t, dup, E.make_splits(t, dup)


def test_split_deterministic(built):
    t, dup, sp = built
    t2, dup2 = E.build_table()
    assert json.dumps(sp, sort_keys=True) == json.dumps(E.make_splits(t2, dup2), sort_keys=True)


def test_committed_splits_match_rebuild(built):
    p = E.OUT / "splits_ajftd_v5.json"
    if not p.exists():
        pytest.skip("splits not written yet")
    assert p.read_text() == json.dumps(built[2], indent=1, sort_keys=True)


def test_no_group_leakage_between_dev_and_locked(built):
    t, _, sp = built
    role = {int(k): v for k, v in sp["fuel_role"].items()}
    t = t.assign(role=t["fuel_id"].map(role))
    assert set(t["role"]) <= {"dev", "locked", "excluded_leak"}
    dev_cells = set(t.loc[t.role == "dev", "cell"])
    locked_cells = set(t.loc[t.role == "locked", "cell"])
    assert dev_cells and locked_cells and not (dev_cells & locked_cells)
    # no raw source/year cell is split, and no identical numeric vector sits on both sides
    cols = E.FEATURES + [E.TARGET_COL[k] for k in E.TARGETS] + ["visc_m20_assumed"]
    keys = {}
    for _, r in t.iterrows():
        if r[cols].notna().sum() >= 6:
            keys.setdefault(E._vec_key(r, cols), set()).add(r["role"])
    assert not [k for k, v in keys.items() if {"dev", "locked"} <= v]
    assert abs(sp["n_locked_cells"] - np.ceil(E.LOCK_FRACTION * sp["n_cells"])) < 1e-9


def test_dedup_one_row_per_fuel(built):
    t, dup, _ = built
    assert t["fuel_id"].is_unique
    assert not set(dup) & set(t["fuel_id"].astype(int))  # merged ids are not kept as rows


def test_splits_refuse_overwrite(tmp_path, built):
    sp = built[2]
    p = tmp_path / "s.json"
    assert E.write_splits_guarded(sp, p) == "written"
    assert E.write_splits_guarded(sp, p) == "unchanged"
    other = dict(sp, seed="different")
    with pytest.raises(FileExistsError):
        E.write_splits_guarded(other, p)
    assert json.loads(p.read_text())["seed"] == sp["seed"]


def test_locked_refuses_second_run(tmp_path):
    (tmp_path / "locked_scores.csv").write_text("x\n")
    with pytest.raises(FileExistsError):
        E.cmd_locked(tmp_path)
    (tmp_path / "locked_scores.csv").unlink()
    (tmp_path / "locked_run.json").write_text("{}")
    with pytest.raises(FileExistsError):
        E.cmd_locked(tmp_path)


def test_bootstrap_verdict_logic():
    rng = np.random.default_rng(1)
    g = np.repeat(np.arange(10), 20)
    ec = rng.uniform(1, 2, 200)
    d, lo, hi, n = E.boot_diff(ec - 0.5, ec, g)
    assert n == 10 and hi < 0 and E.verdict("ridge", lo, hi, n) == "beats const"
    d, lo, hi, n = E.boot_diff(ec, ec, g)
    assert lo == hi == 0 and E.verdict("ridge", lo, hi, n) == "no skill shown"
    assert E.verdict("ridge", -1, -0.5, 3).startswith("not decidable")


def test_flash_correlation_sanity():
    import pandas as pd
    X = pd.DataFrame({"distillation_T10_C": [175.0, np.nan]})
    m = E.Model("api_2b7_t10", "flash").fit(pd.DataFrame({"distillation_T10_C": [1.0]}), [40.0])
    p = m.predict(X)
    assert 45 < p[0] < 60 and p[1] == 40.0
    b = E.Model("butler_t10", "flash").fit(pd.DataFrame({"distillation_T10_C": [1.0]}), [40.0]).predict(X)
    assert 45 < b[0] < 55
