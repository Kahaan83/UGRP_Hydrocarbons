import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from fuelblend.v7.guard import LOCKS, LockedRowError, assert_dev_only

ROOT = Path(__file__).resolve().parents[2]


def test_earlier_lock_id_raises():
    from fuelblend.v7.guard import locked_ids
    some = sorted(locked_ids("splits_locked"))[:1]
    with pytest.raises(LockedRowError):
        assert_dev_only(some, locks=("splits_locked",))


def test_model_table_default_is_base_path():
    code = "import fuelblend.data.model_table as m; print(m.DEDUP)"
    env = {k: v for k, v in os.environ.items() if k != "FUELBLEND_DEDUP_PATH"}; env["PYTHONPATH"] = str(ROOT / "src")
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True).stdout
    assert out.strip() == "data/unified/measurements_dedup.csv"


def test_env_var_switches_path():
    code = "import fuelblend.data.model_table as m; print(m.DEDUP)"
    env = dict(os.environ, FUELBLEND_DEDUP_PATH="data/unified_v9/measurements_dedup.csv", PYTHONPATH=str(ROOT / "src"))
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True).stdout
    assert out.strip() == "data/unified_v9/measurements_dedup.csv"


def test_v9_lock_disjoint_from_dev_and_earlier_locks():
    lock = set(json.loads((ROOT / "results/v9/lock.json").read_text())["ids"])
    assert lock
    assert_dev_only(sorted(lock), locks=tuple(LOCKS))  # new ids do not collide with earlier locks
