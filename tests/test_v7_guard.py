import json

import numpy as np
import pandas as pd
import pytest

from fuelblend.v7 import guard
from fuelblend.v7.guard import LOCKS, LockedRowError, ROOT


def _one(name):
    return sorted(guard.locked_ids(name))[0]


@pytest.mark.parametrize("name", sorted(LOCKS))
def test_guard_refuses_each_lock_file(name):
    with pytest.raises(LockedRowError, match=name):
        guard.assert_dev_only([_one(name)], locks=(name,))


def test_guard_passes_dev_ids():
    ys = json.loads((ROOT / "results/v6/gnn/ysi_holdout/lock.json").read_text())
    pool = pd.read_csv(ROOT / "results/v6/gnn/pool.csv", usecols=["inchikey"])
    dev = pool.loc[~pool.inchikey.isin(ys["locked_inchikeys"]), "inchikey"].iloc[0]
    guard.assert_dev_only([dev], locks=("ysi_holdout",))
    guard.assert_dev_only(["zzz_not_an_id"], locks=tuple(LOCKS))


def test_mask_label_only():
    ik = _one("ysi_holdout")
    df = pd.DataFrame({"inchikey": [ik, "DEV-X"], "ysi": [10.0, 20.0], "flash_C": [50.0, 60.0]})
    out, n = guard.mask_locked(df, "inchikey", locks=("ysi_holdout",), label_cols=["ysi"])
    assert n == 1 and len(out) == 2
    assert np.isnan(out.loc[0, "ysi"]) and out.loc[0, "flash_C"] == 50.0
    assert out.loc[1, "ysi"] == 20.0


def test_drop_counts():
    ids = sorted(guard.locked_ids("flash_oxy_lock"))[:3]
    df = pd.DataFrame({"system_id": ids + ["dev1", "dev2"], "y": range(5)})
    out, n = guard.mask_locked(df, "system_id", locks=("flash_oxy_lock",))
    assert n == 3 and len(out) == 2


def test_prefix_canonical_and_explicit_locks():
    ids = [f"leenhouts:{i}" for i in sorted(guard.locked_ids("mixture_lock")) if ":" not in i][:1]
    with pytest.raises(LockedRowError):
        guard.assert_dev_only(ids, locks=("mixture_lock",))
    with pytest.raises(ValueError):
        guard.assert_dev_only(ids, locks=())
