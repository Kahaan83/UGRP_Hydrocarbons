import json

import numpy as np
import pandas as pd

from fuelblend.v7 import guard, runner


def test_folds_match_v6():
    from sklearn.model_selection import GroupKFold
    comp, folds, tasks = runner.load_track_b_pool()
    pool = pd.read_csv(runner.POOL_CSV)
    c = pool.drop_duplicates("inchikey")[["inchikey", "smiles", "is_hc"]].reset_index(drop=True)
    c["grp"] = c.smiles.map(runner._group)
    ref = np.full(len(c), -1)
    for k, (_, te) in enumerate(GroupKFold(5).split(c, groups=c.grp)):
        ref[te] = k
    assert np.array_equal(folds, ref) and len(comp) == len(c) and "ysi" in tasks


def test_ysi_locked_label_masked():
    comp, folds, tasks = runner.load_track_b_pool()
    locked = comp.inchikey.isin(guard.locked_ids("ysi_holdout"))
    assert locked.sum() == comp.attrs["n_locked_ysi"] > 0
    assert comp.loc[locked, "ysi"].isna().all()
    other = [t for t in tasks if t != "ysi"]
    assert comp.loc[locked, other].notna().any().any()


def test_write_results_meta_keys(tmp_path):
    runner.write_results(tmp_path, summary=pd.DataFrame({"a": [1]}), oof=pd.DataFrame({"b": [2]}), meta={"x": 1},
                         drops={"ysi_holdout": 3})
    m = json.loads((tmp_path / "meta.json").read_text())
    assert {"device", "torch_version", "cuda_version", "git_commit", "guard_drops", "x"} <= set(m)
    assert (tmp_path / "cv_summary.csv").exists() and (tmp_path / "oof.csv").exists()


def test_paired_boot_gain_sign():
    g = np.repeat(np.arange(20), 3)
    lo, hi = runner.paired_boot_gain(np.full(60, 1.0), np.full(60, 2.0), g)
    assert lo > 0 and hi > 0
