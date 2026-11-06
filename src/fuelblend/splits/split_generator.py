"""Deterministic split generator (seed 20261030). Writes JSON manifests with row ids.

S1 warm: all molecules seen (relative to S2/S3 holdout molecules).
S2 one molecule unseen; S3 all molecules unseen (molecule-grouped, strict).
real_fuel rows (no molecule lists) are excluded and handled by the locked-set logic.
"""
import json, os
import numpy as np
import pandas as pd
from .locked_sets import locked_sets

SEED = 20261030

def _mols(row):
    try:
        return json.loads(row) if isinstance(row, str) else (row or [])
    except Exception:
        return []

def split_s1_s2_s3(df, seed=SEED, hold_frac=0.1, exclude_locked=True):
    rng = np.random.default_rng(seed)
    locked_multi = locked_fuel = set()
    if exclude_locked:
        lk_path = os.path.join("splits", "locked.json")
        if os.path.exists(lk_path):
            lk = json.load(open(lk_path))
            locked_multi = set(lk["multi_system_ids_locked"])
            locked_fuel = set(lk["real_fuel_ids_locked"])
        else:
            raise FileNotFoundError("splits/locked.json missing; run locked_sets() first")
    all_mols = set()
    for c in df["comp_inchikeys"]:
        all_mols.update(_mols(c))
    mols = sorted(all_mols)
    rng.shuffle(mols)
    n = len(mols)
    n_hold = max(1, int(n * hold_frac))
    hold_s2 = set(mols[:n_hold])
    hold_s3 = set(mols[n_hold:2*n_hold]) if 2 * n_hold <= n else set(mols[n_hold:])
    idx_s1, idx_s2, idx_s3 = [], [], []
    for i in df.index:
        r = df.loc[i]
        if not _mols(r["comp_inchikeys"]) or r.get("kind") == "real_fuel":
            continue
        if r.get("system_id") in locked_multi:
            continue
        if str(r.get("fuel_id")) in locked_fuel:
            continue
        ms = set(_mols(r["comp_inchikeys"]))
        if r.get("kind") == "pure":
            # hold pure compounds out of S1 so their data never trains S1 models
            if ms <= hold_s3:
                idx_s3.append(i); continue
            if ms <= hold_s2:
                idx_s2.append(i); continue
        if ms <= hold_s3:
            idx_s3.append(i)
        elif ms & hold_s3:
            pass  # drop: mixes S3 molecules with training molecules
        elif ms & hold_s2:
            idx_s2.append(i)
        else:
            idx_s1.append(i)
    return {"S1": idx_s1, "S2": idx_s2, "S3": idx_s3,
            "holdout_molecules_s2": sorted(hold_s2),
            "holdout_molecules_s3": sorted(hold_s3)}

def leakage_check(df, res):
    s2_mols = set(res["holdout_molecules_s2"]); s3_mols = set(res["holdout_molecules_s3"])
    assert s2_mols and s3_mols, "holdout sets empty"
    assert not (s2_mols & s3_mols), "S2/S3 molecule overlap"
    for i in res["S3"]:
        ms = set(_mols(df.loc[i, "comp_inchikeys"]))
        assert ms <= s3_mols, "S3 row contains training molecule"
    for i in res["S2"]:
        ms = set(_mols(df.loc[i, "comp_inchikeys"]))
        assert (ms & s2_mols) or (ms & s3_mols), "S2 row fully seen"
        assert not ((ms & s2_mols) and (ms & s3_mols)), "S2 mixes S2 and S3 molecules"
    for i in res["S1"]:
        ms = set(_mols(df.loc[i, "comp_inchikeys"]))
        assert not (ms & s2_mols) and not (ms & s3_mols), "S1 contains held-out molecule"

def write_manifest(name, payload, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(payload, f, indent=1)

if __name__ == "__main__":
    df = pd.read_csv("data/processed/model_table.csv", low_memory=False)
    locked_sets()  # ensure splits/locked.json exists and is current before splitting
    res = split_s1_s2_s3(df)
    leakage_check(df, res)
    write_manifest("s1s2s3", res, "splits/s1s2s3.json")
    print({k: len(v) if isinstance(v, list) else len(v) for k, v in res.items()})
