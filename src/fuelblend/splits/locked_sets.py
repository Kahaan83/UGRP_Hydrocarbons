"""Locked-set selection (WP2). NEVER used for tuning; access must be logged.

- real fuels: 20% of distinct fuel_id per source (stratified), seed 20261030.
- multi-component records (kind='multi'): split 20% to a dev-lock stratum too.
Outputs splits/locked.json + splits/locked_final.json (dev-lock vs final-lock).
"""
import json, os
import numpy as np
import pandas as pd

SEED = 20261030
LOCK_FRAC = 0.20


def select_locked(df):
    rng = np.random.default_rng(SEED)
    lock_ids = []
    # real fuels: stratify by source, 20% of fuel_ids each
    rf = df[df["kind"] == "real_fuel"]
    for src in rf["sources"].unique():
        ids = sorted(rf[rf["sources"] == src]["fuel_id"].astype(str).unique())
        ids = list(map(str, ids))
        k = max(1, int(np.ceil(len(ids) * LOCK_FRAC)))
        lock_ids += [str(x) for x in rng.choice(len(ids), size=k, replace=False)]
    locked_fuel = sorted(set(ids[i] for i in map(int, [])))  # placeholder replaced below
    return lock_ids


def locked_sets(path="data/unified/measurements_dedup.csv", out_dir="splits"):
    df = pd.read_csv(path, low_memory=False)
    rng = np.random.default_rng(SEED)
    os.makedirs(out_dir, exist_ok=True)

    # real fuels, 20% per source
    rf = df[df["kind"] == "real_fuel"]
    locked_fuel = []
    for src in rf["sources"].unique():
        grp = rf[rf["sources"] == src]
        ids = sorted(set(map(str, grp["fuel_id"].astype(str))))
        k = max(1, int(np.ceil(len(ids) * LOCK_FRAC)))
        pick = rng.choice(len(ids), size=k, replace=False)
        locked_fuel += [str(ids[i]) for i in pick]

    # multi-component systems, 20% per source bucket
    locked_multi = []
    for src in sorted(df[df["kind"] == "multi"]["sources"].unique()):
        grp = df[(df["kind"] == "multi") & (df["sources"] == src)]
        sids = sorted(set(map(str, grp["system_id"])))
        k = max(1, int(np.ceil(len(sids) * LOCK_FRAC)))
        pick = rng.choice(len(sids), size=k, replace=False)
        locked_multi += [str(sids[i]) for i in pick]

    # dev-lock vs final-lock (half/half within each pool)
    fids = sorted(set(locked_fuel)); mids = sorted(set(locked_multi))
    dev_f = set(str(fids[i]) for i in rng.choice(len(fids), size=max(1, len(fids) // 2), replace=False))
    dev_m = set(str(mids[i]) for i in rng.choice(len(mids), size=max(1, len(mids) // 2), replace=False))

    payload = {
        "seed": SEED,
        "real_fuel_ids_locked": sorted(fids),
        "multi_system_ids_locked": sorted(mids),
        "dev_lock": {"real_fuel_ids": sorted(dev_f), "multi_system_ids": sorted(dev_m)},
        "rule": "Locked rows never used for tuning/model selection; access logged in reports/LOCKED_ACCESS.log",
    }
    with open(os.path.join(out_dir, "locked.json"), "w") as f:
        json.dump(payload, f, indent=1)
    return payload


if __name__ == "__main__":
    p = locked_sets()
    print("locked real fuels:", len(p["real_fuel_ids_locked"]),
          "locked multi systems:", len(p["multi_system_ids_locked"]))
