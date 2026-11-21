"""v6 Track A step 1: freezing-point dataset and holdout lock.

Declared before any v6 model is run.  Rules (also written to results/v6/freezing/holdout_lock.json):
 * Start from the v5 T2 freezing entities (results/v5/BASE/delta_freeze_t2_unifac/rows.csv, 58 rows) and add AJFTD fuels
   20, 21, 22, 64, 154, 155 (32 is the same fuel as FuelLib posf10264 and is not added).
 * CV/holdout group = v6 dup_group (data/real_fuels/v6/ajftd_fuel_groups.csv); one entity per group (v5 source priority).
 * Holdout (locked) = v5 S4 locked fuels (role=='locked' in v5 rows) U fuels in the older splits/locked.json
   (real_fuel_ids_locked, which v5 headline did not honour) U AJFTD 20 (also in splits/locked.json).
 * Target for new AJFTD fuels = median of uncensored readings (same rule as s5_freezing_nalkane.freezing_ajftd).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/v6/freezing"
V5ROWS = ROOT / "results/v5/BASE/delta_freeze_t2_unifac/rows.csv"
GROUPS = ROOT / "data/real_fuels/v6/ajftd_fuel_groups.csv"
LONG6 = ROOT / "data/real_fuels/v6/ajftd_fuels_long.csv"
LOCKED_OLD = ROOT / "splits/locked.json"
NEW_IDS = ["20", "21", "22", "64", "154", "155"]
PRIO = {"AV-33-22": 0, "FuelLib": 1, "AJFTD": 2}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    v5 = pd.read_csv(V5ROWS)
    g = pd.read_csv(GROUPS, dtype=str)
    gid = dict(zip(g.fuel_id, g.dup_group))
    posf_groups = set(g.dup_group)
    old_locked = set(json.load(open(LOCKED_OLD))["real_fuel_ids_locked"])

    rows = []
    for r in v5.itertuples():
        src, key = r.entity.split(":", 1)
        if src == "AJFTD":
            grp, fid = gid[key], key
        elif src == "FuelLib":
            n = key.replace("FUELLIB_posf", "")
            grp, fid = (f"posf:{n}" if f"posf:{n}" in posf_groups else f"fuellib:{n}"), key
        else:
            grp, fid = f"av3322:{key}", key
        rows.append(dict(entity=r.entity, source=src, fuel_id=fid, group=grp, freezing_C=float(r.freezing_C),
                         v5_role=r.role, in_old_locked=(src == "AJFTD" and fid in old_locked), comp_source="v5"))
    # new AJFTD fuels: median of uncensored readings
    lg = pd.read_csv(LONG6, low_memory=False, dtype={"fuel_id": str})
    lg = lg[(lg.property == "freezing_point_C") & lg.fuel_id.isin(NEW_IDS)]
    for fid in NEW_IDS:
        d = lg[lg.fuel_id == fid]
        u = d[d.qualifier.isna() & d.value.notna()]
        rows.append(dict(entity=f"AJFTD:{fid}", source="AJFTD", fuel_id=fid, group=gid[fid], freezing_C=float(u.value.median()),
                         v5_role="new", in_old_locked=fid in old_locked, comp_source="v6", n_uncensored=len(u), n_censored=int(d.qualifier.notna().sum())))
    d = pd.DataFrame(rows)
    d["holdout"] = (d.v5_role == "locked") | d.in_old_locked
    d["prio"] = d.source.map(PRIO)
    d = d.sort_values(["group", "prio", "entity"])
    dropped = d[d.duplicated("group", keep="first")]
    d = d[~d.duplicated("group", keep="first")].drop(columns="prio").reset_index(drop=True)
    d["role"] = np.where(d.holdout, "holdout", "dev")
    # a group may not sit in both roles
    assert d.group.is_unique
    d.to_csv(OUT / "dataset.csv", index=False)

    lock = dict(
        rule=__doc__.strip().splitlines()[2:],
        n_total=len(d), n_dev=int((d.role == "dev").sum()), n_holdout=int((d.role == "holdout").sum()),
        holdout_entities=list(d[d.role == "holdout"].entity), holdout_groups=list(d[d.role == "holdout"].group),
        dropped_double_count=list(dropped.entity),
        new_fuels_in_dev=[e for e in d[(d.role == "dev") & (d.v5_role == "new")].entity],
        new_fuels_in_holdout=[e for e in d[(d.role == "holdout") & (d.v5_role == "new")].entity],
        note_fuel_64="freezing = -90 C (single uncensored reading); '<-61' readings are censored and dropped; reported with and without 64.",
        input_sha256=dict(v5_rows=sha(V5ROWS), groups=sha(GROUPS), long_v6=sha(LONG6), old_locked=sha(LOCKED_OLD)),
        dataset_sha256=sha(OUT / "dataset.csv"))
    json.dump(lock, open(OUT / "holdout_lock.json", "w"), indent=1)
    print(d.role.value_counts().to_dict(), "dropped:", list(dropped.entity), "holdout:", lock["holdout_entities"])
    print(d[d.v5_role == "new"][["entity", "group", "freezing_C", "role"]].to_string())


if __name__ == "__main__":
    main()
