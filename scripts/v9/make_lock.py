"""Write results/v9/lock.json once: v9-new SYSTEM ids (no id in v5) with is_v9_locked.  Refuses to overwrite."""
import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.v9.lock import RULE, v9_lock_ids  # noqa: E402
from fuelblend.v7.guard import LOCKS, assert_dev_only  # noqa: E402

OUT = ROOT / "results/v9/lock.json"


def new_system_ids():
    v5 = set(pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=["system_id"]).system_id.astype(str))
    d = pd.read_csv(ROOT / "data/unified_v9/measurements_dedup.csv", low_memory=False)
    return sorted(set(d[d.v9_new].system_id.astype(str)) - v5)


def main():
    if OUT.exists():
        raise SystemExit(f"{OUT} exists; the v9 lock is written once")
    new = new_system_ids()
    ids = sorted(v9_lock_ids(new))
    assert_dev_only(new, locks=tuple(LOCKS))  # new ids must not overlap any earlier lock
    OUT.write_text(json.dumps({"rule": RULE, "ids": ids, "n_new_systems": len(new),
                               "created": str(date.today())}, indent=2))
    print(len(ids), "locked of", len(new))


if __name__ == "__main__":
    main()
