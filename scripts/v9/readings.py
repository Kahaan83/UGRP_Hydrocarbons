"""v9 readings R9-1..R9-4 under the pass rules of docs/v9_predeclaration.md.

--set dev    : data new in v9 that is NOT in the v9 lock (and no id in any earlier lock).
--set v9_lock: the v9 lock, scored once.
v9 added too little data for the model-level bootstraps to be defined (floors of the predeclaration unmet);
each reading therefore reports the counts that decide the rule and the verdict that follows from it.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.v7.guard import LOCKS, assert_dev_only, canon, locked_ids  # noqa: E402
from fuelblend.v9.lock import is_v9_locked  # noqa: E402

FLOOR = 8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["dev", "v9_lock"], required=True)
    a = ap.parse_args()
    lock = set(json.loads((ROOT / "results/v9/lock.json").read_text())["ids"])
    v5 = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=["system_id"]).system_id.astype(str)
    d = pd.read_csv(ROOT / "data/unified_v9/measurements_dedup.csv", low_memory=False)
    new = d[d.v9_new].copy()
    new["system_id"] = new.system_id.astype(str)
    new_sys = set(new.system_id) - set(v5)
    sel = new[new.system_id.isin(lock)] if a.set == "v9_lock" else new[~new.system_id.isin(lock)]
    if a.set == "dev":
        earlier = set().union(*(locked_ids(n) for n in LOCKS))
        sel = sel[~sel.system_id.map(canon).isin(earlier)]  # fuels already scored in v5 to v7 locks stay out of dev
        assert_dev_only(sorted(set(sel.system_id)), locks=tuple(LOCKS))
        assert not (set(sel.system_id) & lock)
    comp = set(pd.read_csv(ROOT / "data/unified_v9/real_fuel_composition_v9.csv", usecols=["fuel_id"]).fuel_id.astype(str))
    rf = sel[(sel.kind == "real_fuel")]
    has_comp = rf[rf.system_id.str.extract(r"^fuel:(\d+)$")[0].isin(comp)]
    frz = sel[sel.property == "freezing_point_C"]
    visc = sel[sel.property.str.startswith("viscosity")]
    cold = has_comp[has_comp.property == "viscosity_kin_cSt"]
    rows = [
        dict(reading="R9-1", n=frz.system_id.nunique(), floor=FLOOR,
             verdict="status: no skill shown" + (" (lock floor of 8 fuels unmet, descriptive only)" if a.set == "v9_lock" else ""),
             note="no verified freezing fuel was added (Edwards 2017 unverifiable)"),
        dict(reading="R9-2", n=visc[visc.system_id.isin(new_sys) & ~visc.kind.eq("real_fuel")].system_id.nunique(),
             floor="CI gain > 0", verdict="more data helps: not shown (retrain not run)",
             note=f"{len(visc)} new viscosity groups; new viscosity systems are real fuels without composition, "
                  "so no new composition-bearing training system exists"),
        dict(reading="R9-3", n=has_comp.system_id.nunique(), floor=FLOOR,
             verdict="status: no skill shown (floor unmet, descriptive only)",
             note="neat fuels with composition: " + ",".join(sorted(has_comp.system_id.unique()))),
        dict(reading="R9-4", n=sel[sel.kind.isin(["binary", "multi", "pure"])].system_id.nunique(), floor="descriptive",
             verdict="descriptive, no verdict change",
             note=f"{len(sel[sel.kind.isin(['binary', 'multi', 'pure'])])} groups"),
    ]
    out = ROOT / "results/v9/readings"
    out.mkdir(parents=True, exist_ok=True)
    for r in rows:
        pd.DataFrame([r]).to_csv(out / f"{r['reading']}_{a.set}.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
