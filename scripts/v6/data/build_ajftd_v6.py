"""Build the v6 AJFTD tables (bulk properties + GCxGC bins) from the 2026-11-20 re-scrape.

The v5 tables in data/real_fuels/v5/ stay frozen (they feed v5 results). This script writes the same tables to
data/real_fuels/v6/ using the unchanged v5 loaders, so the only differences are the fuels added by the re-scrape:

  ajftd_fuels_long.csv, ajftd_fuels_wide.csv, ajftd_build_log.json      (src/fuelblend/data/loaders_ajftd.py)
  ajftd_gcxgc_composition.csv, ajftd_gcxgc_fuels.csv, ajftd_gcxgc_build_log.json
                                                                          (src/fuelblend/data/loaders_ajftd_gcxgc.py)
  ajftd_fuel_groups.csv   one row per fuel_id: posf number parsed from the label and a duplicate group id, so a fuel
                          that appears under two AJFTD ids (e.g. 32 "10264" and 24450 "POSF10264") is never split
                          across CV folds.
  ajftd_v6_vs_v5.json     what changed relative to the frozen v5 tables.

Checks (the script fails loudly if they do not hold):
  * every v5 long row is present unchanged in v6 (the re-scrape only adds fuels);
  * every v5 GCxGC composition row is present unchanged in v6.
AJFTD is login-gated: the raw scrape stays in extras/ (gitignored); only these derived tables are written.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from fuelblend.data import loaders_ajftd as la
from fuelblend.data import loaders_ajftd_gcxgc as lg

ROOT = Path(__file__).resolve().parents[3]
V5 = ROOT / "data" / "real_fuels" / "v5"
V6 = ROOT / "data" / "real_fuels" / "v6"

_POSF = re.compile(r"(?:POSF\s*)?(\d{4,5})")


def posf_number(label: str) -> str | None:
    """'POSF10264' / '10264' / 'CRC003 POSF3962 JETA' -> '10264' / '10264' / '3962'; None if no 4-5 digit id."""
    s = str(label)
    m = re.search(r"POSF\s*(\d{4,5})", s, flags=re.I)
    if m:
        return m.group(1)
    m = re.fullmatch(r"\s*(\d{4,5})\s*", s)
    return m.group(1) if m else None


def fuel_groups(wide: pd.DataFrame, gcxgc_fuels: pd.DataFrame) -> pd.DataFrame:
    ids = pd.concat([wide[["fuel_id", "fuel_label"]], gcxgc_fuels[["fuel_id", "fuel_label"]]]).drop_duplicates("fuel_id")
    ids["posf"] = ids["fuel_label"].map(posf_number)
    ids["dup_group"] = [f"posf:{p}" if isinstance(p, str) else f"id:{int(i)}" for p, i in zip(ids["posf"], ids["fuel_id"])]
    n = ids.groupby("dup_group")["fuel_id"].transform("size")
    ids["n_ids_in_group"] = n
    return ids.sort_values(["dup_group", "fuel_id"]).reset_index(drop=True)


def _rows_preserved(old: pd.DataFrame, new: pd.DataFrame, name: str) -> None:
    cols = list(old.columns)
    sub = new[new["fuel_id"].isin(old["fuel_id"])][cols]
    a = old.sort_values(cols).reset_index(drop=True)
    b = sub.sort_values(cols).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b, check_dtype=False, obj=name)


def main() -> None:
    V6.mkdir(parents=True, exist_ok=True)
    bulk_log = la.build_ajftd_v5(out_dir=V6)

    comp, unparsed = lg.load()
    fuels = lg.fuel_table(comp, wide=V6 / "ajftd_fuels_wide.csv")
    comp.to_csv(V6 / "ajftd_gcxgc_composition.csv", index=False)
    fuels.to_csv(V6 / "ajftd_gcxgc_fuels.csv", index=False)
    ok = fuels[fuels.composition_complete]
    g_log = dict(source=str(lg.RAW.relative_to(ROOT)).replace("\\", "/"), fuels_with_bins=int(len(fuels)),
                 fuels_complete=int(len(ok)), closure_window_pct=list(lg.CLOSURE), bins=int(len(comp)),
                 unparsed_metric_labels=unparsed,
                 complete_with_property={t: int(ok[f"has_{t}"].sum()) for t in lg.TARGETS},
                 note="weight % as printed by AJFTD (AFRL GCxGC); bins are class x carbon number, isomers unresolved")
    (V6 / "ajftd_gcxgc_build_log.json").write_text(json.dumps(g_log, indent=1), encoding="utf-8")

    wide = pd.read_csv(V6 / "ajftd_fuels_wide.csv", low_memory=False)
    groups = fuel_groups(wide, fuels)
    groups.to_csv(V6 / "ajftd_fuel_groups.csv", index=False)

    # --- frozen-v5 consistency checks -------------------------------------------------------------------------
    long5 = pd.read_csv(V5 / "ajftd_fuels_long.csv", low_memory=False)
    long6 = pd.read_csv(V6 / "ajftd_fuels_long.csv", low_memory=False)
    _rows_preserved(long5, long6, "ajftd_fuels_long")
    comp5 = pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")
    comp6 = pd.read_csv(V6 / "ajftd_gcxgc_composition.csv")
    _rows_preserved(comp5, comp6, "ajftd_gcxgc_composition")

    wide5 = pd.read_csv(V5 / "ajftd_fuels_wide.csv", low_memory=False)
    fuels5 = pd.read_csv(V5 / "ajftd_gcxgc_fuels.csv")
    new_gc = sorted(set(fuels.fuel_id) - set(fuels5.fuel_id))
    dup_new = groups[groups.fuel_id.isin(new_gc) & (groups.n_ids_in_group > 1)]
    diff = dict(
        bulk_fuels_v5=int(wide5.fuel_id.nunique()), bulk_fuels_v6=int(wide.fuel_id.nunique()),
        bulk_fuels_added=sorted(int(x) for x in set(wide.fuel_id) - set(wide5.fuel_id)),
        wide_columns_added=[c for c in wide.columns if c not in wide5.columns],
        gcxgc_fuels_v5=int(len(fuels5)), gcxgc_fuels_v6=int(len(fuels)),
        gcxgc_fuels_added=[int(x) for x in new_gc],
        gcxgc_added_with_freezing=[int(x) for x in fuels[fuels.fuel_id.isin(new_gc) & fuels.has_freezing_point_C].fuel_id],
        gcxgc_added_complete=[int(x) for x in fuels[fuels.fuel_id.isin(new_gc) & fuels.composition_complete].fuel_id],
        gcxgc_added_sharing_posf_with_other_ids={
            str(int(r.fuel_id)): sorted(int(x) for x in groups[groups.dup_group == r.dup_group].fuel_id if x != r.fuel_id)
            for r in dup_new.itertuples()},
        checks=["v5 long rows preserved", "v5 gcxgc composition rows preserved"],
        bulk_build_log=bulk_log.get("fuels_with_numeric"),
    )
    (V6 / "ajftd_v6_vs_v5.json").write_text(json.dumps(diff, indent=1), encoding="utf-8")
    print(json.dumps(diff, indent=1))


if __name__ == "__main__":
    main()
