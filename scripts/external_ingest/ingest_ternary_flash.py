#!/usr/bin/env python3
"""ingest_ternary_flash.py - Ingest parsed ternary alkane flash points into canonical measurement schema.

Reads:
  - data/interim/external_staging/ternary_flash_liaw_korean.csv
  - data/interim/external_staging/ternary_flash_li2014.csv

Outputs:
  - data/interim/external_staging/canonical_ternary_flash_measurements.csv

NOTE: Staged output ONLY. DOES NOT modify data/unified/ or alter the production pipeline.
"""

from pathlib import Path
import json
import hashlib
import pandas as pd

STAGING_DIR = Path("data/interim/external_staging")

COMPOUNDS_META = {
    "n_nonane": {"inchikey": "BKIMTYGQJCQEES-UHFFFAOYSA-N", "smiles": "CCCCCCCCC", "mw": 128.26},
    "n_decane": {"inchikey": "DIOQZZXIWGTBCZ-UHFFFAOYSA-N", "smiles": "CCCCCCCCCC", "mw": 142.28},
    "n_undecane": {"inchikey": "GHVNFACUYNNKDA-UHFFFAOYSA-N", "smiles": "CCCCCCCCCCC", "mw": 156.31},
    "n_dodecane": {"inchikey": "SNRUBQQJIBEYMU-UHFFFAOYSA-N", "smiles": "CCCCCCCCCCCC", "mw": 170.33},
}


def compute_system_id(inchikeys: list[str]) -> str:
    sorted_keys = sorted(inchikeys)
    return hashlib.md5("".join(sorted_keys).encode("utf-8")).hexdigest()[:12]


def process_df(df: pd.DataFrame, source_label: str) -> list[dict]:
    records = []
    for _, row in df.iterrows():
        # gather components with x > 0
        items = []
        for col_name in ["n_nonane_x", "n_decane_x", "n_undecane_x", "n_dodecane_x"]:
            if col_name in row and row[col_name] > 0.001:
                key = col_name.replace("_x", "")
                meta = COMPOUNDS_META[key]
                items.append({
                    "inchikey": meta["inchikey"],
                    "smiles": meta["smiles"],
                    "mw": meta["mw"],
                    "x": float(row[col_name]),
                })

        if not items:
            continue

        # Sort by inchikey
        items.sort(key=lambda it: it["inchikey"])

        # Normalize mole fractions
        tot_x = sum(it["x"] for it in items)
        for it in items:
            it["x"] = it["x"] / tot_x
            it["mass_raw"] = it["x"] * it["mw"]

        tot_m = sum(it["mass_raw"] for it in items)
        for it in items:
            it["w"] = it["mass_raw"] / tot_m

        inchikeys = [it["inchikey"] for it in items]
        x_mole = [round(it["x"], 6) for it in items]
        w_mass = [round(it["w"], 6) for it in items]
        n_comp = len(items)
        kind = "pure" if n_comp == 1 else ("binary" if n_comp == 2 else "multi")
        sys_id = compute_system_id(inchikeys)

        records.append({
            "source": source_label,
            "property": "flash_point_C",
            "value": float(row["flash_point_C"]),
            "T_K": float("nan"),
            "P_kPa": 101.325,
            "kind": kind,
            "n_comp": n_comp,
            "is_hc": True,
            "system_id": sys_id,
            "comp_inchikeys": json.dumps(inchikeys),
            "x_mole": json.dumps(x_mole),
            "w_mass": json.dumps(w_mass),
            "fuel_id": "",
            "ref": str(row.get("source", source_label)),
            "raw_value": float(row["flash_point_C"]),
            "raw_unit": "C",
        })
    return records


def main():
    records = []
    p_liaw = STAGING_DIR / "ternary_flash_liaw_korean.csv"
    if p_liaw.exists():
        df_liaw = pd.read_csv(p_liaw)
        records.extend(process_df(df_liaw, "liaw_ternary_korean"))

    p_li = STAGING_DIR / "ternary_flash_li2014.csv"
    if p_li.exists():
        df_li = pd.read_csv(p_li)
        records.extend(process_df(df_li, "li2014_ternary"))

    df_out = pd.DataFrame(records).drop_duplicates(subset=["system_id", "comp_inchikeys", "x_mole", "value"])
    out_file = STAGING_DIR / "canonical_ternary_flash_measurements.csv"
    df_out.to_csv(out_file, index=False)
    print(f"Staged {len(df_out)} canonical ternary/binary flash point records to {out_file}")
    print("Breakdown by kind:")
    print(df_out["kind"].value_counts())


if __name__ == "__main__":
    main()
