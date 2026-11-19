#!/usr/bin/env python3
"""ingest_crc_face1_diesel.py - Ingest parsed CRC FACE-1 diesel properties into canonical real fuel schema.

Reads:
  - data/interim/external_staging/crc_face1_diesel_properties.csv
  - data/interim/external_staging/crc_face1_diesel_composition_table52.csv

Outputs:
  - data/interim/external_staging/canonical_crc_face1_diesel_measurements.csv

NOTE: Staged output ONLY. DOES NOT modify data/unified/ or alter the production pipeline.
"""

from pathlib import Path
import json
import pandas as pd

STAGING_DIR = Path("data/interim/external_staging")


def main():
    p_props = STAGING_DIR / "crc_face1_diesel_properties.csv"
    if not p_props.exists():
        print(f"ERROR: {p_props} not found. Run parse_crc_face1_diesel.py first.")
        return

    df_props = pd.read_csv(p_props)
    records = []

    for _, row in df_props.iterrows():
        fuel_id = str(row["fuel_id"])
        prop = str(row["property"])
        val = float(row["value"])
        unit = str(row["unit"])
        t_k = float(row["T_K"]) if pd.notna(row["T_K"]) else float("nan")

        records.append({
            "source": "crc_face1",
            "property": prop,
            "value": val,
            "T_K": t_k,
            "P_kPa": 101.325,
            "kind": "real_fuel",
            "n_comp": float("nan"),
            "is_hc": True,
            "system_id": f"fuel:{fuel_id}",
            "comp_inchikeys": json.dumps([]),
            "x_mole": json.dumps([]),
            "w_mass": json.dumps([]),
            "fuel_id": fuel_id,
            "ref": f"CRC_FACE-1:{row['source_table']}:{row['lab']}",
            "raw_value": val,
            "raw_unit": unit,
        })

    df_out = pd.DataFrame(records)
    out_file = STAGING_DIR / "canonical_crc_face1_diesel_measurements.csv"
    df_out.to_csv(out_file, index=False)
    print(f"Staged {len(df_out)} canonical diesel measurements to {out_file}")
    print("Breakdown by property:")
    print(df_out["property"].value_counts())


if __name__ == "__main__":
    main()
