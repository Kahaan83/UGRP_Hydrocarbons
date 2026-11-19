#!/usr/bin/env python3
"""ingest_boehm2022_cold_viscosity.py - Ingest parsed Boehm 2022 cold viscosity/density into canonical schema.

Reads:
  - data/interim/external_staging/boehm2022_raw_extracted.csv

Outputs:
  - data/interim/external_staging/canonical_boehm2022_cold_viscosity.csv

NOTE: Staged output ONLY. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import json
import hashlib
import pandas as pd

STAGING_DIR = Path("data/interim/external_staging")

COMPOUNDS_META = {
    "c1": {"name": "2,2,4,4,6,8,8-heptamethylnonane", "inchikey": "HNJBEVLQSNBDRC-UHFFFAOYSA-N", "smiles": "CC(C)CC(C)CC(C)(C)CC(C)(C)C", "mw": 226.44},
    "c2": {"name": "methylcyclohexane", "inchikey": "UAEPNZWRGJTJPN-UHFFFAOYSA-N", "smiles": "CC1CCCCC1", "mw": 98.19},
    "c3": {"name": "hexylbenzene", "inchikey": "QWTDNUCVQCZILF-UHFFFAOYSA-N", "smiles": "CCCCCCc1ccccc1", "mw": 162.27},
}


def compute_system_id(inchikeys: list[str]) -> str:
    sorted_keys = sorted(inchikeys)
    return hashlib.md5("".join(sorted_keys).encode("utf-8")).hexdigest()[:12]


def main():
    p = STAGING_DIR / "boehm2022_raw_extracted.csv"
    if not p.exists():
        print(f"ERROR: {p} not found. Run parse_boehm2022_cold_viscosity.py first.")
        return

    df = pd.read_csv(p)
    records = []

    for _, row in df.iterrows():
        v1 = float(row["c1_heptamethylnonane_pct"])
        v2 = float(row["c2_methylcyclohexane_pct"])
        v3 = float(row["c3_hexylbenzene_pct"])
        t_k = float(row["T_K"])

        items = []
        for key, val in [("c1", v1), ("c2", v2), ("c3", v3)]:
            if val > 0.001:
                meta = COMPOUNDS_META[key]
                items.append({
                    "inchikey": meta["inchikey"],
                    "smiles": meta["smiles"],
                    "mw": meta["mw"],
                    "w": val,
                })

        if not items:
            continue

        items.sort(key=lambda x: x["inchikey"])
        tot_w = sum(it["w"] for it in items)
        for it in items:
            it["w"] = it["w"] / tot_w
            it["n_moles"] = it["w"] / it["mw"]

        tot_moles = sum(it["n_moles"] for it in items)
        for it in items:
            it["x"] = it["n_moles"] / tot_moles

        inchikeys = [it["inchikey"] for it in items]
        x_mole = [round(it["x"], 6) for it in items]
        w_mass = [round(it["w"], 6) for it in items]
        n_comp = len(items)
        kind = "pure" if n_comp == 1 else ("binary" if n_comp == 2 else "multi")
        sys_id = compute_system_id(inchikeys)

        # Kinematic viscosity record
        visc = float(row["viscosity_cSt"])
        if pd.notna(visc) and visc > 0:
            records.append({
                "source": "boehm2022_heat_lab",
                "property": "viscosity_kin_cSt",
                "value": visc,
                "T_K": t_k,
                "P_kPa": 101.325,
                "kind": kind,
                "n_comp": n_comp,
                "is_hc": True,
                "system_id": sys_id,
                "comp_inchikeys": json.dumps(inchikeys),
                "x_mole": json.dumps(x_mole),
                "w_mass": json.dumps(w_mass),
                "fuel_id": "",
                "ref": "Boehm2022:Table:ThisWork",
                "raw_value": visc,
                "raw_unit": "cSt",
            })

        # Density record
        dens = float(row["density_g_cc"])
        if pd.notna(dens) and dens > 0:
            records.append({
                "source": "boehm2022_heat_lab",
                "property": "density_kg_m3",
                "value": round(dens * 1000.0, 2),
                "T_K": t_k,
                "P_kPa": 101.325,
                "kind": kind,
                "n_comp": n_comp,
                "is_hc": True,
                "system_id": sys_id,
                "comp_inchikeys": json.dumps(inchikeys),
                "x_mole": json.dumps(x_mole),
                "w_mass": json.dumps(w_mass),
                "fuel_id": "",
                "ref": "Boehm2022:Table:ThisWork",
                "raw_value": dens,
                "raw_unit": "g/cc",
            })

    df_out = pd.DataFrame(records)
    out_file = STAGING_DIR / "canonical_boehm2022_cold_viscosity.csv"
    df_out.to_csv(out_file, index=False)
    print(f"Staged {len(df_out)} canonical cold viscosity & density measurements -> {out_file}")
    print("Breakdown by property:")
    print(df_out["property"].value_counts())
    print("Breakdown by T_K:")
    print(df_out["T_K"].value_counts().sort_index())


if __name__ == "__main__":
    main()
