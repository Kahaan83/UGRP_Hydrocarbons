#!/usr/bin/env python3
"""ingest_boehm2022_cold_viscosity.py - Ingest parsed Boehm 2022 cold viscosity/density into canonical schema.

Reads:
  - data/interim/external_staging/boehm2022_raw_extracted.csv

Outputs:
  - data/interim/external_staging/canonical_boehm2022_cold_viscosity.csv
  - data/interim/external_staging/canonical_boehm2022_compounds.csv

NOTE: Staged output ONLY. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import json
import hashlib
import pandas as pd
import numpy as np

STAGING_DIR = Path("data/interim/external_staging")

COMPOUNDS_META = {
    "2,2,4,4,6,8,8-heptamethylnonane": {"smiles": "CC(C)(C)CC(C)(C)CC(C)CC(C)(C)C", "inchikey": "VCLJODPNBNEBKW-UHFFFAOYSA-N", "mw": 226.44, "cls": "alkane"},
    "2,2,4,6,6-pentamethylheptane": {"smiles": "CC(C)CC(C)(C)CC(C)(C)C", "inchikey": "NOFQKTWPZFUCOO-UHFFFAOYSA-N", "mw": 170.34, "cls": "alkane"},
    "2,2,4-trimethylpentane": {"smiles": "CC(C)CC(C)(C)C", "inchikey": "NHTMVDHEPJAVLT-UHFFFAOYSA-N", "mw": 114.23, "cls": "alkane"},
    "2,3-dimethylpentane": {"smiles": "CCC(C)C(C)C", "inchikey": "WGECXQBGLLYSFP-UHFFFAOYSA-N", "mw": 100.20, "cls": "alkane"},
    "2,4-dimethylpentane": {"smiles": "CC(C)CC(C)C", "inchikey": "BZHMBWZPUJHVEE-UHFFFAOYSA-N", "mw": 100.20, "cls": "alkane"},
    "butylcyclohexane": {"smiles": "CCCCC1CCCCC1", "inchikey": "GGBJHURWWWLEQH-UHFFFAOYSA-N", "mw": 140.27, "cls": "cycloalkane"},
    "ethylbenzene": {"smiles": "CCc1ccccc1", "inchikey": "YNQLUTRBYVCPMQ-UHFFFAOYSA-N", "mw": 106.17, "cls": "aromatic"},
    "heptylbenzene": {"smiles": "CCCCCCCc1ccccc1", "inchikey": "LBNXAWYDQUGHGX-UHFFFAOYSA-N", "mw": 176.30, "cls": "aromatic"},
    "hexylbenzene": {"smiles": "c1ccccc1CCCCCC", "inchikey": "LTEQMZWBSYACLV-UHFFFAOYSA-N", "mw": 162.28, "cls": "aromatic"},
    "m-xylene": {"smiles": "Cc1cccc(C)c1", "inchikey": "IVSZLXZYQVIEFR-UHFFFAOYSA-N", "mw": 106.17, "cls": "aromatic"},
    "methylcyclohexane": {"smiles": "CC1CCCCC1", "inchikey": "UAEPNZWRGJTJPN-UHFFFAOYSA-N", "mw": 98.19, "cls": "cycloalkane"},
    "n-hexane": {"smiles": "CCCCCC", "inchikey": "VLKZOEOYAKHREP-UHFFFAOYSA-N", "mw": 86.18, "cls": "alkane"},
    "n-octane": {"smiles": "CCCCCCCC", "inchikey": "TVMXDCGIABBOFY-UHFFFAOYSA-N", "mw": 114.23, "cls": "alkane"},
    "n-pentane": {"smiles": "CCCCC", "inchikey": "OFBQJSOFQDEBGM-UHFFFAOYSA-N", "mw": 72.15, "cls": "alkane"},
    "o-xylene": {"smiles": "Cc1ccccc1C", "inchikey": "CTQNGGLPUBDAKN-UHFFFAOYSA-N", "mw": 106.17, "cls": "aromatic"},
    "trans-decalin": {"smiles": "C1CC[C@@H]2CCCC[C@H]2C1", "inchikey": "NNBZCPXTIHJBJL-MGCOHNPYSA-N", "mw": 138.25, "cls": "cycloalkane"},
    "2-methylpentane": {"smiles": "CCCC(C)C", "inchikey": "AFABGHUZZDYHJO-UHFFFAOYSA-N", "mw": 86.18, "cls": "alkane"},
    "ethylcyclohexane": {"smiles": "CCC1CCCCC1", "inchikey": "IIEWJVIFRVWJOD-UHFFFAOYSA-N", "mw": 112.22, "cls": "cycloalkane"},
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
    seen_compounds = {}

    for _, row in df.iterrows():
        sec = row["section"]
        comp_dict = json.loads(row["composition_json"])
        t_k = float(row["T_K"])
        dens_gcc = float(row["density_g_cc"]) if pd.notna(row["density_g_cc"]) else None
        visc_cst = float(row["viscosity_cSt"]) if pd.notna(row["viscosity_cSt"]) else None

        # Build items
        items = []
        is_real_fuel = False
        fuel_id = None

        for name, vol_pct in comp_dict.items():
            if name in COMPOUNDS_META:
                meta = COMPOUNDS_META[name]
                seen_compounds[meta["inchikey"]] = {
                    "inchikey": meta["inchikey"],
                    "smiles": meta["smiles"],
                    "mw": meta["mw"],
                    "cls": meta["cls"],
                    "name": name,
                    "is_hc": True
                }
                items.append({
                    "name": name,
                    "inchikey": meta["inchikey"],
                    "mw": meta["mw"],
                    "v": vol_pct,
                })
            else:
                # Real fuel / POSF lump
                is_real_fuel = True
                fuel_id = f"POSF:{name}"
                lump_key = f"LUMP:{name}"
                items.append({
                    "name": name,
                    "inchikey": lump_key,
                    "mw": 150.0,
                    "v": vol_pct,
                })

        if not items:
            continue

        items.sort(key=lambda x: x["inchikey"])
        inchikeys = [it["inchikey"] for it in items]
        n_comp = len(items)

        # Volume fractions to mass and mole fractions
        tot_v = sum(it["v"] for it in items)
        v_fracs = [it["v"] / tot_v for it in items]
        w_fracs = v_fracs  # Approx volume ~ mass for similar hydrocarbons
        moles_unnorm = [w / it["mw"] for w, it in zip(w_fracs, items)]
        tot_moles = sum(moles_unnorm)
        x_fracs = [m / tot_moles for m in moles_unnorm]

        if is_real_fuel:
            kind = "real_fuel"
            system_id = f"fuel:{fuel_id}"
        elif n_comp == 1:
            kind = "pure"
            system_id = compute_system_id(inchikeys)
        elif n_comp == 2:
            kind = "binary"
            system_id = compute_system_id(inchikeys)
        else:
            kind = "multi"
            system_id = compute_system_id(inchikeys)

        base_rec = {
            "source": "boehm2022",
            "T_K": t_k,
            "P_kPa": 101.325,
            "kind": kind,
            "n_comp": n_comp,
            "is_hc": True,
            "system_id": system_id,
            "comp_inchikeys": json.dumps(inchikeys),
            "x_mole": json.dumps([round(f, 6) for f in x_fracs]),
            "w_mass": json.dumps([round(f, 6) for f in w_fracs]),
            "fuel_id": fuel_id,
            "ref": "Boehm2022:Table:ThisWork",
        }

        # 1. Kinematic Viscosity
        if visc_cst is not None and not np.isnan(visc_cst):
            r_visc = base_rec.copy()
            r_visc.update({
                "property": "viscosity_kin_cSt",
                "value": round(visc_cst, 6),
                "raw_value": visc_cst,
                "raw_unit": "cSt",
            })
            records.append(r_visc)

        # 2. Density (convert g/cc -> kg/m3)
        if dens_gcc is not None and not np.isnan(dens_gcc):
            r_dens = base_rec.copy()
            r_dens.update({
                "property": "density_kg_m3",
                "value": round(dens_gcc * 1000.0, 3),
                "raw_value": dens_gcc,
                "raw_unit": "g/cc",
            })
            records.append(r_dens)

    canon_df = pd.DataFrame(records)
    out_csv = STAGING_DIR / "canonical_boehm2022_cold_viscosity.csv"
    canon_df.to_csv(out_csv, index=False)
    print(f"Ingested {len(canon_df)} canonical rows to {out_csv}")
    print(canon_df.groupby(["property", "kind"])["value"].count().unstack(fill_value=0))

    comp_df = pd.DataFrame(list(seen_compounds.values()))
    comp_csv = STAGING_DIR / "canonical_boehm2022_compounds.csv"
    comp_df.to_csv(comp_csv, index=False)
    print(f"Saved {len(comp_df)} compounds to {comp_csv}")

if __name__ == "__main__":
    main()
