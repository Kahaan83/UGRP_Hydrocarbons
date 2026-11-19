#!/usr/bin/env python3
"""ingest_leenhouts_mon_multicomponent.py - Ingest 292 multicomponent MON blend rows from Leenhouts dataset.

Source: Leenhouts et al. / Kuzhagaliyeva et al. (Zenodo 19914760).
File: data/interim/leenhouts_dcn/datasets/fuel_ignition_numbers/processed_data/mon_mix_exp.csv

Outputs:
  - data/interim/external_staging/canonical_leenhouts_mon_mix_measurements.csv
  - data/interim/external_staging/canonical_leenhouts_mon_compounds.csv

NOTE: Staged output ONLY. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import json
import hashlib
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

RAW_CSV = Path("data/interim/leenhouts_dcn/datasets/fuel_ignition_numbers/processed_data/mon_mix_exp.csv")
STAGING_DIR = Path("data/interim/external_staging")


def compute_system_id(inchikeys: list[str]) -> str:
    sorted_keys = sorted(inchikeys)
    return hashlib.md5("".join(sorted_keys).encode("utf-8")).hexdigest()[:12]


def main():
    if not RAW_CSV.exists():
        print(f"ERROR: {RAW_CSV} not found.")
        return

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(RAW_CSV)
    print(f"Reading {len(df)} MON mixture rows from {RAW_CSV}...")

    # Cache for InChI -> (InChIKey, MW, is_hc)
    cache = {}

    records = []
    unique_compounds = {}

    for idx, row in df.iterrows():
        val = float(row["value"])
        if pd.isna(val):
            continue

        items = []
        is_hc_blend = True

        for i in range(85):
            c_col = f"component_inchi_{i}"
            f_col = f"fraction_{i}"
            if c_col in row and pd.notna(row[c_col]) and f_col in row and pd.notna(row[f_col]):
                inchi = str(row[c_col]).strip()
                x_frac = float(row[f_col])
                if x_frac <= 0.0:
                    continue

                if inchi not in cache:
                    mol = Chem.MolFromInchi(inchi)
                    if mol:
                        ik = Chem.MolToInchiKey(mol)
                        mw = float(ExactMolWt(mol))
                        smi = Chem.MolToSmiles(mol)
                        is_hc = not any(atom.GetSymbol() not in ["C", "H"] for atom in mol.GetAtoms())
                        cache[inchi] = (ik, smi, mw, is_hc)
                    else:
                        cache[inchi] = (None, None, None, None)

                ik, smi, mw, is_hc = cache[inchi]
                if ik:
                    if not is_hc:
                        is_hc_blend = False
                    items.append({
                        "inchikey": ik,
                        "smiles": smi,
                        "mw": mw,
                        "x": x_frac,
                    })
                    if ik not in unique_compounds:
                        unique_compounds[ik] = {
                            "inchikey": ik,
                            "smiles": smi,
                            "inchi": inchi,
                            "mw": mw,
                            "is_hc": is_hc,
                        }

        if not items:
            continue

        # Sort items by InChIKey
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
        kind = "binary" if n_comp == 2 else "multi"
        sys_id = compute_system_id(inchikeys)

        records.append({
            "source": "leenhouts_dcn",
            "property": "mon",
            "value": val,
            "T_K": float("nan"),
            "P_kPa": 101.325,
            "kind": kind,
            "n_comp": n_comp,
            "is_hc": is_hc_blend,
            "system_id": sys_id,
            "comp_inchikeys": json.dumps(inchikeys),
            "x_mole": json.dumps(x_mole),
            "w_mass": json.dumps(w_mass),
            "fuel_id": f"leenhouts_mon_mix_{idx}",
            "ref": "Kuzhagaliyeva2022:mon_mix_exp",
            "raw_value": val,
            "raw_unit": "ON",
        })

    df_out = pd.DataFrame(records)
    out_meas = STAGING_DIR / "canonical_leenhouts_mon_mix_measurements.csv"
    df_out.to_csv(out_meas, index=False)
    print(f"Staged {len(df_out)} canonical MON mixture rows -> {out_meas}")
    print(f"Number of unique multi-component systems: {df_out['system_id'].nunique()}")

    df_cmps = pd.DataFrame(list(unique_compounds.values()))
    out_cmps = STAGING_DIR / "canonical_leenhouts_mon_compounds.csv"
    df_cmps.to_csv(out_cmps, index=False)
    print(f"Staged {len(df_cmps)} compounds -> {out_cmps}")


if __name__ == "__main__":
    main()
