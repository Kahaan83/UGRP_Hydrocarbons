#!/usr/bin/env python3
"""ingest_sjtu_ai4fuel.py - Ingest parsed SJTU AI4fuel data into canonical schema format.

Reads:
  - data/interim/external_staging/ai4fuel/*.csv

Outputs:
  - data/interim/external_staging/canonical_ai4fuel_measurements.csv
  - data/interim/external_staging/canonical_ai4fuel_compounds.csv

NOTE: Staged output ONLY. DOES NOT modify data/unified/ or alter the production pipeline.
"""

from pathlib import Path
import json
import hashlib
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

IN_DIR = Path("data/interim/external_staging/ai4fuel")
STAGING_DIR = Path("data/interim/external_staging")

# Properties to include in canonical measurements
CANONICAL_PROPS = {
    "ron": ("ron", float("nan"), 101.325, "ON", 1.0),
    "mon": ("mon", float("nan"), 101.325, "ON", 1.0),
    "cetane_number_D613": ("cetane_number_D613", float("nan"), 101.325, "-", 1.0),
    "dcn": ("dcn", float("nan"), 101.325, "-", 1.0),
    "flash_point_C": ("flash_point_C", float("nan"), 101.325, "C", 1.0),
    "melting_point_C": ("freezing_point_C", float("nan"), 101.325, "C", 1.0),
    "viscosity_dyn_mPas": ("viscosity_dyn_mPas", 293.15, 101.325, "mPa*s", 1.0),
    "density_g_cm3": ("density_kg_m3", 293.15, 101.325, "kg/m3", 1000.0),
    "vapor_pressure_kPa": ("vapor_pressure_kPa", 298.15, 101.325, "kPa", 1.0),
    "lfl_vol_pct": ("lfl_vol_pct", 298.15, 101.325, "vol%", 1.0),
    "lhv_MJ_kg": ("heat_of_combustion_net_MJ_kg", 298.15, 101.325, "MJ/kg", 1.0),
}


def compute_system_id(inchikey: str) -> str:
    return hashlib.md5(inchikey.encode("utf-8")).hexdigest()[:12]


def main():
    if not IN_DIR.exists():
        print(f"ERROR: {IN_DIR} not found. Run parse_sjtu_ai4fuel.py first.")
        return

    mol_cache = {}  # smiles -> (inchikey, inchi, formula, mw, is_hc)
    measurements = []
    compounds = {}

    csv_files = list(IN_DIR.glob("ai4fuel_*.csv"))
    print(f"Ingesting {len(csv_files)} AI4fuel property files...")

    for csv_path in csv_files:
        prop_file_key = csv_path.stem.replace("ai4fuel_", "")
        if prop_file_key not in CANONICAL_PROPS:
            continue

        canon_name, def_t_k, def_p_kpa, unit_str, scale_factor = CANONICAL_PROPS[prop_file_key]
        df = pd.read_csv(csv_path)

        for _, row in df.iterrows():
            smi = str(row["smiles"]).strip()
            if not smi or smi == "nan":
                continue

            if smi not in mol_cache:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    mol_cache[smi] = (None, None, None, None, None)
                else:
                    ik = Chem.MolToInchiKey(mol)
                    inchi = Chem.MolToInchi(mol)
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    mw = float(ExactMolWt(mol))
                    is_hc = not any(atom.GetSymbol() not in ["C", "H"] for atom in mol.GetAtoms())
                    mol_cache[smi] = (ik, inchi, formula, mw, is_hc)

            ik, inchi, formula, mw, is_hc = mol_cache[smi]
            if not ik:
                continue

            # Record compound
            if ik not in compounds:
                compounds[ik] = {
                    "inchikey": ik,
                    "smiles": smi,
                    "inchi": inchi,
                    "formula": formula,
                    "mw": mw,
                    "is_hc": is_hc,
                    "cls": str(row.get("class", "other")),
                    "name": str(row.get("name", "")),
                    "cas": str(row.get("cas", "")),
                }

            val = float(row["value"]) * scale_factor
            sys_id = compute_system_id(ik)

            measurements.append({
                "source": "sjtu_ai4fuel",
                "property": canon_name,
                "value": round(val, 4),
                "T_K": def_t_k,
                "P_kPa": def_p_kpa,
                "kind": "pure",
                "n_comp": 1,
                "is_hc": is_hc,
                "system_id": sys_id,
                "comp_inchikeys": json.dumps([ik]),
                "x_mole": json.dumps([1.0]),
                "w_mass": json.dumps([1.0]),
                "fuel_id": "",
                "ref": "AI4fuel:NatureCommun2026",
                "raw_value": float(row["value"]),
                "raw_unit": unit_str,
            })

    df_meas = pd.DataFrame(measurements)
    out_meas = STAGING_DIR / "canonical_ai4fuel_measurements.csv"
    df_meas.to_csv(out_meas, index=False)
    print(f"Staged {len(df_meas)} canonical measurements -> {out_meas}")
    print("Breakdown by property:")
    print(df_meas["property"].value_counts())
    print("Hydrocarbons vs non-hydrocarbons:")
    print(df_meas["is_hc"].value_counts())

    df_comp = pd.DataFrame(list(compounds.values()))
    out_comp = STAGING_DIR / "canonical_ai4fuel_compounds.csv"
    df_comp.to_csv(out_comp, index=False)
    print(f"Staged {len(df_comp)} unique compounds -> {out_comp}")


if __name__ == "__main__":
    main()
