#!/usr/bin/env python3
"""ingest_kaust_jameel.py - Ingest parsed KAUST octane data into canonical schema format.

Reads:
  - data/interim/external_staging/kaust_jameel_table3_pure.csv
  - data/interim/external_staging/kaust_jameel_table4_blends.csv
  - data/interim/external_staging/kaust_jameel_table5_face.csv
  - data/unified/compounds.csv (for reference / lookup)

Outputs:
  - data/interim/external_staging/kaust_jameel_canonical_measurements.csv
  - data/interim/external_staging/kaust_jameel_new_compounds.csv

NOTE: Staged output ONLY. DOES NOT write to data/unified/ or alter the production pipeline.
"""

from pathlib import Path
import json
import hashlib
import ast
import pandas as pd
from rdkit import Chem

STAGING_DIR = Path("data/interim/external_staging")
COMPOUNDS_REF = Path("data/unified/compounds.csv")

# Standard component name to SMILES map for blends
COMPONENT_SMILES = {
    "n-heptane": "CCCCCCC",
    "i-octane": "CC(C)CC(C)(C)C",
    "toluene": "Cc1ccccc1",
    "1,2,4-trimethylbenzene": "Cc1ccc(C)c(C)c1",
    "cyclopentane": "C1CCCC1",
    "1-hexene": "CCCCC=C",
    "ethanol": "CCO",
}

# Densities at 20 °C (kg/m3) and MW (g/mol) for volume to mass conversion
COMPONENT_DENSITY_MW = {
    "CCCCCCC": {"density": 684.0, "mw": 100.20},
    "CC(C)CC(C)(C)C": {"density": 692.0, "mw": 114.23},
    "Cc1ccccc1": {"density": 867.0, "mw": 92.14},
    "Cc1ccc(C)c(C)c1": {"density": 876.0, "mw": 120.19},
    "C1CCCC1": {"density": 751.0, "mw": 70.13},
    "CCCCC=C": {"density": 673.0, "mw": 84.16},
    "CCO": {"density": 789.0, "mw": 46.07},
}


def get_inchikey_from_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToInchiKey(mol)
    return ""


def compute_system_id(inchikeys: list[str]) -> str:
    sorted_keys = sorted(inchikeys)
    return hashlib.md5("".join(sorted_keys).encode("utf-8")).hexdigest()[:12]


def ingest_pure():
    p = STAGING_DIR / "kaust_jameel_table3_pure.csv"
    if not p.exists():
        return [], []
    df = pd.read_csv(p)
    records = []
    compounds = []

    for _, row in df.iterrows():
        name = row["name"].strip()
        # Clean common names for RDKit parsing
        clean_name = name.lower()
        # Direct SMILES lookup for known standard compounds
        smiles = None
        if clean_name in ["propane"]:
            smiles = "CCC"
        elif clean_name in ["n-butane", "butane"]:
            smiles = "CCCC"
        elif clean_name in ["n-pentane", "pentane"]:
            smiles = "CCCCC"
        elif clean_name in ["n-hexane", "hexane"]:
            smiles = "CCCCCC"
        elif clean_name in ["n-heptane", "heptane"]:
            smiles = "CCCCCCC"
        elif clean_name in ["n-octane", "octane"]:
            smiles = "CCCCCCCC"
        elif clean_name in ["n-nonane", "nonane"]:
            smiles = "CCCCCCCCC"
        elif clean_name in ["n-decane", "decane"]:
            smiles = "CCCCCCCCCC"
        elif clean_name in ["benzene"]:
            smiles = "c1ccccc1"
        elif clean_name in ["toluene"]:
            smiles = "Cc1ccccc1"
        elif clean_name in ["cyclohexane"]:
            smiles = "C1CCCCC1"
        elif clean_name in ["ethanol"]:
            smiles = "CCO"

        # If not in simple list, try to resolve via RDKit if it's already a valid SMILES or IUPAC
        mol = None
        if smiles:
            mol = Chem.MolFromSmiles(smiles)

        if not mol:
            # Generate placeholder or try IUPAC name conversion if possible
            # We record compound name
            continue

        ik = Chem.MolToInchiKey(mol)
        sys_id = compute_system_id([ik])

        # Add RON record
        if pd.notna(row["ron"]):
            records.append({
                "source": "kaust_jameel2018",
                "property": "ron",
                "value": float(row["ron"]),
                "T_K": float("nan"),
                "P_kPa": 101.325,
                "kind": "pure",
                "n_comp": 1,
                "is_hc": not any(atom.GetSymbol() not in ["C", "H"] for atom in mol.GetAtoms()),
                "system_id": sys_id,
                "comp_inchikeys": json.dumps([ik]),
                "x_mole": json.dumps([1.0]),
                "w_mass": json.dumps([1.0]),
                "fuel_id": "",
                "ref": "Jameel2018:Table3",
                "raw_value": float(row["ron"]),
                "raw_unit": "ON",
            })

        # Add MON record
        if pd.notna(row["mon"]):
            records.append({
                "source": "kaust_jameel2018",
                "property": "mon",
                "value": float(row["mon"]),
                "T_K": float("nan"),
                "P_kPa": 101.325,
                "kind": "pure",
                "n_comp": 1,
                "is_hc": not any(atom.GetSymbol() not in ["C", "H"] for atom in mol.GetAtoms()),
                "system_id": sys_id,
                "comp_inchikeys": json.dumps([ik]),
                "x_mole": json.dumps([1.0]),
                "w_mass": json.dumps([1.0]),
                "fuel_id": "",
                "ref": "Jameel2018:Table3",
                "raw_value": float(row["mon"]),
                "raw_unit": "ON",
            })

    return records


def ingest_blends():
    p = STAGING_DIR / "kaust_jameel_table4_blends.csv"
    if not p.exists():
        return []
    df = pd.read_csv(p)
    records = []

    for _, row in df.iterrows():
        comp_dict_raw = row["composition_vol_pct"]
        if isinstance(comp_dict_raw, str):
            try:
                comp_dict = ast.literal_eval(comp_dict_raw)
            except Exception:
                continue
        else:
            comp_dict = comp_dict_raw

        if not comp_dict:
            continue

        # Convert vol fractions to mass and mole fractions
        smiles_vol = []
        for name, vol in comp_dict.items():
            if vol > 0 and name in COMPONENT_SMILES:
                smiles_vol.append((COMPONENT_SMILES[name], vol))

        if not smiles_vol:
            continue

        # Sort by InChIKey
        items = []
        is_hc_all = True
        for smi, vol in smiles_vol:
            ik = get_inchikey_from_smiles(smi)
            mol = Chem.MolFromSmiles(smi)
            if any(atom.GetSymbol() not in ["C", "H"] for atom in mol.GetAtoms()):
                is_hc_all = False
            prop = COMPONENT_DENSITY_MW[smi]
            # mass proportional to vol * density
            mass_raw = vol * prop["density"]
            items.append({
                "inchikey": ik,
                "smiles": smi,
                "vol": vol,
                "mass_raw": mass_raw,
                "mw": prop["mw"],
            })

        items.sort(key=lambda x: x["inchikey"])
        total_mass = sum(x["mass_raw"] for x in items)
        for x in items:
            x["w"] = x["mass_raw"] / total_mass
            x["n_moles"] = x["w"] / x["mw"]
        total_moles = sum(x["n_moles"] for x in items)
        for x in items:
            x["x"] = x["n_moles"] / total_moles

        inchikeys = [x["inchikey"] for x in items]
        w_mass = [round(x["w"], 6) for x in items]
        x_mole = [round(x["x"], 6) for x in items]
        sys_id = compute_system_id(inchikeys)
        n_comp = len(items)
        kind = "binary" if n_comp == 2 else "multi"

        # RON
        if pd.notna(row["ron"]):
            records.append({
                "source": "kaust_jameel2018",
                "property": "ron",
                "value": float(row["ron"]),
                "T_K": float("nan"),
                "P_kPa": 101.325,
                "kind": kind,
                "n_comp": n_comp,
                "is_hc": is_hc_all,
                "system_id": sys_id,
                "comp_inchikeys": json.dumps(inchikeys),
                "x_mole": json.dumps(x_mole),
                "w_mass": json.dumps(w_mass),
                "fuel_id": str(row["blend_name"]),
                "ref": f"Jameel2018:Table4:{row['source_ref']}",
                "raw_value": float(row["ron"]),
                "raw_unit": "ON",
            })

        # MON
        if pd.notna(row["mon"]):
            records.append({
                "source": "kaust_jameel2018",
                "property": "mon",
                "value": float(row["mon"]),
                "T_K": float("nan"),
                "P_kPa": 101.325,
                "kind": kind,
                "n_comp": n_comp,
                "is_hc": is_hc_all,
                "system_id": sys_id,
                "comp_inchikeys": json.dumps(inchikeys),
                "x_mole": json.dumps(x_mole),
                "w_mass": json.dumps(w_mass),
                "fuel_id": str(row["blend_name"]),
                "ref": f"Jameel2018:Table4:{row['source_ref']}",
                "raw_value": float(row["mon"]),
                "raw_unit": "ON",
            })

    return records


def main():
    print("Ingesting KAUST Jameel data into canonical schema...")
    pure_records = ingest_pure()
    blend_records = ingest_blends()
    all_records = pure_records + blend_records

    df_out = pd.DataFrame(all_records)
    out_file = STAGING_DIR / "kaust_jameel_canonical_measurements.csv"
    df_out.to_csv(out_file, index=False)
    print(f"Staged {len(df_out)} canonical octane measurements to {out_file}")
    print("Breakdown by property:")
    print(df_out["property"].value_counts())
    print("Breakdown by kind:")
    print(df_out["kind"].value_counts())


if __name__ == "__main__":
    main()
