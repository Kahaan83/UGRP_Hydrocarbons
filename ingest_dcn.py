"""Ingest DCN data from Leenhouts benchmark (T15).
Pure: 229 -> 80 HC
Mix: 484 -> 346 HC (213 binary, 74 ternary, 59 quaternary)
Outputs data/processed/dcn_hc.csv
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from common import mol_from, is_hydrocarbon, mol_weight, mole_to_mass, mass_to_mole

DIR = Path(__file__).resolve().parent
IN_DIR = DIR / "data/interim/leenhouts_benchmark/datasets/datasets/fuel_ignition_numbers/processed_data"
OUT_FILE = DIR / "data/processed/dcn_hc.csv"


def build_dcn():
    pure_df = pd.read_csv(IN_DIR / "dcn_pure_exp.csv")
    mix_df = pd.read_csv(IN_DIR / "dcn_mix_exp.csv")

    rows = []

    # 1. Pure compounds
    for _, r in pure_df.iterrows():
        inchi = str(r["component_inchi_0"]).strip()
        m = mol_from(inchi=inchi)
        if not is_hydrocarbon(m):
            continue
        smiles = Chem.MolToSmiles(m)
        inchikey = Chem.MolToInchiKey(m)
        rows.append({
            "kind": "pure",
            "n_comp": 1,
            "smiles": json.dumps([smiles]),
            "inchi": json.dumps([inchi]),
            "inchikey": json.dumps([inchikey]),
            "x_mole": json.dumps([1.0]),
            "w_mass": json.dumps([1.0]),
            "dcn": float(r["value"]),
            "unit": "dimensionless",
            "tier": "A",
            "source": "leenhouts_benchmark:dcn_pure_exp",
        })

    # 2. Mixtures
    for _, r in mix_df.iterrows():
        inchis, fracs = [], []
        for j in range(4):
            col_i, col_f = f"component_inchi_{j}", f"fraction_{j}"
            if col_i in r and pd.notna(r[col_i]) and str(r[col_i]).strip():
                inchis.append(str(r[col_i]).strip())
                fracs.append(float(r[col_f]))

        mols = [mol_from(inchi=inc) for inc in inchis]
        if any(m is None or not is_hydrocarbon(m) for m in mols):
            continue

        smiles = [Chem.MolToSmiles(m) for m in mols]
        inchikeys = [Chem.MolToInchiKey(m) for m in mols]
        mws = [mol_weight(m) for m in mols]

        ftype = str(r["fraction_type"]).lower()
        if "molar" in ftype or "mole" in ftype:
            tot_x = sum(fracs)
            x = [v / tot_x for v in fracs]
            w = mole_to_mass(x, mws).tolist()
        elif "mass" in ftype or "weight" in ftype:
            tot_w = sum(fracs)
            w = [v / tot_w for v in fracs]
            x = mass_to_mole(w, mws).tolist()
        else:
            raise ValueError(f"Unknown fraction_type: {ftype}")

        rows.append({
            "kind": "mixture",
            "n_comp": len(inchis),
            "smiles": json.dumps(smiles),
            "inchi": json.dumps(inchis),
            "inchikey": json.dumps(inchikeys),
            "x_mole": json.dumps([round(float(v), 6) for v in x]),
            "w_mass": json.dumps([round(float(v), 6) for v in w]),
            "dcn": float(r["value"]),
            "unit": "dimensionless",
            "tier": "A",
            "source": "leenhouts_benchmark:dcn_mix_exp",
        })

    out_df = pd.DataFrame(rows)
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_FILE, index=False)
    print(f"Wrote {len(out_df)} records to {OUT_FILE}")
    print(f"  Pure HC: {(out_df['kind'] == 'pure').sum()}")
    print(f"  Mix HC:  {(out_df['kind'] == 'mixture').sum()}")
    print("  Mix by n_comp:", out_df[out_df["kind"] == "mixture"]["n_comp"].value_counts().to_dict())
    return out_df


if __name__ == "__main__":
    build_dcn()
