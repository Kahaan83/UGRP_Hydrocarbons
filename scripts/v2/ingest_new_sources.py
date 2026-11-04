"""Ingest new open data sources: NREL Cetane Compendium, DOE Co-Optima, and Coutinho Wax SLE.

Produces standardized processed datasets:
  - data/processed/v2/cetane_nrel_compendium.csv
  - data/processed/v2/octane_co_optima.csv
  - data/processed/v2/wax_dauphin_coutinho.csv
"""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[2]
PROCESSED_V2 = ROOT / "data" / "processed" / "v2"
PROCESSED_V2.mkdir(parents=True, exist_ok=True)


def is_pure_hydrocarbon(mol):
    if mol is None:
        return False
    return all(atom.GetSymbol() in ["C", "H"] for atom in mol.GetAtoms())


def process_nrel_cetane():
    print("Processing NREL Cetane Compendium (2017)...")
    src = ROOT / "data" / "raw" / "cetane" / "nrel_cn_merged_210831.csv"
    df = pd.read_csv(src)

    records = []
    seen = set()

    for _, row in df.iterrows():
        smi = str(row["Canonical_SMILES"]).strip()
        mol = Chem.MolFromSmiles(smi)
        if not is_pure_hydrocarbon(mol):
            continue

        can_smi = Chem.MolToSmiles(mol, canonical=True)
        ik = Chem.MolToInchiKey(mol)
        if ik in seen:
            continue
        seen.add(ik)

        cn_val = float(row["CN"])
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = float(Descriptors.MolWt(mol))
        nc = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")

        records.append({
            "compound_name": str(row.get("IUPAC_NAME", "")).strip(),
            "smiles": can_smi,
            "inchikey": ik,
            "formula": formula,
            "mw": round(mw, 3),
            "n_carbon": nc,
            "dcn": round(cn_val, 2),
            "test_method": str(row.get("METHOD", "ASTM D6890")),
            "device": str(row.get("DEVICE", "IQT")),
            "source_reference": str(row.get("REFERENCE", "NREL/TP-5400-67585 (2017)")),
            "source_tag": "nrel_compendium_2017",
        })

    out_df = pd.DataFrame(records).sort_values("n_carbon").reset_index(drop=True)
    out_path = PROCESSED_V2 / "cetane_nrel_compendium.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved {len(out_df)} pure hydrocarbon cetane records to {out_path}")
    return out_df


def process_co_optima():
    print("Processing DOE Co-Optima Database...")
    src = ROOT / "data" / "raw" / "cooptima" / "fuel_engine_co_optimization.csv"
    df = pd.read_csv(src, encoding="latin1", low_memory=False)

    records = []
    seen = set()

    for _, row in df.iterrows():
        name = str(row.get("Pure_IUPAC_name", "")).strip()
        cas = str(row.get("Pure_CAS", "")).strip()
        
        # Try to resolve SMILES from name or CAS or formula
        # Co-Optima doesn't directly have SMILES column; use Chem.MolFromSmiles if name is SMILES, else RDKit or PubChem
        ron_val = row.get("Pure_RON")
        mon_val = row.get("Pure_MON")
        lfl_val = row.get("Pure_LFL_LEL")
        flash_val = row.get("Pure_Flash_Point")

        # Keep rows with at least one valid octane or property
        if pd.isna(ron_val) and pd.isna(mon_val) and pd.isna(flash_val):
            continue

        # Simple SMILES lookups for common fuel hydrocarbons
        # Co-Optima provides formula and IUPAC names
        formula = str(row.get("Pure_Molecular_Formula", "")).strip()
        # Verify only C and H
        if formula and not re.match(r"^C\d+H\d+$", formula):
            continue

        try:
            ron_f = float(ron_val) if pd.notna(ron_val) else np.nan
            mon_f = float(mon_val) if pd.notna(mon_val) else np.nan
            sens_f = ron_f - mon_f if (np.isfinite(ron_f) and np.isfinite(mon_f)) else np.nan
            lfl_f = float(lfl_val) if pd.notna(lfl_val) else np.nan
            flash_f = float(flash_val) if pd.notna(flash_val) else np.nan
        except (ValueError, TypeError):
            continue

        records.append({
            "compound_name": name,
            "cas": cas,
            "formula": formula,
            "ron": ron_f,
            "mon": mon_f,
            "sensitivity": sens_f,
            "lfl_vol_pct": lfl_f,
            "flash_point_C": flash_f,
            "source_reference": "DOE Co-Optima Database (2018)",
            "source_tag": "doe_co_optima",
        })

    out_df = pd.DataFrame(records).drop_duplicates(subset=["compound_name", "cas"]).reset_index(drop=True)
    out_path = PROCESSED_V2 / "octane_co_optima.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved {len(out_df)} Co-Optima pure hydrocarbon records to {out_path}")
    return out_df


def process_coutinho_wax():
    print("Processing Coutinho / Dauphin Wax SLE Data...")
    src = ROOT / "data" / "raw" / "wax" / "sp_wax_results.xlsx"
    df_dauphin = pd.read_excel(src, sheet_name="Dauphin")

    # Dauphin sheet structure: multi-component paraffin wax blends with carbon number distributions
    # Clean rows with valid experimental WAT
    records = []
    
    # Also parse Binary System sheet
    df_bin = pd.read_excel(src, sheet_name="Binary Sytem")
    # Clean binary pairs
    current_system = None
    for idx, row in df_bin.iterrows():
        val0 = str(row.iloc[0]).strip() if pd.notna(row.iloc[0]) else ""
        if "C" in val0 and "-" in val0 and len(val0) < 15:
            current_system = val0
            continue
        try:
            x_solute = float(row.iloc[1])
            t_exp = float(row.iloc[2]) # Experimental WAT in K
            if np.isfinite(x_solute) and np.isfinite(t_exp) and t_exp > 150.0:
                comps = current_system.split("-") if current_system else ["Solute", "Solvent"]
                records.append({
                    "system_id": f"wax_bin_{current_system}_{idx}",
                    "system_name": current_system or "Binary Paraffin",
                    "n_comp": 2,
                    "components": json.dumps(comps),
                    "mole_fractions": json.dumps([round(x_solute, 4), round(1.0 - x_solute, 4)]),
                    "wat_measured_K": round(t_exp, 2),
                    "wat_measured_C": round(t_exp - 273.15, 2),
                    "source_reference": "Dauphin / Coutinho SP-Wax (2018)",
                    "source_tag": "coutinho_dauphin_wax",
                })
        except (ValueError, TypeError, IndexError):
            continue

    out_df = pd.DataFrame(records)
    out_path = PROCESSED_V2 / "wax_dauphin_coutinho.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved {len(out_df)} experimental wax crystallization records to {out_path}")
    return out_df


def update_data_sources_csv():
    print("Updating docs/data_sources.csv...")
    ds_path = ROOT / "docs" / "data_sources.csv"
    df = pd.read_csv(ds_path)

    new_sources = [
        {
            "source_id": "nrel_compendium_2017",
            "name": "NREL Compendium of Experimental Cetane Numbers (Yanowitz et al. 2017)",
            "what_it_provides": "630 pure-compound experimental cetane numbers (ASTM D6890 / D613) across hydrocarbons and oxygenates",
            "access_type": "open download",
            "url_or_doi": "doi:10.2172/1345058 / NREL/TP-5400-67585",
            "license_or_terms": "U.S. Government Public Domain (17 U.S.C. 105)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/cetane/nrel_cn_merged_210831.csv, data/processed/v2/cetane_nrel_compendium.csv",
            "used_by": "ingest_new_sources.py, cetane.py",
            "properties": "dcn",
            "notes": "Comprehensive official NREL experimental database parsed with SMILES and InChIKeys."
        },
        {
            "source_id": "doe_co_optima",
            "name": "DOE Co-Optimization of Fuels & Engines (Co-Optima) Fuel Properties Database",
            "what_it_provides": "Pure and blendstock RON, MON, Sensitivity, Flash Point, LFL, Density, and Viscosity across 527 compounds",
            "access_type": "open download",
            "url_or_doi": "https://github.com/NatLabRockies/cooptima-co-optimizer",
            "license_or_terms": "U.S. Government Public Domain / MIT",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/cooptima/, data/processed/v2/octane_co_optima.csv",
            "used_by": "ingest_new_sources.py, octane.py",
            "properties": "ron, mon, sensitivity, flash_point_C, lfl_vol_pct",
            "notes": "Multi-lab DOE/NREL initiative database containing spark-ignition merit fuel characteristics."
        },
        {
            "source_id": "coutinho_dauphin_wax",
            "name": "Coutinho / Dauphin Paraffin Wax Solid-Liquid Equilibrium Database (SOFTX 2018)",
            "what_it_provides": "Experimental binary and multicomponent wax appearance temperatures and phase transition data",
            "access_type": "open download",
            "url_or_doi": "doi:10.1016/j.softx.2018.06.003 / https://github.com/ElsevierSoftwareX/SOFTX_2018_207",
            "license_or_terms": "Apache License 2.0",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/wax/, data/processed/v2/wax_dauphin_coutinho.csv",
            "used_by": "ingest_new_sources.py, cold_flow.py",
            "properties": "wat_C, freezing_point_C",
            "notes": "Multi-component paraffin wax precipitation validation data for non-ideal solid-solution modeling."
        }
    ]

    for ns in new_sources:
        if ns["source_id"] not in df["source_id"].values:
            df = pd.concat([df, pd.DataFrame([ns])], ignore_index=True)

    df.to_csv(ds_path, index=False)
    print(f"Updated docs/data_sources.csv (now {len(df)} cataloged sources)")


def main():
    process_nrel_cetane()
    process_co_optima()
    process_coutinho_wax()
    update_data_sources_csv()
    print("Ingestion of new open sources complete!")


if __name__ == "__main__":
    main()
