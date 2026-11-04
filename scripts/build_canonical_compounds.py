"""
scripts/build_canonical_compounds.py - Standardize identifiers for all compounds (T1.3).

Standardizes all compound identifiers across the repository:
- Canonical SMILES (isomericSmiles=True)
- Standard InChI and 27-character InChIKey via RDKit
- CAS Registry Number where available
- Compound common / systematic name where available
- Chemical formula, molecular weight, carbon count, and hydrocarbon classification
- Preserves existing Antoine and LFL parameters from data/unified/compounds.csv

Outputs:
- data/compounds_canonical.csv
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi, rdMolDescriptors

# Suppress RDKit stereo warnings during InChI generation
RDLogger.DisableLog("rdApp.*")

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
from common import is_hydrocarbon

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger("build_canonical_compounds")

INCHIKEY_REGEX = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CAS_REGEX = re.compile(r"^[1-9]\d{1,6}-\d{2}-\d$")
CAS_SEARCH_REGEX = re.compile(r"\b[1-9]\d{1,6}-\d{2}-\d\b")

VALID_CLASSES = {
    "n_paraffin",
    "iso_paraffin",
    "cycloalkane",
    "aromatic",
    "alkene_alkyne",
    "oxygenate",
    "other",
}

CANONICAL_COLUMNS = [
    "inchikey",
    "smiles_canonical",
    "inchi",
    "cas",
    "name",
    "formula",
    "mw",
    "is_hc",
    "cls",
    "n_carbon",
    "antoine_A",
    "antoine_B",
    "antoine_C",
    "antoine_Tlo",
    "antoine_Thi",
    "antoine_source",
    "antoine_covers_window",
    "antoine_unit",
    "lfl_vol_pct",
    "lfl_T_C",
    "lfl_source",
]


def classify_compound(m: Chem.Mol) -> str:
    """Classify molecule into standard fuel categories."""
    if m is None:
        return "other"
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    if is_hydrocarbon(m):
        if any(a.GetIsAromatic() for a in m.GetAtoms()):
            return "aromatic"
        if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
            return "alkene_alkyne"
        if m.GetRingInfo().NumRings() > 0:
            return "cycloalkane"
        branched = any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms())
        return "iso_paraffin" if branched else "n_paraffin"
    if "O" in syms and all(s in ("C", "H", "O") for s in syms):
        return "oxygenate"
    return "other"


def load_enrichment_maps(base_dir: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Build CAS and Name lookup maps from local repository sources."""
    cas_map: Dict[str, str] = {}
    name_map: Dict[str, str] = {}

    # 1. compounds.csv (curated root reference table)
    root_csv = base_dir / "compounds.csv"
    if root_csv.exists():
        df_root = pd.read_csv(root_csv)
        for _, r in df_root.iterrows():
            m = Chem.MolFromSmiles(str(r["smiles"]))
            if m:
                ik = inchi.MolToInchiKey(m)
                if pd.notna(r.get("cas")) and str(r["cas"]).strip():
                    cas_val = str(r["cas"]).strip()
                    if CAS_REGEX.match(cas_val):
                        cas_map[ik] = cas_val
                if pd.notna(r.get("name")) and str(r["name"]).strip():
                    name_map[ik] = str(r["name"]).strip()

    # 2. data/processed/critical_constants.csv (high-coverage CAS)
    crit_csv = base_dir / "data" / "processed" / "critical_constants.csv"
    if crit_csv.exists():
        df_crit = pd.read_csv(crit_csv)
        for _, r in df_crit.iterrows():
            ik = str(r["inchikey"])
            if ik not in cas_map and pd.notna(r.get("cas")) and str(r["cas"]).strip():
                cas_val = str(r["cas"]).strip()
                if CAS_REGEX.match(cas_val):
                    cas_map[ik] = cas_val

    # 3. data/antoine_sources.csv
    ant_csv = base_dir / "data" / "antoine_sources.csv"
    if ant_csv.exists():
        df_ant = pd.read_csv(ant_csv)
        for _, r in df_ant.iterrows():
            ik = str(r["inchikey"])
            if ik not in cas_map and pd.notna(r.get("cas")) and str(r["cas"]).strip():
                cas_val = str(r["cas"]).strip()
                if CAS_REGEX.match(cas_val):
                    cas_map[ik] = cas_val
            if ik not in name_map and pd.notna(r.get("name")) and str(r["name"]).strip():
                name_map[ik] = str(r["name"]).strip()

    # 4. data/raw/sooting/ysi_compounds.csv
    ysi_csv = base_dir / "data" / "raw" / "sooting" / "ysi_compounds.csv"
    if ysi_csv.exists():
        df_ysi = pd.read_csv(ysi_csv)
        for _, r in df_ysi.iterrows():
            ik = str(r["inchikey"])
            if ik not in cas_map and pd.notna(r.get("cas")) and str(r["cas"]).strip():
                cas_val = str(r["cas"]).strip()
                if CAS_REGEX.match(cas_val):
                    cas_map[ik] = cas_val
            if ik not in name_map and pd.notna(r.get("name")) and str(r["name"]).strip():
                name_map[ik] = str(r["name"]).strip()

    # 5. data/processed/ait_compounds.csv
    ait_csv = base_dir / "data" / "processed" / "ait_compounds.csv"
    if ait_csv.exists():
        df_ait = pd.read_csv(ait_csv)
        for _, r in df_ait.iterrows():
            ik = str(r["inchikey"])
            if ik not in cas_map and pd.notna(r.get("cas")) and str(r["cas"]).strip():
                cas_val = str(r["cas"]).strip()
                if CAS_REGEX.match(cas_val):
                    cas_map[ik] = cas_val

    # 6. data/raw/chew2024/si_viscosity_dataset.csv
    chew_csv = base_dir / "data" / "raw" / "chew2024" / "si_viscosity_dataset.csv"
    if chew_csv.exists():
        df_chew = pd.read_csv(chew_csv, usecols=["Name", "CANON_SMILES"])
        for _, r in df_chew.iterrows():
            m = Chem.MolFromSmiles(str(r["CANON_SMILES"]))
            if m:
                ik = inchi.MolToInchiKey(m)
                if ik not in name_map and pd.notna(r["Name"]) and str(r["Name"]).strip():
                    name_map[ik] = str(r["Name"]).strip()

    # 7. data/processed/thermoml_hc.csv
    thermoml_csv = base_dir / "data" / "processed" / "thermoml_hc.csv"
    if thermoml_csv.exists():
        df_t = pd.read_csv(thermoml_csv, usecols=["inchis", "names"])
        for _, r in df_t.iterrows():
            try:
                inchis = json.loads(r["inchis"])
                names = json.loads(r["names"])
                for ic, nm in zip(inchis, names):
                    if nm and str(nm).strip():
                        ik = inchi.InchiToInchiKey(ic)
                        if ik and ik not in name_map:
                            name_map[ik] = str(nm).strip()
            except Exception:
                pass

    # 8. data/raw/pubchem/json/ (RecordTitle and CAS regex fallback)
    cid_map_file = base_dir / "data" / "raw" / "pubchem" / "inchikey_cid_map.json"
    if cid_map_file.exists():
        with open(cid_map_file, encoding="utf-8") as f:
            ik_cid = json.load(f)

        for ik, cid in ik_cid.items():
            jp = base_dir / "data" / "raw" / "pubchem" / "json" / f"{cid}.json"
            if jp.exists():
                try:
                    with open(jp, encoding="utf-8") as f:
                        data = json.load(f)
                    if ik not in name_map:
                        title = data.get("Record", {}).get("RecordTitle")
                        if title and str(title).strip():
                            name_map[ik] = str(title).strip()
                    if ik not in cas_map:
                        text = json.dumps(data)
                        matches = CAS_SEARCH_REGEX.findall(text)
                        if matches:
                            cas_map[ik] = matches[0]
                except Exception:
                    pass

    return cas_map, name_map


def build_canonical_compounds(
    unified_csv: Path,
    root_csv: Path,
    output_csv: Path,
) -> Tuple[pd.DataFrame, List[Dict]]:
    """Build standardized, deduplicated compounds dataframe."""
    base_dir = unified_csv.resolve().parent.parent.parent
    cas_map, name_map = load_enrichment_maps(base_dir)

    df_unified = pd.read_csv(unified_csv)
    df_root = pd.read_csv(root_csv) if root_csv.exists() else pd.DataFrame()

    registry: Dict[str, Dict] = {}
    merged_duplicates: List[Dict] = []
    order: List[str] = []

    # Process data/unified/compounds.csv
    for _, row in df_unified.iterrows():
        s = row.get("smiles")
        m = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
        if m is None and pd.notna(row.get("inchi")):
            m = Chem.MolFromInchi(str(row["inchi"]))

        if m is None:
            raise ValueError(f"Could not parse molecule for row: {row.to_dict()}")

        smiles_can = Chem.MolToSmiles(m, isomericSmiles=True)
        ik = inchi.MolToInchiKey(m)
        ic = inchi.MolToInchi(m)

        if not INCHIKEY_REGEX.match(ik):
            raise ValueError(f"Invalid InChIKey generated: '{ik}' for SMILES '{smiles_can}'")

        rec = {
            "inchikey": ik,
            "smiles_canonical": smiles_can,
            "inchi": ic,
            "cas": cas_map.get(ik, ""),
            "name": name_map.get(ik, ""),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "mw": round(float(Descriptors.MolWt(m)), 3),
            "is_hc": bool(is_hydrocarbon(m)),
            "cls": classify_compound(m),
            "n_carbon": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C"),
            "antoine_A": row.get("antoine_A"),
            "antoine_B": row.get("antoine_B"),
            "antoine_C": row.get("antoine_C"),
            "antoine_Tlo": row.get("antoine_Tlo"),
            "antoine_Thi": row.get("antoine_Thi"),
            "antoine_source": row.get("antoine_source"),
            "antoine_covers_window": row.get("antoine_covers_window"),
            "antoine_unit": row.get("antoine_unit"),
            "lfl_vol_pct": row.get("lfl_vol_pct"),
            "lfl_T_C": row.get("lfl_T_C"),
            "lfl_source": row.get("lfl_source"),
        }

        if ik in registry:
            # Duplicate within unified compounds
            merged_duplicates.append({
                "inchikey": ik,
                "source": "unified/compounds.csv",
                "source_name": "",
                "source_smiles": s,
                "canonical_smiles": smiles_can,
            })
            # Merge non-empty fields
            for k, v in rec.items():
                if (pd.isna(registry[ik][k]) or registry[ik][k] == "") and pd.notna(v) and v != "":
                    registry[ik][k] = v
        else:
            registry[ik] = rec
            order.append(ik)

    # Process compounds.csv (root)
    if not df_root.empty:
        for _, row in df_root.iterrows():
            s = row.get("smiles")
            m = Chem.MolFromSmiles(str(s)) if pd.notna(s) else None
            if m is None:
                continue

            smiles_can = Chem.MolToSmiles(m, isomericSmiles=True)
            ik = inchi.MolToInchiKey(m)

            if ik in registry:
                merged_duplicates.append({
                    "inchikey": ik,
                    "source": "compounds.csv",
                    "source_name": row.get("name", ""),
                    "source_smiles": s,
                    "canonical_smiles": smiles_can,
                })
                # Prioritize curated name and CAS from compounds.csv
                if pd.notna(row.get("name")) and str(row["name"]).strip():
                    registry[ik]["name"] = str(row["name"]).strip()
                if pd.notna(row.get("cas")) and str(row["cas"]).strip():
                    cas_val = str(row["cas"]).strip()
                    if CAS_REGEX.match(cas_val):
                        registry[ik]["cas"] = cas_val
            else:
                # Compound present in root compounds.csv but not in unified
                ic = inchi.MolToInchi(m)
                rec = {
                    "inchikey": ik,
                    "smiles_canonical": smiles_can,
                    "inchi": ic,
                    "cas": str(row.get("cas", "")).strip() or cas_map.get(ik, ""),
                    "name": str(row.get("name", "")).strip() or name_map.get(ik, ""),
                    "formula": rdMolDescriptors.CalcMolFormula(m),
                    "mw": round(float(Descriptors.MolWt(m)), 3),
                    "is_hc": bool(is_hydrocarbon(m)),
                    "cls": classify_compound(m),
                    "n_carbon": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C"),
                    "antoine_A": np.nan,
                    "antoine_B": np.nan,
                    "antoine_C": np.nan,
                    "antoine_Tlo": np.nan,
                    "antoine_Thi": np.nan,
                    "antoine_source": np.nan,
                    "antoine_covers_window": False,
                    "antoine_unit": np.nan,
                    "lfl_vol_pct": np.nan,
                    "lfl_T_C": np.nan,
                    "lfl_source": np.nan,
                }
                registry[ik] = rec
                order.append(ik)

    records = [registry[ik] for ik in order]
    df_canon = pd.DataFrame(records)[CANONICAL_COLUMNS]

    # Ensure output directory exists and write CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df_canon.to_csv(output_csv, index=False)
    logger.info("Wrote %d canonical compounds to %s", len(df_canon), output_csv)

    return df_canon, merged_duplicates


def main():
    parser = argparse.ArgumentParser(description="Standardize compound identifiers (T1.3)")
    parser.add_argument(
        "--unified-csv",
        type=Path,
        default=BASE_DIR / "data" / "unified" / "compounds.csv",
        help="Path to data/unified/compounds.csv",
    )
    parser.add_argument(
        "--root-csv",
        type=Path,
        default=BASE_DIR / "compounds.csv",
        help="Path to root compounds.csv",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=BASE_DIR / "data" / "compounds_canonical.csv",
        help="Path for output canonical compounds CSV",
    )
    args = parser.parse_args()

    df_canon, duplicates = build_canonical_compounds(
        unified_csv=args.unified_csv,
        root_csv=args.root_csv,
        output_csv=args.output_csv,
    )

    print("\n--- Canonical Compounds Summary (T1.3) ---")
    print(f"Total Unique Compounds : {len(df_canon)}")
    print(f"Duplicates Merged      : {len(duplicates)}")
    for d in duplicates[:5]:
        print(f"  - [{d['source']}] {d['source_name'] or d['source_smiles']} -> InChIKey {d['inchikey']}")
    if len(duplicates) > 5:
        print(f"  ... and {len(duplicates) - 5} more duplicates logged.")

    print("\nClasses Distribution:")
    for cls, count in df_canon["cls"].value_counts().items():
        print(f"  {cls:<15}: {count}")

    cas_count = df_canon["cas"].replace("", np.nan).notna().sum()
    name_count = df_canon["name"].replace("", np.nan).notna().sum()
    print("\nMetadata Coverage:")
    print(f"  CAS Available : {cas_count} / {len(df_canon)} ({cas_count / len(df_canon) * 100:.1f}%)")
    print(f"  Name Available: {name_count} / {len(df_canon)} ({name_count / len(df_canon) * 100:.1f}%)")
    print(f"  Antoine Sets  : {df_canon['antoine_source'].notna().sum()} / {len(df_canon)}")
    print(f"  LFL Records   : {df_canon['lfl_source'].notna().sum()} / {len(df_canon)}")


if __name__ == "__main__":
    main()
