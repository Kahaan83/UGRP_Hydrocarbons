"""
tests/test_compounds_canonical.py - Test suite for standardized canonical compounds (T1.3).

Verifies:
1. Every row has a valid 27-character InChIKey matching regex ^[A-Z]{14}-[A-Z]{10}-[A-Z]$.
2. All InChIKeys in data/compounds_canonical.csv are unique (no duplicates).
3. Every compound from data/unified/compounds.csv is mapped.
4. SMILES strings are valid RDKit molecules and recompute to the same InChIKey.
5. All root compounds.csv compounds are mapped.
6. Required standard columns and preserved Antoine/LFL columns exist.
7. Class labels belong to the standard set.
"""

from __future__ import annotations

from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

# Suppress RDKit stereo warnings
RDLogger.DisableLog("rdApp.*")

BASE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE_DIR))
from common import is_hydrocarbon

INCHIKEY_REGEX = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

EXPECTED_COLUMNS = [
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

VALID_CLASSES = {
    "n_paraffin",
    "iso_paraffin",
    "cycloalkane",
    "aromatic",
    "alkene_alkyne",
    "oxygenate",
    "other",
}


def load_canonical_df() -> pd.DataFrame:
    canon_path = BASE_DIR / "data" / "compounds_canonical.csv"
    assert canon_path.exists(), f"File not found: {canon_path}"
    return pd.read_csv(canon_path)


def test_schema_and_columns():
    """Verify all expected columns are present."""
    df = load_canonical_df()
    for col in EXPECTED_COLUMNS:
        assert col in df.columns, f"Missing required column: {col}"
    assert len(df) >= 1038, f"Expected at least 1038 compounds, found {len(df)}"


def test_valid_inchikeys():
    """Requirement 1: Every row has a valid 27-character InChIKey."""
    df = load_canonical_df()
    for idx, key in enumerate(df["inchikey"]):
        assert isinstance(key, str), f"Row {idx}: InChIKey is not a string ({key})"
        assert len(key) == 27, f"Row {idx}: InChIKey length is {len(key)} != 27 ({key})"
        assert INCHIKEY_REGEX.match(key), f"Row {idx}: InChIKey format invalid ({key})"


def test_unique_inchikeys():
    """Requirement 2: All InChIKeys in data/compounds_canonical.csv are unique."""
    df = load_canonical_df()
    keys = df["inchikey"].tolist()
    unique_keys = set(keys)
    assert len(keys) == len(unique_keys), (
        f"Found {len(keys) - len(unique_keys)} duplicate InChIKeys in canonical table!"
    )


def test_all_unified_compounds_mapped():
    """Requirement 3: Every compound from data/unified/compounds.csv is mapped."""
    df_canon = load_canonical_df()
    canon_keys = set(df_canon["inchikey"])

    unified_path = BASE_DIR / "data" / "unified" / "compounds.csv"
    assert unified_path.exists(), f"Unified compounds missing: {unified_path}"
    df_uni = pd.read_csv(unified_path)

    missing = set(df_uni["inchikey"]) - canon_keys
    assert not missing, f"Missing {len(missing)} compounds from unified table: {missing}"


def test_all_root_compounds_mapped():
    """Requirement: Every compound from compounds.csv (root) is mapped."""
    df_canon = load_canonical_df()
    canon_keys = set(df_canon["inchikey"])

    root_path = BASE_DIR / "compounds.csv"
    if root_path.exists():
        df_root = pd.read_csv(root_path)
        for _, row in df_root.iterrows():
            m = Chem.MolFromSmiles(str(row["smiles"]))
            assert m is not None, f"Invalid SMILES in root compounds.csv: {row['smiles']}"
            ik = inchi.MolToInchiKey(m)
            assert ik in canon_keys, (
                f"Root compound '{row.get('name', row['smiles'])}' with key {ik} not mapped!"
            )


def test_valid_rdkit_smiles():
    """Requirement 4: Smiles are valid RDKit molecules and recompute to the same InChIKey."""
    df = load_canonical_df()
    for idx, row in df.iterrows():
        s = row["smiles_canonical"]
        assert isinstance(s, str) and len(s) > 0, f"Row {idx}: empty SMILES"
        m = Chem.MolFromSmiles(s)
        assert m is not None, f"Row {idx}: RDKit could not parse SMILES '{s}'"

        # Recomputed InChIKey must match
        recomputed_key = inchi.MolToInchiKey(m)
        assert recomputed_key == row["inchikey"], (
            f"Row {idx}: Key mismatch for SMILES '{s}'. "
            f"Expected {row['inchikey']}, got {recomputed_key}"
        )


def test_valid_classes():
    """Verify class labels belong to the standard set and is_hc is accurate."""
    df = load_canonical_df()
    for idx, row in df.iterrows():
        c = row["cls"]
        assert c in VALID_CLASSES, f"Row {idx}: Invalid class '{c}'"

        m = Chem.MolFromSmiles(row["smiles_canonical"])
        expected_hc = is_hydrocarbon(m)
        assert bool(row["is_hc"]) == expected_hc, (
            f"Row {idx}: is_hc flag mismatch ({row['is_hc']} vs {expected_hc})"
        )


def test_preserved_physics_parameters():
    """Verify Antoine and LFL columns from unified compounds are preserved."""
    df_canon = load_canonical_df().set_index("inchikey")
    df_uni = pd.read_csv(BASE_DIR / "data" / "unified" / "compounds.csv").set_index("inchikey")

    # Check non-null Antoine coverage is fully preserved
    uni_antoine = df_uni["antoine_source"].dropna()
    for ik, src in uni_antoine.items():
        assert ik in df_canon.index, f"Key {ik} missing from canonical index"
        assert df_canon.loc[ik, "antoine_source"] == src, (
            f"Antoine source mismatch for {ik}: {df_canon.loc[ik, 'antoine_source']} vs {src}"
        )

    # Check non-null LFL coverage is fully preserved
    uni_lfl = df_uni["lfl_vol_pct"].dropna()
    for ik, val in uni_lfl.items():
        assert ik in df_canon.index, f"Key {ik} missing from canonical index"
        canon_val = df_canon.loc[ik, "lfl_vol_pct"]
        assert np.isclose(canon_val, val), (
            f"LFL mismatch for {ik}: {canon_val} vs {val}"
        )


if __name__ == "__main__":
    print("Running test_schema_and_columns()...")
    test_schema_and_columns()
    print("Running test_valid_inchikeys()...")
    test_valid_inchikeys()
    print("Running test_unique_inchikeys()...")
    test_unique_inchikeys()
    print("Running test_all_unified_compounds_mapped()...")
    test_all_unified_compounds_mapped()
    print("Running test_all_root_compounds_mapped()...")
    test_all_root_compounds_mapped()
    print("Running test_valid_rdkit_smiles()...")
    test_valid_rdkit_smiles()
    print("Running test_valid_classes()...")
    test_valid_classes()
    print("Running test_preserved_physics_parameters()...")
    test_preserved_physics_parameters()
    print("\nALL 8 TESTS PASSED SUCCESSFULLY!")
