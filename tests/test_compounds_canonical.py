"""Unit tests for standardized canonical compound registry (T1.3)."""
from pathlib import Path
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
CANON_CSV = ROOT / "data" / "compounds_canonical.csv"
UNIFIED_CSV = ROOT / "data" / "unified" / "compounds.csv"


def test_compounds_canonical_schema():
    assert CANON_CSV.exists(), f"Missing {CANON_CSV}"
    df = pd.read_csv(CANON_CSV)

    # 1. Total row count matches unified
    df_uni = pd.read_csv(UNIFIED_CSV)
    assert len(df) == len(df_uni)

    # 2. InChIKey validity & uniqueness
    assert df["inchikey"].nunique() == len(df), "InChIKeys must be unique"
    for k in df["inchikey"]:
        assert isinstance(k, str) and len(k) == 27 and k[14] == "-" and k[25] == "-", f"Invalid InChIKey: {k}"

    # 3. Canonical SMILES validity with RDKit
    for s in df["smiles_canonical"].dropna():
        if s.strip():
            m = Chem.MolFromSmiles(s)
            assert m is not None, f"Unparseable SMILES: {s}"

    # 4. CAS coverage
    cas_count = df["cas"].notna().sum()
    assert cas_count >= 900, f"Expected at least 900 compounds with CAS numbers, found {cas_count}"
