"""Structure and physical-range checks for data/raw/sle_constants/tm_dhfus_verified.csv."""
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

CSV = Path(__file__).resolve().parents[1] / "data" / "raw" / "sle_constants" / "tm_dhfus_verified.csv"
COLS = ["inchikey", "name", "smiles", "tm_K", "dhfus_kJ_mol", "dhfus_kind", "source_tm", "source_dh", "doi_or_url",
        "coutinho_tm_K", "coutinho_dhfus_kJ_mol", "agrees", "note"]
A_LIST = {"n-decane": "DIOQZVSQGTUSAI-UHFFFAOYSA-N", "n-undecane": "RSJKGSCJYJTIGS-UHFFFAOYSA-N",
          "n-dodecane": "SNRUBQQJIBEYMU-UHFFFAOYSA-N", "n-tridecane": "IIYFAKIEWZDVMP-UHFFFAOYSA-N",
          "n-tetradecane": "BGHCVCJVXZWKCC-UHFFFAOYSA-N", "butylcyclohexane": "GGBJHURWWWLEQH-UHFFFAOYSA-N",
          "trans-decalin": "NNBZCPXTIHJBJL-MGCOHNPYSA-N", "1,2,4-trimethylbenzene": "GWHJZXXIDMPWGX-UHFFFAOYSA-N"}


@pytest.fixture(scope="module")
def df():
    return pd.read_csv(CSV)


def test_parses_with_columns(df):
    assert all(c in df.columns for c in COLS)
    assert len(df) >= 39
    assert df.inchikey.is_unique


def test_inchikeys_valid(df):
    for _, r in df.iterrows():
        assert Chem.MolToInchiKey(Chem.MolFromSmiles(r.smiles)) == r.inchikey, r["name"]


def test_a_list_complete(df):
    sub = df[df.inchikey.isin(A_LIST.values())]
    assert set(sub.inchikey) == set(A_LIST.values())
    assert sub.tm_K.notna().all() and sub.dhfus_kJ_mol.notna().all()
    assert (sub.list == "A").all()


def test_physical_range(df):
    ok = df.dropna(subset=["tm_K"])
    assert ((ok.tm_K > 100) & (ok.tm_K < 550)).all()
    h = df.dropna(subset=["dhfus_kJ_mol"])
    assert ((h.dhfus_kJ_mol > 0.1) & (h.dhfus_kJ_mol < 100)).all()
    t = df.dropna(subset=["dhfus_total_incl_ss_kJ_mol"])
    assert (t.dhfus_total_incl_ss_kJ_mol > t.dhfus_kJ_mol).all()
    assert (df.dropna(subset=["dhfus_kJ_mol"]).dhfus_kind == "fusion only").all()
    assert df[df.tm_K.notna()].source_tm.str.len().gt(10).all()
    assert df[df.dhfus_kJ_mol.notna()].source_dh.str.len().gt(10).all()
