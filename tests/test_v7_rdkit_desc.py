import numpy as np
from fuelblend.features.rdkit_desc import rdkit_2d
from rdkit.Chem import Descriptors


def test_rdkit_2d_shape_and_names():
    df = rdkit_2d(["CCO", "c1ccccc1", "CCCCCC"])
    assert df.shape == (3, len(Descriptors.descList))
    assert list(df.columns) == [n for n, _ in Descriptors.descList]


def test_invalid_smiles_row_nan():
    df = rdkit_2d(["CCO", "not_a_smiles"])
    assert df.iloc[1].isna().all() and df.iloc[0].notna().any()
