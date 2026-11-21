"""RDKit 2D descriptors (all Descriptors.descList). NaN kept; invalid SMILES give an all-NaN row. Constant-column dropping is done inside each fold by the caller."""
from __future__ import annotations
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")


def rdkit_2d(smiles: list[str]) -> pd.DataFrame:
    names = [n for n, _ in Descriptors.descList]
    rows = []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is None:
            rows.append([np.nan] * len(names)); continue
        r = []
        for _, f in Descriptors.descList:
            try: r.append(float(f(m)))
            except Exception: r.append(np.nan)
        rows.append(r)
    df = pd.DataFrame(rows, columns=names)
    return df.replace([np.inf, -np.inf], np.nan)
