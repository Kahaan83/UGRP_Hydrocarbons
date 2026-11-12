"""Loader for transcribed measured n-alkane mixture flash points (Li et al. 2014; Ha and Lee 2016).

Reads data/raw/flash_lab_li2014/li2014_flash_mixtures.csv (one row per component per point) and returns one row per measured
point in the layout of results/v4/F/flash_lab_loso.csv (system_id, n_comp, measured_C, method, source, comp_smiles, x_mole).
Mole fractions are renormalised to sum to 1 (as printed they sum to 0.99-1.01 in a few rows); the printed values stay in the raw CSV.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "flash_lab_li2014" / "li2014_flash_mixtures.csv"


def _inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    return inchi.MolToInchiKey(mol)


def system_id_from_smiles(smiles) -> str:
    return hashlib.md5("|".join(sorted(_inchikey(s) for s in smiles)).encode()).hexdigest()[:12]


def load_li2014(path: Path | str = RAW) -> pd.DataFrame:
    raw = pd.read_csv(path)
    out = []
    for label, g in raw.groupby("system_label", sort=False):
        x = g.x_mole.to_numpy(float)
        smi = list(g.component_smiles)
        out.append(dict(system_id=system_id_from_smiles(smi), n_comp=len(g), measured_C=float(g.measured_C.iloc[0]),
                        method=g.method.iloc[0], source="li2014" if "Li2014" in label else "ha2016",
                        comp_smiles=json.dumps(smi), x_mole=json.dumps([round(float(v), 6) for v in x / x.sum()]),
                        system_label=label, doi=g.doi.iloc[0], table=g.table.iloc[0], page=int(g.page.iloc[0])))
    return pd.DataFrame(out)
