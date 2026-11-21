"""Shared helpers for the external-ingest scripts (staging only; never writes to production paths)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# Production helper: md5 of "|".join(sorted InChIKeys), first 12 hex chars. Imported, not copied.
from fuelblend.data.thermoml_sle import system_id  # noqa: E402

STAGING_DIR = REPO / "data" / "interim" / "external_staging"
OUT_DIR = REPO / "results" / "external_benchmarks"
UNIFIED_MEAS = REPO / "data" / "unified" / "measurements.csv"
UNIFIED_COMP = REPO / "data" / "unified" / "compounds.csv"


def mol_info(smiles: str):
    """Return (inchikey, average MolWt, is_hc) from SMILES, or None if RDKit cannot parse it."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    ik = Chem.MolToInchiKey(m)
    is_hc = all(a.GetSymbol() in ("C", "H") for a in m.GetAtoms())
    return ik, float(Descriptors.MolWt(m)), is_hc


def fuel_system_id(source: str, sample: str) -> str:
    """Real-fuel convention used in production loaders: fuel:<source>:<sample>."""
    return f"fuel:{source}:{sample}"
