"""Pure-compound feature extraction (WP4).

One feature row per pure compound at a stated temperature. Vector layout (19 features, < 25 total):

    [mol_wt, heavy_atom_count, n_aromatic_rings, ring_count, formal_charge,
     hbond_donor, hbond_acceptor, rotatable_bonds, labute_asa, tpsa, n_carbon,
     cls_<c> for c in CLASSES,
     logT, T_K, inv_T_K]

where logT = log(T_K), inv_T_K = 1000 / T_K. T_fallback_K is used when T_K is NaN
(flash/freezing pure rows carry the measured value as the temperature).
"""
from __future__ import annotations

import functools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]
DESCRIPTOR = ["mol_wt", "heavy_atom_count", "n_aromatic_rings", "ring_count", "formal_charge",
               "hbond_donor", "hbond_acceptor", "rotatable_bonds", "labute_asa", "tpsa", "n_carbon"]
T_FEATURES = ["logT", "T_K", "inv_T_K"]
FEATURE_NAMES = list(DESCRIPTOR) + ["cls_" + c for c in CLASSES] + list(T_FEATURES)
N_FEATURES = len(FEATURE_NAMES)
CACHE_COLUMNS = ["inchikey", "T_K"] + [n for n in FEATURE_NAMES if n != "T_K"]
DEFAULT_T_K = 298.15


def smiles_to_mol(smiles):
    if not smiles:
        return None
    try:
        return Chem.MolFromSmiles(str(smiles))
    except Exception:
        return None


def _to_float(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


@functools.lru_cache(maxsize=4096)
def descriptor_vector(smiles) -> dict:
    """11 RDKit descriptor values; NaN entries if the SMILES do not parse."""
    m = smiles_to_mol(smiles)
    if m is None:
        return {k: float("nan") for k in DESCRIPTOR}
    return {
        "mol_wt": float(Descriptors.MolWt(m)),
        "heavy_atom_count": float(Descriptors.HeavyAtomCount(m)),
        "n_aromatic_rings": float(Descriptors.NumAromaticRings(m)),
        "ring_count": float(m.GetRingInfo().NumRings()),
        "formal_charge": float(sum(a.GetFormalCharge() for a in m.GetAtoms())),
        "hbond_donor": float(Descriptors.NumHDonors(m)),
        "hbond_acceptor": float(Descriptors.NumHAcceptors(m)),
        "rotatable_bonds": float(Descriptors.NumRotatableBonds(m)),
        "labute_asa": float(Descriptors.LabuteASA(m)),
        "tpsa": float(Descriptors.TPSA(m)),
        "n_carbon": float(sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")),
    }


def class_vector(cls, classes=CLASSES) -> np.ndarray:
    v = np.zeros(len(classes), dtype=float)
    c = str(cls) if cls is not None else ""
    if c in classes:
        v[classes.index(c)] = 1.0
    return v


def pure_feature_row(inchikey, smiles, T_K, cls, classes=CLASSES, T_fallback_K=DEFAULT_T_K) -> np.ndarray:
    """19-dim feature row for one pure compound (inchikey/smiles not strictly needed for math)."""
    d = descriptor_vector(smiles)
    t = _to_float(T_K)
    if math.isnan(t) or math.isinf(t):
        t = float(T_fallback_K)
    row = [float(d[k]) for k in DESCRIPTOR]
    row += [float(x) for x in class_vector(cls, classes)]
    row.append(math.log(t))
    row.append(t)
    row.append(1000.0 / t)
    return np.asarray(row, dtype=float)


def pure_feature_frame(rows, classes=CLASSES, T_fallback_K=DEFAULT_T_K) -> np.ndarray:
    """rows: list of dicts with keys inchikey/smiles/T_K/cls. Returns (n, N_FEATURES) array."""
    out = []
    for r in rows:
        out.append(pure_feature_row(r.get("inchikey"), r.get("smiles"), r.get("T_K"),
                                    r.get("cls"), classes=classes, T_fallback_K=T_fallback_K))
    X = np.asarray(out, dtype=float)
    if X.size == 0:
        return X.reshape(0, N_FEATURES)
    if X.ndim == 1:
        return X.reshape(1, -1)
    return X


def cache_key(inchikey, T_K) -> tuple:
    return (str(inchikey), round(float(T_K), 6))


def default_cache_path() -> Path:
    return Path("data/runs/pure_feature_cache.csv")


def load_cache(path=None) -> pd.DataFrame:
    p = Path(path) if path is not None else default_cache_path()
    if not p.exists():
        return pd.DataFrame(columns=CACHE_COLUMNS)
    return pd.read_csv(p, low_memory=False)


def save_cache(df: pd.DataFrame, path=None) -> None:
    p = Path(path) if path is not None else default_cache_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    df[CACHE_COLUMNS].to_csv(p, index=False)


def lookup_cache(df: pd.DataFrame, inchikey, T_K) -> np.ndarray | None:
    k = cache_key(inchikey, T_K)
    sub = df[(df["inchikey"] == k[0]) & (df["T_K"] == k[1])]
    if len(sub) == 0:
        return None
    return sub.iloc[0, 2:].to_numpy(dtype=float)

def drop_class_columns(X: np.ndarray) -> np.ndarray:
    """Return X without the cls_* one-hot block (ablation: without class one-hot)."""
    if not isinstance(X, np.ndarray):
        X = np.asarray(X, dtype=float)
    return np.asarray(X[:, [i for i, n in enumerate(FEATURE_NAMES) if not n.startswith("cls_")]], dtype=float)
