"""Identity consistency: does a stored InChIKey belong to the stored SMILES? (v4 Task 2)

The key is recomputed from the SMILES with RDKit (``Chem.MolToInchiKey``) and compared with the
stored one.

Comparison rule. The FULL 27-character key is compared, because a swapped key (e.g. the
hexadecane key stored beside a decane SMILES) differs in the connectivity block, but so does
a wrong stereo or protonation layer, and a join on the full key breaks in both cases. The two
are told apart so the caller can triage them:

``key_mismatch``        connectivity block (first 14 characters) differs: a different molecule.
``layer_only``          connectivity block equal, stereo/protonation/other layers differ
                        (``block1_match`` is True). Usually benign for a hydrocarbon property
                        model that ignores stereochemistry; still a failed full-key join.
``nonstandard_key``     SMILES is present but the stored key is not InChIKey-shaped and not a
                        lump marker (e.g. the ``SURR_<smiles>`` placeholders of the surrogate
                        tables): identity cannot be verified or joined.
``unparseable_smiles``  RDKit cannot parse the SMILES.

Skipped (never reported): blank/NaN SMILES, lump rows (key starting with ``LUMP`` -- the real
tables use ``LUMP_<name>`` -- or a true ``is_lumped`` column value).
"""
from __future__ import annotations

import re
from typing import Iterable

import pandas as pd

_KEY_RE = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]")
COLUMNS = ["row_index", "smiles", "stored_inchikey", "recomputed_inchikey",
           "block1_match", "kind"]


def _is_lump(key: str, lumped) -> bool:
    if key.upper().startswith("LUMP"):
        return True
    return lumped is True or str(lumped).strip().lower() == "true"


def inchikey_mismatches(df: pd.DataFrame, smiles_col: str = "smiles", key_col: str = "inchikey",
                        *, id_cols: Iterable[str] = ()) -> pd.DataFrame:
    """Rows of ``df`` whose ``key_col`` is not the RDKit InChIKey of ``smiles_col``.

    Returns columns ``row_index`` (the index label in ``df``), any ``id_cols`` requested,
    ``smiles``, ``stored_inchikey``, ``recomputed_inchikey``, ``block1_match``, ``kind``
    (see module docstring). Lumps and blank SMILES are skipped.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    id_cols = list(id_cols)
    lumped_col = df["is_lumped"] if "is_lumped" in df.columns else pd.Series(False, index=df.index)
    cache: dict[str, str | None] = {}
    rows = []
    id_vals = {c: df[c].tolist() for c in id_cols if c in df.columns}
    for pos, (idx, smi, key, lumped) in enumerate(zip(df.index, df[smiles_col], df[key_col], lumped_col)):
        if not isinstance(smi, str) or not smi.strip():
            continue
        key = "" if key is None or (isinstance(key, float) and pd.isna(key)) else str(key).strip()
        if _is_lump(key, lumped):
            continue
        smi = smi.strip()
        if smi not in cache:
            mol = Chem.MolFromSmiles(smi)
            cache[smi] = Chem.MolToInchiKey(mol) if mol is not None else None
        new = cache[smi]
        if new is None:
            kind, b1 = "unparseable_smiles", False
        elif not _KEY_RE.fullmatch(key):
            kind, b1 = "nonstandard_key", False
        elif new == key:
            continue
        else:
            b1 = new[:14] == key[:14]
            kind = "layer_only" if b1 else "key_mismatch"
        rec = {"row_index": idx, "smiles": smi, "stored_inchikey": key,
               "recomputed_inchikey": new if new is not None else "", "block1_match": b1, "kind": kind}
        for c in id_cols:
            rec[c] = id_vals[c][pos] if c in id_vals else None
        rows.append(rec)
    cols = ["row_index", *id_cols, *COLUMNS[1:]]
    return pd.DataFrame(rows, columns=cols)
