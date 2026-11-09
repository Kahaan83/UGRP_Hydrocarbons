"""v4 Task 2: InChIKey vs SMILES identity consistency.

The NASA / CRC / surrogate composition tables pair InChIKeys with the wrong SMILES.
These tests run on the REAL tables plus a small synthetic unit test of the checker.
"""
from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

pytest.importorskip("rdkit")

from fuelblend.data.identities import inchikey_mismatches  # noqa: E402

# real keys (computed with RDKit): heptane / hexadecane
HEPTANE = ("CCCCCCC", "IMNFDUFMRHMDMM-UHFFFAOYSA-N")
HEXADECANE = ("CCCCCCCCCCCCCCCC", "DCAYPVUWAIABOU-UHFFFAOYSA-N")


def test_synthetic_unit():
    df = pd.DataFrame({
        "smiles": [HEPTANE[0], HEXADECANE[0], "", None, "CCCCCCC", "CCCCCCC", "not a smiles"],
        "inchikey": [HEPTANE[1], HEPTANE[1],            # row 1: deliberate mismatch
                     "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",     # blank SMILES -> skipped
                     "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N",     # NaN SMILES -> skipped
                     "LUMP_UNRESOLVED_TAIL",            # lump key -> skipped
                     "IMNFDUFMRHMDMM-UHFFFAOYSB-N",     # same connectivity block, other layer
                     "IMNFDUFMRHMDMM-UHFFFAOYSA-N"],    # unparseable SMILES
    })
    out = inchikey_mismatches(df)
    assert list(out["row_index"]) == [1, 5, 6]
    r = out.set_index("row_index")
    assert r.loc[1, "kind"] == "key_mismatch" and not r.loc[1, "block1_match"]
    assert r.loc[1, "recomputed_inchikey"] == HEXADECANE[1]
    assert r.loc[5, "kind"] == "layer_only" and r.loc[5, "block1_match"]
    assert r.loc[6, "kind"] == "unparseable_smiles"
    assert inchikey_mismatches(df.iloc[[0]]).empty        # matching key is not flagged
    # lump flag column is honoured too
    d2 = pd.DataFrame({"smiles": ["CCCCCCC"], "inchikey": [HEXADECANE[1]], "is_lumped": [True]})
    assert inchikey_mismatches(d2).empty


def test_nasa_diesel_identity_mismatch_detected():
    df = pd.read_csv(ROOT / "data/real_fuels/v3/compositions_v3.csv")
    nasa = df[df["system_id"].astype(str).str.contains("NASA_LFP")]
    assert len(nasa) > 0
    out = inchikey_mismatches(nasa)
    assert len(out) >= 3
    assert (out["kind"] == "key_mismatch").sum() >= 3
    # the specific swaps reported by the Task 1 audit
    pairs = set(zip(out["smiles"], out["stored_inchikey"]))
    assert ("CCCCCCCCCC", "DCAYPVUWAIABOU-UHFFFAOYSA-N") in pairs        # decane carries hexadecane key
    assert ("CCCCCCCCCCCCC", "NHTMVDHEPJAVLT-UHFFFAOYSA-N") in pairs     # tridecane carries isooctane key
    assert ("CCCCCCCCCCCCCCCC", "IMNFDUFMRHMDMM-UHFFFAOYSA-N") in pairs   # hexadecane carries heptane key


def test_canonical_registry_clean():
    reg = ROOT / "data/unified/compounds.csv"          # registry used by src/fuelblend (untracked, rebuildable)
    if not reg.exists():
        pytest.skip("data/unified/compounds.csv not present (gitignored, rebuildable)")
    out = inchikey_mismatches(pd.read_csv(reg))
    assert len(out) == 0, out.head().to_string()
    anch = pd.read_csv(ROOT / "data/processed/v2/compounds_anchored_vp.csv")   # tracked copy used by v3 flash code
    assert len(inchikey_mismatches(anch)) == 0
