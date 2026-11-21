"""kaust_registry.py - curated name->SMILES table for KAUST Table 3 and the blend components.

Stereochemistry is deliberately omitted (the source names are mostly cis/trans-unspecified), so the
InChIKeys are the flat (connectivity-only) keys. Every InChIKey is computed by RDKit, never typed.
OCR typos in the extracted Table 3 names are corrected here (see NAME_FIXES).
"""
from __future__ import annotations

import json

import pandas as pd

from _common import STAGING_DIR, UNIFIED_MEAS, REPO, mol_info

# Table 3 serial number -> SMILES (128 entries).
TABLE3_SMILES = {
    1: "CCC", 2: "CCCC", 3: "CCCCC", 4: "CC(C)CCC", 5: "CCC(C)CC", 6: "CC(C)(C)CC", 7: "CC(C)C(C)C",
    8: "CCC(CC)CC", 9: "CC(C)(C)CCC", 10: "CC(C)C(C)CC", 11: "CC(C)CC(C)C", 12: "CCC(C)(C)CC",
    13: "CC(C)(C)CCCC", 14: "CC(C)C(C)CCC", 15: "CC(C)CC(C)CC", 16: "CCC(C)(C)CCC", 17: "CCC(C)C(C)CC",
    18: "CC(C)C(CC)CC", 19: "CCC(C)(CC)CC", 20: "CC(C)(C)C(C)CC", 21: "CC(C)CC(C)(C)C",
    22: "CC(C)C(C)(C)CC", 23: "CC(C)C(C)C(C)C", 24: "CCC(CC)(CC)CC", 25: "CC(C)(C)C(CC)CC",
    26: "CC(C)C(CC)C(C)C", 27: "CCC(C)(C)C(C)(C)C", 28: "CCC(C)(C)CC(C)CC", 29: "CCCC(C)(C)C(C)(C)C",
    30: "C=C", 31: "CC=C", 32: "C=CCC", 33: "CC=CC", 34: "C=CCCC", 35: "C=C(C)CC", 36: "CC=C(C)C",
    37: "C=CCCCC", 38: "C=C(C)CCC", 39: "C=CC(C)CC", 40: "C=CCC(C)C", 41: "CC(C)=CCC", 42: "CC=CC(C)C",
    43: "C=C(CC)CC", 44: "C=C(C)C(C)C", 45: "C=CC(C)(C)C", 46: "CC(C)=C(C)C", 47: "C=CCCCCC",
    48: "C=C(C)CCCC", 49: "C=CC(C)CCC", 50: "C=CCC(C)CC", 51: "C=CCCC(C)C", 52: "CC(C)=CCCC",
    53: "CC=C(C)CCC", 54: "C=CC(CC)CC", 55: "CC=C(CC)CC", 56: "C=C(C)C(C)CC", 57: "C=C(C)CC(C)C",
    58: "C=CC(C)(C)CC", 59: "C=CC(C)C(C)C", 60: "C=CCC(C)(C)C", 61: "CCC(C)=C(C)C", 62: "CC(C)=CC(C)C",
    63: "CC=C(C)C(C)C", 64: "CC=CC(C)(C)C", 65: "C=C(CC)C(C)C", 66: "C=C(C)C(C)(C)C",
    67: "C=C(C)CCCCC", 68: "C=C(C)C(C)CCC", 69: "CCCC(C)=C(C)C", 70: "CCC=CC(C)(C)C",
    71: "C=C(C)C(C)(C)CC", 72: "C=C(C)CC(C)(C)C", 73: "CC(C)=CC(C)(C)C", 74: "C=CC(C)=C",
    75: "C=CCCC=C", 76: "C1=CCCC1", 77: "CC1=CCCC1", 78: "CCC1=CCCC1", 79: "CCC1C=CCC1", 80: "C1=CCCCC1",
    81: "CC1=CCCCC1", 82: "CCC1=CCCCC1", 83: "C1CCCC1", 84: "CC1CCCC1", 85: "CCC1CCCC1",
    86: "CC1(C)CCCC1", 87: "CC1CCC(C)C1", 88: "CCCC1CCCC1", 89: "CC(C)C1CCCC1", 90: "CCC1CCC(C)C1",
    91: "CC1CCC(C)(C)C1", 92: "C1CCCCC1", 93: "CC1CCCCC1", 94: "CCC1CCCCC1", 95: "CC1(C)CCCCC1",
    96: "CC1CCCCC1C", 97: "CC1CCCC(C)C1", 98: "CC1CCC(C)CC1", 99: "CC(C)C1CCCCC1", 100: "c1ccccc1",
    101: "Cc1ccccc1", 102: "CCc1ccccc1", 103: "Cc1ccccc1C", 104: "Cc1cccc(C)c1", 105: "Cc1ccc(C)cc1",
    106: "CCCc1ccccc1", 107: "CC(C)c1ccccc1", 108: "CCc1ccccc1C", 109: "CCc1cccc(C)c1",
    110: "CCc1ccc(C)cc1", 111: "Cc1cccc(C)c1C", 112: "Cc1cc(C)cc(C)c1", 113: "CCCCc1ccccc1",
    114: "CC(C)Cc1ccccc1", 115: "CCC(C)c1ccccc1", 116: "CC(C)(C)c1ccccc1", 117: "CCCc1ccccc1C",
    118: "CCCc1cccc(C)c1", 119: "CC(C)c1ccccc1C", 120: "Cc1ccc(C(C)C)cc1", 121: "CCc1cccc(CC)c1",
    122: "CCc1ccc(CC)cc1", 123: "CCc1cccc(C)c1C", 124: "CCc1ccc(C)cc1C", 125: "CCc1cc(C)cc(C)c1",
    126: "CCc1cc(C)ccc1C", 127: "Cc1ccc(C)c(C)c1C", 128: "CCO",
}

# Blend component name (as in the parsed Table 4/5) -> (SMILES, Table 3 s_no or None, +0 row name or None)
# Pure RON/MON are never typed: they come from the +0 row (Table 4) when one exists, otherwise Table 3.
# n-heptane (0/0) is not in Table 3: PRF definition (primary reference fuel scale zero point).
COMPONENTS = {
    "n-heptane": {"smiles": "CCCCCCC", "t3": None, "zero_row": None, "prf_definition": (0.0, 0.0)},
    "i-octane": {"smiles": "CC(C)CC(C)(C)C", "t3": 21, "zero_row": None},
    "toluene": {"smiles": "Cc1ccccc1", "t3": 101, "zero_row": None},
    "ethanol": {"smiles": "CCO", "t3": 128, "zero_row": None},
    "1,2,4-trimethylbenzene": {"smiles": "Cc1ccc(C)c(C)c1", "t3": None, "zero_row": "TMB+0"},
    "cyclopentane": {"smiles": "C1CCCC1", "t3": 83, "zero_row": "CP+0"},
    "1-hexene": {"smiles": "C=CCCCC", "t3": 37, "zero_row": "Hex+0"},
    # Commercial diisobutylene is a mixture of the 1- and 2-isomers; no pure DIB row exists in the
    # tables, so the main isomer (2,4,4-trimethyl-1-pentene, Table 3 row 72) is used as a proxy.
    "diisobutylene": {"smiles": "C=C(C)CC(C)(C)C", "t3": 72, "zero_row": None, "proxy": True},
}


def component_inchikeys() -> dict:
    return {k: mol_info(v["smiles"])[0] for k, v in COMPONENTS.items()}


def load_tables():
    t3 = pd.read_csv(STAGING_DIR / "kaust_jameel_table3_pure.csv")
    t4 = pd.read_csv(STAGING_DIR / "kaust_jameel_table4_blends.csv")
    t5 = pd.read_csv(STAGING_DIR / "kaust_jameel_table5_face.csv")
    return t3, t4, t5


def pure_props(t3=None, t4=None) -> dict:
    """Pure RON/MON per component, taken from measured rows (see COMPONENTS). Returns name -> dict."""
    if t3 is None or t4 is None:
        t3, t4, _ = load_tables()
    out = {}
    for name, c in COMPONENTS.items():
        if c.get("zero_row"):
            r = t4[t4["blend_name"] == c["zero_row"]].iloc[0]
            ron, mon, src = float(r["ron"]), float(r["mon"]), f"Table 4 row {c['zero_row']}"
        elif c.get("prf_definition"):
            ron, mon = c["prf_definition"]
            src = "PRF scale definition"
        else:
            r = t3[t3["s_no"] == c["t3"]].iloc[0]
            ron, mon, src = float(r["ron"]), float(r["mon"]), f"Table 3 row {c['t3']}"
        out[name] = {"ron": ron, "mon": mon, "source": src}
    return out


def _density_sources():
    """Densities (kg/m3 near 293 K): production median where available, else AI4fuel (CRC handbook)."""
    ik = component_inchikeys()
    m = pd.read_csv(UNIFIED_MEAS, usecols=["property", "value", "T_K", "comp_inchikeys"])
    m = m[(m["property"] == "density_kg_m3") & m["T_K"].between(292.0, 294.5)]
    ai = pd.read_csv(STAGING_DIR / "ai4fuel" / "ai4fuel_density_g_cm3.csv")
    from rdkit import Chem
    ai["ik"] = [Chem.MolToInchiKey(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else "" for s in ai["smiles"]]
    out = {}
    for name, key in ik.items():
        x = m[m["comp_inchikeys"] == json.dumps([key])]
        if len(x):
            out[name] = (float(x["value"].median()), "data/unified (median, 292-294.5 K)")
            continue
        y = ai[ai["ik"] == key]
        if len(y):
            out[name] = (float(y["value"].median()) * 1000.0, "AI4fuel CRC handbook (15-30 C)")
            continue
        out[name] = (None, "none")
    return out


# Documented fallback when neither production nor AI4fuel has the compound (flagged in the report).
UNSOURCED_DENSITY_FALLBACK = {"diisobutylene": 715.0}


def densities() -> dict:
    """name -> (density kg/m3, source string)."""
    d = _density_sources()
    for k, (v, s) in list(d.items()):
        if v is None:
            d[k] = (UNSOURCED_DENSITY_FALLBACK[k], "UNSOURCED literal fallback")
    return d
