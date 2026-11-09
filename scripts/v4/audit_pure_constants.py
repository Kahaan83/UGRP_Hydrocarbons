#!/usr/bin/env python
"""v4 Task 4 (F1): audit listed pure-component constants against reference libraries.

Listed sources audited
  * data/processed/v3/purchasable_compounds.csv   (tag "purchasable_compounds.csv")
  * CANDIDATE_REGISTRY in src/fuelblend/design/inverse_design_v2.py (tag "inverse_design_v2.CANDIDATE_REGISTRY")

Reference sources
  * chemicals / thermo (Tm, Tb, MW, Hfus, flash via chemicals.T_flash, rho at 15 C and Psat at 37.8 C via
    thermo.Chemical), resolved from the canonical SMILES with chemicals.identifiers.search_chemical
  * RON/MON: DOE Co-Optima database file data/processed/v2/octane_co_optima.csv (joined by CAS number)
  * DCN: NREL cetane compendium data/processed/v2/cetane_nrel_compendium.csv (joined by canonical SMILES)

Units of reference_value(): Tm/Tb/flash in deg C, rho15 in kg/m3, MW in g/mol, dHfus in kJ/mol,
Psat37_8 in kPa (vapour pressure at 37.8 C = 310.95 K), RON/MON/DCN dimensionless.

Flag thresholds (|listed - reference|): Tm, Tb, flash > 2 C; rho15 > 1 %; RON, MON, DCN > 3 ON;
Psat37_8 > 5 %; MW > 1 %; dHfus > 10 %.  A row with no reference is flag="no_reference" (never silently ok).

pure_props_v4.csv rule (one row per canonical smiles and property): the reference value if one exists and the
listed value was flagged; the listed value if it agrees with the reference (source notes the verification);
the listed value with source="unverified" when there is no reference at all.
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

PURCH = ROOT / "data/processed/v3/purchasable_compounds.csv"
COOPTIMA = ROOT / "data/processed/v2/octane_co_optima.csv"
NREL = ROOT / "data/processed/v2/cetane_nrel_compendium.csv"
OUT_DIFFS = ROOT / "results/v4/P1/pure_constant_diffs.csv"
OUT_PROPS = ROOT / "data/processed/v4/pure_props_v4.csv"

T_PSAT = 37.8 + 273.15
T_RHO = 15.0 + 273.15
# prop -> (kind, threshold): "abs" threshold in the property unit, "rel" is a fraction of the reference
THRESH = {"Tm": ("abs", 2.0), "Tb": ("abs", 2.0), "flash": ("abs", 2.0), "rho15": ("rel", 0.01),
          "RON": ("abs", 3.0), "MON": ("abs", 3.0), "DCN": ("abs", 3.0), "Psat37_8": ("rel", 0.05),
          "MW": ("rel", 0.01), "dHfus": ("rel", 0.10)}

SRC_CHEM = "chemicals {v} ({fn})"
SRC_THERMO = "thermo {v} Chemical.{fn}"
SRC_COOPTIMA = "DOE Co-Optima database (2018), data/processed/v2/octane_co_optima.csv"
SRC_NREL = "NREL cetane compendium (Yanowitz 2017), data/processed/v2/cetane_nrel_compendium.csv"


def canon(smiles) -> str | None:
    from rdkit import Chem
    if not isinstance(smiles, str) or not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m) if m is not None else None


@lru_cache(maxsize=None)
def _resolve(csmiles: str):
    # chemicals' own SMILES parser rejects aromatic (lower-case) SMILES such as benzene; resolve via InChI
    from chemicals.identifiers import search_chemical
    from rdkit import Chem
    mol = Chem.MolFromSmiles(csmiles)
    for query in ((Chem.MolToInchi(mol) if mol is not None else None), csmiles):
        if not query:
            continue
        try:
            return search_chemical(query)
        except Exception:
            continue
    return None


@lru_cache(maxsize=None)
def _cooptima():
    return pd.read_csv(COOPTIMA) if COOPTIMA.exists() else None


@lru_cache(maxsize=None)
def _nrel():
    if not NREL.exists():
        return None
    d = pd.read_csv(NREL)
    d["_c"] = d.smiles.map(canon)
    return d


def _finite(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def reference_value(smiles: str, prop: str):
    """Reference (value, source) for one pure compound and property, or None if unavailable."""
    if prop not in THRESH:
        return None
    cs = canon(smiles)
    if cs is None:
        return None
    import chemicals
    import thermo
    chem = _resolve(cs)
    if chem is None:
        return None
    cas = chem.CASs
    cv = chemicals.__version__
    try:
        if prop == "Tm":
            v = _finite(chemicals.Tm(cas))
            return (v - 273.15, SRC_CHEM.format(v=cv, fn="Tm")) if v is not None else None
        if prop == "Tb":
            v = _finite(chemicals.Tb(cas))
            return (v - 273.15, SRC_CHEM.format(v=cv, fn="Tb")) if v is not None else None
        if prop == "MW":
            v = _finite(chem.MW)
            return (v, SRC_CHEM.format(v=cv, fn="search_chemical.MW")) if v is not None else None
        if prop == "dHfus":
            v = _finite(chemicals.Hfus(cas))
            return (v / 1000.0, SRC_CHEM.format(v=cv, fn="Hfus")) if v is not None else None
        if prop == "flash":
            v = _finite(chemicals.T_flash(cas))
            if v is not None:
                return (v - 273.15, SRC_CHEM.format(v=cv, fn="T_flash"))
            df = _cooptima()
            if df is not None:
                s = df[(df.cas == cas) & df.flash_point_C.notna()]
                if len(s):
                    return (float(s.flash_point_C.iloc[0]), SRC_COOPTIMA)
            return None
        if prop == "rho15":
            v = _finite(thermo.Chemical(cas, T=T_RHO).rhol)
            return (v, SRC_THERMO.format(v=thermo.__version__, fn="rhol @288.15K")) if v and v > 0 else None
        if prop == "Psat37_8":
            v = _finite(thermo.Chemical(cas, T=T_PSAT).Psat)
            return (v / 1000.0, SRC_THERMO.format(v=thermo.__version__, fn="Psat @310.95K")) if v and v > 0 else None
        if prop in ("RON", "MON"):
            df = _cooptima()
            if df is None:
                return None
            s = df[(df.cas == cas) & df[prop.lower()].notna()]
            return (float(s[prop.lower()].iloc[0]), SRC_COOPTIMA) if len(s) else None
        if prop == "DCN":
            df = _nrel()
            if df is None:
                return None
            s = df[(df._c == cs) & df.dcn.notna()]
            return (float(s.dcn.median()), SRC_NREL) if len(s) else None
    except Exception:
        return None
    return None


def is_flagged(prop: str, abs_delta, rel_delta) -> bool:
    kind, thr = THRESH[prop]
    d = abs_delta if kind == "abs" else rel_delta
    return d is not None and abs(d) > thr


def _listed_rows() -> list:
    rows = []
    p = pd.read_csv(PURCH)
    pmap = {"mw": "MW", "density_15C_kg_m3": "rho15", "tb_C": "Tb", "tm_C": "Tm", "flash_point_C": "flash",
            "dcn": "DCN", "ron": "RON", "mon": "MON"}
    for _, r in p.iterrows():
        for col, prop in pmap.items():
            v = _finite(r[col])
            if v is not None:
                rows.append(dict(compound=r.compound_name, smiles=r.smiles, prop=prop, listed_value=v,
                                 listed_source="purchasable_compounds.csv"))
    from fuelblend.design.inverse_design_v2 import CANDIDATE_REGISTRY as R
    rmap = {"mw": "MW", "rho_15": "rho15", "tm_C": "Tm", "tb_C": "Tb", "hfus_kJ_mol": "dHfus",
            "psat_378_kpa": "Psat37_8", "dcn": "DCN", "ron": "RON", "mon": "MON"}
    for name, d in R.items():
        for key, prop in rmap.items():
            v = _finite(d.get(key))
            if v is not None:
                rows.append(dict(compound=name, smiles=d["smiles"], prop=prop, listed_value=v,
                                 listed_source="inverse_design_v2.CANDIDATE_REGISTRY"))
    return rows


def audit() -> pd.DataFrame:
    out = []
    for r in _listed_rows():
        ref = reference_value(r["smiles"], r["prop"])
        r = dict(r)
        kind = THRESH[r["prop"]][0]
        if ref is None:
            r.update(ref_value=np.nan, ref_source="", delta=np.nan, flag="no_reference")
        else:
            absd = r["listed_value"] - ref[0]
            reld = absd / ref[0] if ref[0] else None
            delta = absd if kind == "abs" else reld
            r.update(ref_value=ref[0], ref_source=ref[1], delta=delta,
                     flag="flagged" if is_flagged(r["prop"], absd, reld) else "ok")
        r["delta_kind"] = "absolute" if kind == "abs" else "relative"
        out.append(r)
    cols = ["compound", "smiles", "prop", "listed_value", "listed_source", "ref_value", "ref_source",
            "delta", "flag", "delta_kind"]
    return pd.DataFrame(out)[cols]


def build_props(diffs: pd.DataFrame) -> pd.DataFrame:
    d = diffs.copy()
    d["_cs"] = d.smiles.map(canon)
    out = []
    for (cs, prop), g in d.groupby(["_cs", "prop"], sort=False):
        # registry row first when both tables list the compound
        g = g.sort_values("listed_source", key=lambda s: s.ne("inverse_design_v2.CANDIDATE_REGISTRY"),
                          kind="stable")
        fl = g[g.flag == "flagged"]
        first = fl.iloc[0] if len(fl) else g.iloc[0]
        if first.flag == "no_reference":
            val, src = first.listed_value, "unverified"
            note = f"listed in {first.listed_source}, no reference found"
        elif first.flag == "flagged":
            val, src = first.ref_value, first.ref_source
            note = f"reference replaces flagged listed value {first.listed_value:g} ({first.listed_source})"
        else:
            val, src = first.listed_value, f"{first.listed_source} (verified vs {first.ref_source})"
            note = "listed value within threshold of reference"
        out.append(dict(compound=first.compound, smiles=cs, prop=prop, value=float(val), source=src, note=note))
    return pd.DataFrame(out)


def run(write: bool = True) -> dict:
    diffs = audit()
    props = build_props(diffs)
    if write:
        OUT_DIFFS.parent.mkdir(parents=True, exist_ok=True)
        OUT_PROPS.parent.mkdir(parents=True, exist_ok=True)
        diffs.to_csv(OUT_DIFFS, index=False)
        props.to_csv(OUT_PROPS, index=False)
    return {"diffs": diffs, "props": props}


def main() -> int:
    res = run(write=True)
    d = res["diffs"]
    print(f"audited rows: {len(d)}  compounds: {d.smiles.map(canon).nunique()}")
    print(d.groupby(["prop", "flag"]).size().unstack(fill_value=0).to_string())
    print(f"files: {OUT_DIFFS}\n       {OUT_PROPS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
