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

pure_props_v4.csv rule (R16, one row per canonical smiles and property): `value` is ALWAYS the listed value;
the reference never overwrites it. status: verified (within threshold) | disputed (flagged, or two listings
disagree; listed_value_alt records the second) | unverified (no reference, or Psat reference < 0.01 kPa where the
relative metric is noise) | not_comparable (mixture/commercial/surrogate) | lookup_error (library raised; error
class in ref_source) | reference_only (pure compound, no listed value, value filled from the reference).
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


# Exceptions the libraries raise for "cannot compute this property for this compound" at the
# call level; anything outside this tuple is a bug and propagates.
LOOKUP_ERRORS = (ValueError, KeyError, TypeError, AttributeError, IndexError, ZeroDivisionError,
                 OverflowError, NotImplementedError, RuntimeError)
# Psat: below this reference value (kPa) a relative comparison is noise (R16)
PSAT_MIN_REF_KPA = 0.01
# R16: isomer mixtures / commercial grades / stereo-unspecified compounds would be compared against a pure
# reference meaninglessly. Names containing these words, names in MIXTURE_NAMES, or a listed MW > 3 % off
# the SMILES MW make a row kind="mixture_or_surrogate".
MIXTURE_WORDS = ("commercial", "mixed", "surrogate")
# Only stereo-unspecified decalin (merges cis and trans, CAS 91-17-8 is the isomer mix) is listed by bare name.
# Bare isododecane and heptamethylnonane are single-isomer SMILES whose listed MW matches, so they are compared.
MIXTURE_NAMES = {"decalin"}


def _lookup(smiles: str, prop: str):
    """Return (value, source, state); state in ok | no_reference | lookup_error:<ExceptionClass>."""
    if prop not in THRESH:
        return None, "", "no_reference"
    cs = canon(smiles)
    if cs is None:
        return None, "", "no_reference"
    import chemicals
    import thermo
    chem = _resolve(cs)
    if chem is None:
        return None, "", "no_reference"
    cas = chem.CASs
    cv = chemicals.__version__
    try:
        if prop == "Tm":
            v = _finite(chemicals.Tm(cas))
            r = (v - 273.15, SRC_CHEM.format(v=cv, fn="Tm")) if v is not None else None
        elif prop == "Tb":
            v = _finite(chemicals.Tb(cas))
            r = (v - 273.15, SRC_CHEM.format(v=cv, fn="Tb")) if v is not None else None
        elif prop == "MW":
            v = _finite(chem.MW)
            r = (v, SRC_CHEM.format(v=cv, fn="search_chemical.MW")) if v is not None else None
        elif prop == "dHfus":
            v = _finite(chemicals.Hfus(cas))
            r = (v / 1000.0, SRC_CHEM.format(v=cv, fn="Hfus")) if v is not None else None
        elif prop == "flash":
            v = _finite(chemicals.T_flash(cas))
            r = (v - 273.15, SRC_CHEM.format(v=cv, fn="T_flash")) if v is not None else None
            if r is None:
                df = _cooptima()
                if df is not None:
                    s = df[(df.cas == cas) & df.flash_point_C.notna()]
                    if len(s):
                        r = (float(s.flash_point_C.iloc[0]), SRC_COOPTIMA)
        elif prop == "rho15":
            v = _finite(thermo.Chemical(cas, T=T_RHO).rhol)
            r = (v, SRC_THERMO.format(v=thermo.__version__, fn="rhol @288.15K")) if v and v > 0 else None
        elif prop == "Psat37_8":
            v = _finite(thermo.Chemical(cas, T=T_PSAT).Psat)
            r = (v / 1000.0, SRC_THERMO.format(v=thermo.__version__, fn="Psat @310.95K")) if v and v > 0 else None
        elif prop in ("RON", "MON"):
            df = _cooptima()
            s = df[(df.cas == cas) & df[prop.lower()].notna()] if df is not None else None
            r = (float(s[prop.lower()].iloc[0]), SRC_COOPTIMA) if s is not None and len(s) else None
        elif prop == "DCN":
            df = _nrel()
            s = df[(df._c == cs) & df.dcn.notna()] if df is not None else None
            r = (float(s.dcn.median()), SRC_NREL) if s is not None and len(s) else None
        else:
            raise ValueError(f"no lookup implemented for {prop}")
    except LOOKUP_ERRORS as e:
        return None, f"lookup_error:{type(e).__name__}", f"lookup_error:{type(e).__name__}"
    return (r[0], r[1], "ok") if r else (None, "", "no_reference")


def reference_value(smiles: str, prop: str):
    """Reference (value, source) for one pure compound and property, or None if unavailable
    (no data or a lookup error; use _lookup() to tell the two apart)."""
    v, src, state = _lookup(smiles, prop)
    return (v, src) if state == "ok" else None


def is_flagged(prop: str, abs_delta, rel_delta, ref_value=None) -> bool:
    kind, thr = THRESH[prop]
    if prop == "Psat37_8" and ref_value is not None and ref_value < PSAT_MIN_REF_KPA:
        return False
    d = abs_delta if kind == "abs" else rel_delta
    return d is not None and abs(d) > thr


def is_mixture(name: str, smiles: str, listed_mw) -> bool:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    n = str(name).lower()
    if any(w in n for w in MIXTURE_WORDS) or n in MIXTURE_NAMES:
        return True
    m = Chem.MolFromSmiles(smiles)
    lm = _finite(listed_mw)
    return bool(m is not None and lm and abs(lm - Descriptors.MolWt(m)) / Descriptors.MolWt(m) > 0.03)


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
                                 listed_source="purchasable_compounds.csv", listed_mw=r.mw))
    from fuelblend.design.inverse_design_v2 import CANDIDATE_REGISTRY as R
    rmap = {"mw": "MW", "rho_15": "rho15", "tm_C": "Tm", "tb_C": "Tb", "hfus_kJ_mol": "dHfus",
            "psat_378_kpa": "Psat37_8", "dcn": "DCN", "ron": "RON", "mon": "MON"}
    for name, d in R.items():
        for key, prop in rmap.items():
            v = _finite(d.get(key))
            if v is not None:
                rows.append(dict(compound=name, smiles=d["smiles"], prop=prop, listed_value=v,
                                 listed_source="inverse_design_v2.CANDIDATE_REGISTRY", listed_mw=d.get("mw")))
    return rows


def audit() -> pd.DataFrame:
    out = []
    for r in _listed_rows():
        r = dict(r)
        kind = THRESH[r["prop"]][0]
        r["kind"] = "mixture_or_surrogate" if is_mixture(r["compound"], r["smiles"], r.pop("listed_mw")) else "pure"
        if r["kind"] != "pure":
            r.update(ref_value=np.nan, ref_source="", delta=np.nan, flag="not_comparable")
        else:
            val, src, state = _lookup(r["smiles"], r["prop"])
            if state != "ok":
                r.update(ref_value=np.nan, ref_source=src,
                         delta=np.nan, flag="no_reference" if state == "no_reference" else "lookup_error")
            else:
                absd = r["listed_value"] - val
                reld = absd / val if val else None
                if r["prop"] == "Psat37_8" and val < PSAT_MIN_REF_KPA:
                    flag = "not_assessed"
                else:
                    flag = "flagged" if is_flagged(r["prop"], absd, reld, val) else "ok"
                r.update(ref_value=val, ref_source=src, delta=absd if kind == "abs" else reld, flag=flag)
        r["delta_kind"] = "absolute" if kind == "abs" else "relative"
        out.append(r)
    cols = ["compound", "smiles", "prop", "listed_value", "listed_source", "ref_value", "ref_source",
            "delta", "flag", "delta_kind", "kind"]
    return pd.DataFrame(out)[cols]


_STATUS = {"ok": "verified", "flagged": "disputed", "no_reference": "unverified", "not_assessed": "unverified",
           "lookup_error": "lookup_error", "not_comparable": "not_comparable"}
REGISTRY_TAG = "inverse_design_v2.CANDIDATE_REGISTRY"


def build_props(diffs: pd.DataFrame) -> pd.DataFrame:
    """One row per (canonical smiles, prop). value is ALWAYS the listed value (R16, never auto-substituted),
    except status=reference_only (no listed value anywhere, pure compound, reference exists)."""
    d = diffs.copy()
    d["_cs"] = d.smiles.map(canon)
    out = []
    for (cs, prop), g in d.groupby(["_cs", "prop"], sort=False):
        g = g.sort_values("listed_source", key=lambda s: s.ne(REGISTRY_TAG), kind="stable")
        a = g.iloc[0]
        row = dict(compound=a.compound, smiles=cs, prop=prop, value=a.listed_value, listed_value=a.listed_value,
                   listed_source=a.listed_source, listed_value_alt=np.nan, listed_source_alt="",
                   ref_value=a.ref_value, ref_source=a.ref_source, delta=a.delta, kind=a.kind,
                   status=_STATUS[a.flag], source=a.listed_source, note="")
        if len(g) > 1 and a.kind == "pure":
            b = g.iloc[1]
            diff = b.listed_value - a.listed_value
            kd, thr = THRESH[prop]
            dis = abs(diff) > thr if kd == "abs" else abs(diff) > thr * abs(a.listed_value or 1)
            if dis or b.flag == "flagged":
                row.update(listed_value_alt=b.listed_value, listed_source_alt=b.listed_source,
                           status="disputed",
                           note="two listed values disagree" if dis else "second listing flagged against reference")
        if row["status"] == "disputed" and not row["note"]:
            row["note"] = "listed value differs from reference beyond threshold; value kept as listed"
        out.append(row)
    props = pd.DataFrame(out)
    # reference_only: pure compounds with a reference but no listed value for that property
    extra = []
    pure_cs = set(d._cs) - set(d.loc[d.kind != "pure", "_cs"])
    for cs, g in d[d._cs.isin(pure_cs)].groupby("_cs"):
        have = set(g.prop)
        for prop in THRESH:
            if prop in have:
                continue
            val, src, state = _lookup(cs, prop)
            if state == "ok" and not (prop == "Psat37_8" and val < PSAT_MIN_REF_KPA):
                extra.append(dict(compound=g.compound.iloc[0], smiles=cs, prop=prop, value=val, listed_value=np.nan,
                                  listed_source="", listed_value_alt=np.nan, listed_source_alt="", ref_value=val,
                                  ref_source=src, delta=np.nan, kind="pure", status="reference_only",
                                  source=src, note="no listed value; filled from reference"))
    if extra:
        props = pd.concat([props, pd.DataFrame(extra)], ignore_index=True)
    return props


def run(write: bool = True, out_diffs: Path = OUT_DIFFS, out_props: Path = OUT_PROPS) -> dict:
    diffs = audit()
    props = build_props(diffs)
    if write:
        out_diffs.parent.mkdir(parents=True, exist_ok=True)
        out_props.parent.mkdir(parents=True, exist_ok=True)
        diffs.to_csv(out_diffs, index=False)
        props.to_csv(out_props, index=False)
    return {"diffs": diffs, "props": props}


def main() -> int:
    res = run(write=True)
    d, p = res["diffs"], res["props"]
    print(f"audited rows: {len(d)}  compounds: {d.smiles.map(canon).nunique()}")
    print(d.groupby(["prop", "flag"]).size().unstack(fill_value=0).to_string())
    print("lookup_error rows:", int((d.flag == "lookup_error").sum()))
    print(p.status.value_counts().to_string())
    print(f"files: {OUT_DIFFS}\n       {OUT_PROPS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
