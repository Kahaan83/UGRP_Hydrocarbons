#!/usr/bin/env python
"""Extract solid-liquid equilibrium (liquidus / eutectic) temperatures of hydrocarbon MIXTURES from ThermoML archives.

    python -m fuelblend.data.thermoml_sle                       # both default archives -> data/processed/v5/thermoml_sle_hc.csv
    python -m fuelblend.data.thermoml_sle ARCHIVE.tgz [ARCHIVE2.tgz ...] [--out CSV] [--max-pressure-kpa 1000]

How ThermoML encodes SLE (checked on real records, e.g. 10.1016/j.fluid.2007.01.029, benzene + cyclohexane):
  * property  PhaseTransition / ePropName "Solid-liquid equilibrium temperature, K", ePropPhase Liquid
  * <PhaseID> blocks: one "Liquid" and one "Crystal" with RegNum/nOrgNum = the component that crystallises
    (several Crystal phases = several solids; a Crystal PhaseID without RegNum = solid not identified)
  * the composition variable (Variable > VariableID > VariableType > eComponentComposition, with VarPhaseID/eVarPhase
    = Liquid) is the LIQUID composition at the liquidus; the last component of a multi-component set is implied
  * eutectic points are separate datasets: "Eutectic temperature, K" (no composition variable, ePropPhase Crystal)
    and, usually in another dataset of the same file, "Eutectic composition: mole|mass fraction" for ONE component
    (nOrgNum under Property-MethodID); we pair them when the pairing is unique within the file and component set.
Pure-component "Normal melting temperature, K" rows (n_comp == 1) are not mixtures and are only counted.

Output columns: file, doi, n_comp, inchis, names, smiles, x_liquid, T_K, P_kPa, solid_phase_desc, point_kind
(liquidus | eutectic | other), is_hc, system_id, plus crystal_inchikeys, comp_scale (original scale), prop_name, prop_phase.
Mass fractions are converted to mole fractions with RDKit molecular weights; volume fractions are not convertible and dropped.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import io
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from fuelblend.data.thermoml_extract import (
    _assign, any_phase, first_text, is_hydrocarbon, iter_xml, kids, to_float, type_of,
)

RDLogger.DisableLog("rdApp.*")

DEFAULT_ARCHIVES = [
    "extras/fuelblend_data/data/raw/thermoml/ThermoML.v2020-09-30.tgz",
    "data/raw/thermoml_2020/thermoml_jced_2020_unlisted.tgz",
]
OUT = Path("data/processed/v5/thermoml_sle_hc.csv")

P_SLE = "Solid-liquid equilibrium temperature, K"
P_EUT_T = "Eutectic temperature, K"
P_EUT_X = {"Eutectic composition: mole fraction": "mole fraction", "Eutectic composition: mass fraction": "mass fraction"}
P_MELT = "Normal melting temperature, K"
WANTED_PROPS = (P_SLE, P_EUT_T, *P_EUT_X, P_MELT)
QUICK_BYTES = tuple(p.encode() for p in (P_SLE, P_EUT_T, "Eutectic composition"))

_MOL: dict = {}


def mol_info(inchi):
    """(inchikey, mw, smiles, is_hc) from a standard InChI, cached."""
    if inchi not in _MOL:
        m = Chem.MolFromInchi(inchi) if inchi else None
        if m is None:
            _MOL[inchi] = (None, None, None, False)
        else:
            _MOL[inchi] = (Chem.InchiToInchiKey(inchi), Descriptors.MolWt(m), Chem.MolToSmiles(m), is_hydrocarbon(m))
    return _MOL[inchi]


def system_id(iks) -> str:
    """Same rule as data/unified: md5 of the sorted InChIKeys joined by '|', first 12 hex chars."""
    return hashlib.md5("|".join(sorted(iks)).encode()).hexdigest()[:12]


def _methodname_org(p):
    """nOrgNum directly under Property-MethodID (the component an eutectic-composition value refers to)."""
    for pm in kids(p, "Property-MethodID"):
        for r in kids(pm, "RegNum"):
            n = first_text(r, "nOrgNum")
            if n:
                return n
    return None


def parse_file(root, fname: str, stats: collections.Counter, max_p: float = 1000.0, drops: list | None = None):
    """Return (compounds, doi, rows).  Each row: dict with ids, comp (liquid composition by id, on scale `scale`), T, P, ..."""
    doi = first_text(root, "sDOI")
    drops = drops if drops is not None else []
    compounds = {}
    for c in kids(root, "Compound"):
        compounds[first_text(c, "nOrgNum")] = dict(inchi=first_text(c, "sStandardInChI"), name=first_text(c, "sCommonName"))
    rows, eut_x = [], []
    for d in kids(root, "PureOrMixtureData"):
        ids = [first_text(c, "nOrgNum") for c in kids(d, "Component")]
        props = {}
        for p in kids(d, "Property"):
            pn = first_text(p, "ePropName") or ""
            if pn in WANTED_PROPS:
                props[first_text(p, "nPropNumber")] = (pn, first_text(p, "ePropPhase"), _methodname_org(p))
        if not props:
            continue
        crystal_ids = []
        for ph in kids(d, "PhaseID"):
            t = first_text(ph, "ePhase")
            if t and t.strip().lower().startswith("crystal"):
                rn = [first_text(r, "nOrgNum") for r in kids(ph, "RegNum")]
                crystal_ids.extend(rn if rn else ["?"])
        defs, consts = {}, []
        for v in kids(d, "Variable"):
            vid = kids(v, "VariableID")
            kind, text = type_of(vid[0]) if vid else (None, None)
            defs[first_text(v, "nVarNumber")] = (kind, text, first_text(v, "nOrgNum"), any_phase(v, "VarPhaseID"))
        for c in kids(d, "Constraint"):
            cid = kids(c, "ConstraintID")
            kind, text = type_of(cid[0]) if cid else (None, None)
            consts.append((kind, text, first_text(c, "nOrgNum"), to_float(first_text(c, "nConstraintValue")),
                           any_phase(c, "ConstraintPhaseID")))
        nonliquid_comp = any(k in ("eComponentComposition", "eSolventComposition") and ph is not None
                             and not ph.strip().lower().startswith("liquid") for k, _, _, ph in
                             [(v[0], v[1], v[2], v[3]) for v in defs.values()] + [(c[0], c[1], c[2], c[4]) for c in consts])
        for nv in kids(d, "NumValues"):
            state, comp = {}, {}
            for kind, text, org, val, ph in consts:
                _assign(state, comp, kind, text, org, val, ph)
            for vv in kids(nv, "VariableValue"):
                kind, text, org, ph = defs.get(first_text(vv, "nVarNumber"), (None, None, None, None))
                _assign(state, comp, kind, text, org, to_float(first_text(vv, "nVarValue")), ph)
            for pv in kids(nv, "PropertyValue"):
                key = first_text(pv, "nPropNumber")
                if key not in props:
                    continue
                pname, pphase, eorg = props[key]
                val = to_float(first_text(pv, "nPropValue"))
                if val is None:
                    stats["no numeric value"] += 1
                    continue
                n = len(ids)
                if pname == P_MELT:
                    stats["pure melting temperature rows (n_comp==1; not extracted)" if n == 1 else "melting temp. with n_comp>1 (not extracted)"] += 1
                    continue
                if pname in P_EUT_X:                      # composition of the eutectic: one component's fraction
                    if eorg is not None and n >= 2:
                        eut_x.append(dict(ids=tuple(sorted(ids)), org=eorg, x=val, scale=P_EUT_X[pname]))
                    continue
                P = state.get("P")
                if P is not None and P > max_p:
                    drops.append(("pressure > %g kPa" % max_p, ids))
                    continue
                if n == 1:
                    stats["pure component (not a mixture)"] += 1
                    continue
                r = dict(file=fname, doi=doi, prop=pname, prop_phase=pphase, value=val, T_K=state.get("T"), P_kPa=P, ids=ids,
                         crystal_ids=crystal_ids, kind="liquidus" if pname == P_SLE else "eutectic", comp=None, scale=None)
                if pname == P_SLE:
                    if state.get("bad_comp"):
                        drops.append(("unsupported composition scale (molality/ratio/...)", ids))
                        continue
                    sc = state.get("scales", ())
                    if len(sc) != 1:
                        drops.append(("no composition variable / scale" if not sc else "mixed composition scales", ids))
                        continue
                    c = dict(comp)
                    if len(c) == n - 1:
                        missing = [i for i in ids if i not in c]
                        rest = 1.0 - sum(c.values())
                        if len(missing) != 1 or rest < -1e-6:
                            drops.append(("inconsistent composition (sum > 1)", ids))
                            continue
                        c[missing[0]] = max(rest, 0.0)
                    elif len(c) != n:
                        drops.append(("composition not given for the Liquid phase (non-liquid composition ignored)" if nonliquid_comp
                                      else "composition incomplete (neither n nor n-1 fractions)", ids))
                        continue
                    if set(c) != set(ids):
                        drops.append(("composition ids do not match components", ids))
                        continue
                    r["comp"], r["scale"] = c, next(iter(sc))
                rows.append(r)
    # pair eutectic temperatures with a uniquely identifiable eutectic composition of the same component set
    for r in rows:
        if r["kind"] != "eutectic":
            continue
        cands = [e for e in eut_x if e["ids"] == tuple(sorted(r["ids"]))]
        if len(r["ids"]) == 2 and len(cands) == 1:
            e = cands[0]
            other = [i for i in r["ids"] if i != e["org"]]
            if len(other) == 1 and 0.0 <= e["x"] <= 1.0:
                r["comp"], r["scale"] = {e["org"]: e["x"], other[0]: 1.0 - e["x"]}, e["scale"]
    return compounds, doi, rows


def to_output(rows, compounds, stats):
    """Apply component identity / hydrocarbon / scale conversion; return list of output dicts."""
    out = []
    for r in rows:
        ids = r["ids"]
        info = [mol_info(compounds.get(i, {}).get("inchi")) for i in ids]
        if any(x[0] is None for x in info):
            stats["component without usable InChI (dropped)"] += 1
            continue
        is_hc = all(x[3] for x in info)
        if not is_hc:
            stats["not all-hydrocarbon (dropped from main CSV)"] += 1
            continue
        if r["kind"] == "eutectic" and r["comp"] is None:
            stats["eutectic temperature without a unique eutectic composition (kept, x_liquid empty)"] += 1
        x = None
        if r["comp"] is not None:
            if r["scale"] == "mole fraction":
                x = [r["comp"][i] for i in ids]
            elif r["scale"] == "mass fraction":
                w = [r["comp"][i] / m[1] for i, m in zip(ids, info)]
                s = sum(w)
                x = [v / s for v in w] if s > 0 else None
            else:
                stats["volume fraction composition (not convertible, dropped)"] += 1
                continue
            if x is not None:
                s = sum(x)
                if abs(s - 1.0) > 5e-3:
                    stats["mole fractions do not sum to 1 (dropped)"] += 1
                    continue
        cr = [i for i in r["crystal_ids"]]
        parts = [(compounds.get(i, {}).get("name") or "unidentified") if i != "?" else "unidentified" for i in cr]
        desc = ("Crystal: " + "; ".join(parts)) if cr else "no crystal phase stated"
        if r["kind"] == "eutectic":
            desc = "eutectic (" + desc + ")"
        iks = [m[0] for m in info]
        out.append(dict(
            file=r["file"], doi=r["doi"], n_comp=len(ids),
            inchis=json.dumps([compounds[i]["inchi"] for i in ids]),
            names=json.dumps([compounds[i]["name"] for i in ids]),
            smiles=json.dumps([m[2] for m in info]),
            x_liquid=json.dumps(x) if x is not None else "",
            T_K=r["value"],
            P_kPa=r["P_kPa"], solid_phase_desc=desc, point_kind=r["kind"], is_hc=True, system_id=system_id(iks),
            inchikeys=json.dumps(iks),
            crystal_inchikeys=json.dumps([mol_info(compounds.get(i, {}).get("inchi"))[0] if i != "?" else None for i in cr]),
            comp_scale=r["scale"], prop_name=r["prop"], prop_phase=r["prop_phase"]))
    return out


def extract(archives, max_p: float = 1000.0):
    stats: collections.Counter = collections.Counter()
    out, n_files, n_bad, n_sle_files = [], 0, 0, 0
    for arch in archives:
        for fname, blob in iter_xml(Path(arch)):
            n_files += 1
            if not any(q in blob for q in QUICK_BYTES) and b"Normal melting temperature" not in blob:
                continue
            try:
                root = ET.parse(io.BytesIO(blob)).getroot()
                drops: list = []
                comps, doi, rows = parse_file(root, fname, stats, max_p, drops)
            except Exception:  # noqa: BLE001 - one odd file must not kill the run
                n_bad += 1
                continue
            for reason, ids in drops:                      # parse-stage drops: count hydrocarbon systems only
                ok = bool(ids) and all(mol_info(comps.get(i, {}).get("inchi"))[3] for i in ids)
                stats[reason + ("" if ok else " [non-HC]")] += 1
            if rows:
                n_sle_files += 1
            out.extend(to_output(rows, comps, stats))
    stats["files read"], stats["files unparseable"], stats["files with SLE rows"] = n_files, n_bad, n_sle_files
    df = pd.DataFrame(out)
    if len(df):
        key = ["file", "point_kind", "inchikeys", "x_liquid", "T_K", "P_kPa", "solid_phase_desc"]
        n0 = len(df)
        df = df.drop_duplicates(key).reset_index(drop=True)
        stats["exact duplicate rows removed (same file/point/x/T/P)"] = n0 - len(df)
    return df, stats


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("archives", nargs="*")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--max-pressure-kpa", type=float, default=1000.0)
    a = ap.parse_args(argv)
    archives = a.archives or [p for p in DEFAULT_ARCHIVES if Path(p).exists()]
    if not archives:
        sys.exit("no ThermoML archive found")
    df, stats = extract(archives, a.max_pressure_kpa)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"wrote {a.out}: {len(df)} rows")
    for k, v in stats.most_common():
        print(f"  {v:>8}  {k}")


if __name__ == "__main__":
    main()
