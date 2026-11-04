#!/usr/bin/env python
"""Extract hydrocarbon-only liquid viscosity / density / vapour-pressure records from the
NIST ThermoML archive (.tgz of ThermoML XML files).

    python thermoml_extract.py                                  # data/raw/thermoml/*.tgz
    python thermoml_extract.py path/to/archive.tgz --props density viscosity vapor
    python thermoml_extract.py --max-pressure-kpa 500           # default 1000 (atmospheric-ish liquids)
    python thermoml_extract.py --debug                          # also print the first parsed rows

Output: data/processed/thermoml_hc.csv  (one row per data point)
  file, prop, prop_phase, value, T_K, P_kPa, n_comp, inchis, names, comp_type, comp_values

Structure facts this parser relies on (checked against the official ThermoML.xsd and a real record):
  * fixed conditions are  Constraint > ConstraintID > ConstraintType > eTemperature|ePressure|eComponentComposition
    (NOT a direct ConstraintType child - the old version of this script missed every fixed T/P)
  * varied conditions are Variable > VariableID > VariableType > e...   with the component in RegNum/nOrgNum
  * composition variables carry a phase (VarPhaseID/ePhase): vapour-phase mole fractions of VLE data are IGNORED
  * each property carries a phase (PropPhaseID/ePropPhase): only 'Liquid' densities/viscosities are kept and
    crystal-phase (sublimation) vapour pressures are dropped
Only these exact ThermoML property names are extracted (anything else, e.g. 'Excess viscosity', is ignored):
  Mass density, kg/m3 | Viscosity, Pa*s | Kinematic viscosity, m2/s | Vapor or sublimation pressure, kPa
"""
from __future__ import annotations
import argparse, collections, io, json, sys, tarfile, traceback
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from common import is_hydrocarbon, mol_from

OUT = Path("data/processed/thermoml_hc.csv")

WANTED = {                       # exact ePropName -> keyword family used by --props
    "Mass density, kg/m3": "density",
    "Viscosity, Pa*s": "viscosity",
    "Kinematic viscosity, m2/s": "viscosity",
    "Vapor or sublimation pressure, kPa": "vapor",
}
SCALES = {"mole fraction", "mass fraction", "volume fraction"}


def ln(el) -> str:
    return el.tag.split("}")[-1]


def kids(el, name):
    return [c for c in el if ln(c) == name]


def first_text(el, name):
    for e in el.iter():
        if ln(e) == name and e.text and e.text.strip():
            return e.text.strip()
    return None


def type_of(container):
    """(tag, text) of the first enum-like descendant (local name starts with 'e'), e.g. ('eTemperature', 'Temperature, K')."""
    if container is None:
        return None, None
    for e in container.iter():
        if ln(e).startswith("e") and e.text and e.text.strip():
            return ln(e), e.text.strip()
    return None, None


def phase_of(parent, wrapper):
    """Text of ePhase inside e.g. VarPhaseID / ConstraintPhaseID, or None."""
    for w in kids(parent, wrapper):
        t = first_text(w, "ePhase")
        if t:
            return t
    return None


_HC_CACHE: dict = {}


def hc_inchi(inchi) -> bool:
    if inchi not in _HC_CACHE:
        _HC_CACHE[inchi] = bool(inchi) and is_hydrocarbon(mol_from(inchi=inchi))
    return _HC_CACHE[inchi]


def to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def is_liquid(phase) -> bool:
    return phase is None or phase.strip().lower().startswith("liquid")


def _assign(state, comp, kind, text, org, val, phase):
    """Record one fixed/varied condition into state (T, P) or comp (liquid composition)."""
    if val is None or not kind:
        return
    if kind == "eTemperature":
        state["T"] = val
    elif kind == "ePressure":
        state["P"] = val
    elif kind in ("eComponentComposition", "eSolventComposition"):
        scale = (text or "").lower().split(":")[-1].strip()
        if not is_liquid(phase):
            return                                   # vapour-phase composition of VLE data: ignore
        if scale not in SCALES or org is None:
            state["bad_comp"] = True                 # molality, molarity, ratios ... not supported
            return
        comp[org] = val
        state.setdefault("scales", set()).add(scale)


def parse_file(root, fname: str, wanted: dict[str, str], max_p: float, stats: collections.Counter):
    compounds = {}
    for c in kids(root, "Compound"):
        compounds[first_text(c, "nOrgNum")] = dict(inchi=first_text(c, "sStandardInChI"), name=first_text(c, "sCommonName"))

    rows = []
    for d in kids(root, "PureOrMixtureData"):
        ids = [first_text(c, "nOrgNum") for c in kids(d, "Component")]
        defs, consts = {}, []
        for v in kids(d, "Variable"):
            vid = kids(v, "VariableID")
            kind, text = type_of(vid[0]) if vid else (None, None)
            defs[first_text(v, "nVarNumber")] = (kind, text, first_text(v, "nOrgNum"), phase_of(v, "VarPhaseID"))
        for c in kids(d, "Constraint"):
            cid = kids(c, "ConstraintID")                     # <- the real container (old code looked for ConstraintType directly)
            kind, text = type_of(cid[0]) if cid else (None, None)
            consts.append((kind, text, first_text(c, "nOrgNum"), to_float(first_text(c, "nConstraintValue")),
                           phase_of(c, "ConstraintPhaseID")))
        props = {}
        for p in kids(d, "Property"):
            props[first_text(p, "nPropNumber")] = (first_text(p, "ePropName") or "", first_text(p, "ePropPhase"))

        for nv in kids(d, "NumValues"):
            state, comp = {}, {}
            for kind, text, org, val, ph in consts:
                _assign(state, comp, kind, text, org, val, ph)
            for vv in kids(nv, "VariableValue"):
                kind, text, org, ph = defs.get(first_text(vv, "nVarNumber"), (None, None, None, None))
                _assign(state, comp, kind, text, org, to_float(first_text(vv, "nVarValue")), ph)
            for pv in kids(nv, "PropertyValue"):
                pname, pphase = props.get(first_text(pv, "nPropNumber"), ("", None))
                fam = wanted.get(pname)
                if fam is None:
                    stats["other property (ignored)"] += 1
                    continue
                val = to_float(first_text(pv, "nPropValue"))
                if val is None:
                    stats["no numeric value"] += 1
                    continue
                if fam == "vapor":
                    if pphase and "crystal" in pphase.lower():
                        stats["vapour pressure of a solid (dropped)"] += 1
                        continue
                else:
                    if not is_liquid(pphase):
                        stats[f"{fam}: non-liquid phase (dropped)"] += 1
                        continue
                    if state.get("P") is not None and state["P"] > max_p:
                        stats[f"{fam}: pressure > {max_p:g} kPa (dropped)"] += 1
                        continue
                n = len(ids)
                c = dict(comp)
                if n == 1:
                    c, scale = {ids[0]: 1.0}, "pure"
                else:
                    if state.get("bad_comp") or len(state.get("scales", ())) != 1:
                        stats["mixture: unsupported/ambiguous composition scale"] += 1
                        continue
                    scale = next(iter(state["scales"]))
                    if len(c) == n - 1:                          # last component implied
                        missing = [i for i in ids if i not in c]
                        rest = 1.0 - sum(c.values())
                        if len(missing) != 1 or rest < -1e-6:
                            stats["mixture: inconsistent composition"] += 1
                            continue
                        c[missing[0]] = max(rest, 0.0)
                    elif len(c) != n:
                        stats["mixture: composition incomplete"] += 1
                        continue
                rows.append(dict(file=fname, prop=pname, prop_phase=pphase, value=val, T_K=state.get("T"),
                                 P_kPa=state.get("P"), ids=ids, comp=c, comp_type=scale))
    return compounds, rows


def iter_xml(path: Path):
    if path.suffix in (".tgz", ".gz") or path.name.endswith(".tar.gz"):
        with tarfile.open(path) as tar:
            for m in tar:
                if m.isfile() and m.name.lower().endswith((".xml", ".thermoml")):
                    yield m.name, tar.extractfile(m).read()
    else:
        yield path.name, path.read_bytes()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", nargs="?")
    ap.add_argument("--props", nargs="+", default=["viscosity", "density", "vapor"], choices=["viscosity", "density", "vapor"])
    ap.add_argument("--max-pressure-kpa", type=float, default=1000.0)
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args(argv)
    wanted = {k: v for k, v in WANTED.items() if v in a.props}
    path = Path(a.archive) if a.archive else next(Path("data/raw/thermoml").glob("*.tgz"), None)
    if path is None or not path.exists():
        sys.exit("No ThermoML archive found. Run fetch.py thermoml or place the .tgz in data/raw/thermoml/")

    out, n_files, n_bad = [], 0, 0
    stats: collections.Counter = collections.Counter()
    errors: list[str] = []
    for fname, blob in iter_xml(path):
        n_files += 1
        if n_files % 2000 == 0:
            print(f"  ...{n_files} files, {len(out)} hydrocarbon rows so far", flush=True)
        try:
            root = ET.parse(io.BytesIO(blob)).getroot()
            comps, rows = parse_file(root, fname, wanted, a.max_pressure_kpa, stats)
        except ET.ParseError:
            n_bad += 1
            continue
        except Exception as e:  # noqa: BLE001 - one odd file must not kill a 20k-file run
            n_bad += 1
            if len(errors) < 5:
                errors.append(f"{fname}: {type(e).__name__}: {e}\n{traceback.format_exc(limit=2)}")
            continue
        for r in rows:
            ok = bool(r["ids"]) and all(hc_inchi(comps.get(i, {}).get("inchi")) for i in r["ids"])
            if not ok:
                stats["not all-hydrocarbon / no InChI"] += 1
                continue
            out.append(dict(
                file=r["file"], prop=r["prop"], prop_phase=r["prop_phase"], value=r["value"], T_K=r["T_K"], P_kPa=r["P_kPa"],
                n_comp=len(r["ids"]), inchis=json.dumps([comps[i]["inchi"] for i in r["ids"]]),
                names=json.dumps([comps[i]["name"] for i in r["ids"]]), comp_type=r["comp_type"],
                comp_values=json.dumps([r["comp"].get(i) for i in r["ids"]])))
    df = pd.DataFrame(out)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"files read: {n_files}  unparseable/errors: {n_bad}  hydrocarbon rows kept: {len(df)}")
    if len(df):
        print(df.groupby(["prop", "n_comp"]).size().to_string())
        have_T = df.T_K.notna().mean()
        print(f"rows with temperature: {have_T:.1%}   with pressure: {df.P_kPa.notna().mean():.1%}")
    print("\nrows NOT kept, by reason:")
    for k, v in stats.most_common():
        print(f"  {v:>9}  {k}")
    if errors:
        print("\nfirst file errors:\n" + "\n".join(errors))
    if a.debug and len(df):
        print("\nfirst rows:\n" + df.head(8).to_string())


if __name__ == "__main__":
    main()
