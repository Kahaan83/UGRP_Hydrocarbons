#!/usr/bin/env python
"""Pure-component melting / triple-point temperatures and enthalpies of fusion of HYDROCARBONS from ThermoML archives (v5 WS3).

    python -m fuelblend.data.thermoml_fusion [ARCHIVE.tgz ...] [--out CSV]

Properties (ThermoML ePropName, pure-component datasets only, n_comp == 1):
  * "Normal melting temperature, K"              -> prop "tm"
  * "Triple point temperature, K"                -> prop "ttp"
  * "Molar enthalpy of transition or fusion, kJ/mol" -> prop "dhfus".
  All three are kept ONLY when the dataset's phases are exactly one Crystal phase plus Liquid (ambient "Air at 1 atmosphere" / "Gas"
  ignored): ThermoML also files crystal-crystal transitions under "Triple point temperature" and "Molar enthalpy of transition"
  (phases Crystal 2 | Crystal 1), which are dropped.
Pressure > 1000 kPa dropped.  Tier name of the output: ``thermoml_measured``.  ``summarize()`` gives per-compound medians, ranges and
``spread`` flags (Tm range > 5 K, dHfus range > 15 % of the median).
"""
from __future__ import annotations

import argparse
import collections
import io
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.data.thermoml_extract import first_text, iter_xml, kids, to_float, type_of, _assign
from fuelblend.data.thermoml_sle import mol_info

DEFAULT_ARCHIVES = ["data/raw/thermoml_archive_2020/ThermoML.v2020-09-30.tgz",
                    "extras/fuelblend_data/data/raw/thermoml/ThermoML.v2020-09-30.tgz"]
OUT = Path("data/processed/v5/thermoml_fusion_pure_hc.csv")
PROPS = {"Normal melting temperature, K": "tm", "Triple point temperature, K": "ttp",
         "Molar enthalpy of transition or fusion, kJ/mol": "dhfus"}
TIER = "thermoml_measured"


def _phases(d) -> list[str]:
    return [(first_text(ph, "ePhase") or "").strip() for ph in kids(d, "PhaseID")]


def is_fusion_phases(phases) -> bool:
    """Exactly one crystal phase plus liquid.  ThermoML also lists the surrounding atmosphere ("Air at 1 atmosphere") or "Gas" as a
    phase of the dataset; those are allowed and ignored (checked on real records, e.g. 10.1016/j.fluid.2008.02.016).  Anything else
    (a second crystal, a liquid crystal, no liquid) is not fusion."""
    low = [p.lower() for p in phases]
    n_cr = sum(p.startswith("crystal") for p in low)
    n_liq = sum(p == "liquid" for p in low)
    other = [p for p in low if not (p.startswith("crystal") or p == "liquid" or p in ("air at 1 atmosphere", "gas"))]
    return n_cr == 1 and n_liq == 1 and not other


def parse_file(root, fname: str, stats: collections.Counter, max_p: float = 1000.0) -> list[dict]:
    doi = first_text(root, "sDOI")
    compounds = {first_text(c, "nOrgNum"): dict(inchi=first_text(c, "sStandardInChI"), name=first_text(c, "sCommonName"))
                 for c in kids(root, "Compound")}
    out = []
    for d in kids(root, "PureOrMixtureData"):
        ids = [first_text(c, "nOrgNum") for c in kids(d, "Component")]
        props = {}
        for p in kids(d, "Property"):
            pn = first_text(p, "ePropName") or ""
            if pn in PROPS:
                props[first_text(p, "nPropNumber")] = PROPS[pn]
        if not props:
            continue
        if len(ids) != 1:
            stats["mixture dataset (not pure)"] += 1
            continue
        phases = _phases(d)
        consts = []
        for c in kids(d, "Constraint"):
            cid = kids(c, "ConstraintID")
            kind, text = type_of(cid[0]) if cid else (None, None)
            consts.append((kind, text, first_text(c, "nOrgNum"), to_float(first_text(c, "nConstraintValue"))))
        for nv in kids(d, "NumValues"):
            state, comp = {}, {}
            for kind, text, org, val in consts:
                _assign(state, comp, kind, text, org, val, None)
            for pv in kids(nv, "PropertyValue"):
                key = first_text(pv, "nPropNumber")
                if key not in props:
                    continue
                prop = props[key]
                val = to_float(first_text(pv, "nPropValue"))
                if val is None:
                    stats["no numeric value"] += 1
                    continue
                P = state.get("P")
                if P is not None and P > max_p:
                    stats["pressure > max"] += 1
                    continue
                if not is_fusion_phases(phases):           # e.g. a "triple point" or enthalpy between Crystal 2 and Crystal 1
                    stats[f"{prop} dataset not exactly one Crystal + Liquid (solid-solid or other, dropped)"] += 1
                    continue
                out.append(dict(file=fname, doi=doi, org=ids[0], prop=prop, value=val, P_kPa=P, phases="|".join(phases),
                                inchi=compounds.get(ids[0], {}).get("inchi"), name=compounds.get(ids[0], {}).get("name")))
    return out


def extract(archives, max_p: float = 1000.0):
    stats: collections.Counter = collections.Counter()
    rows = []
    for arch in archives:
        for fname, blob in iter_xml(Path(arch)):
            stats["files read"] += 1
            if not any(k.encode() in blob for k in PROPS):
                continue
            try:
                root = ET.parse(io.BytesIO(blob)).getroot()
                rows.extend(parse_file(root, fname, stats, max_p))
            except Exception:  # noqa: BLE001 - one odd file must not stop the run
                stats["files unparseable"] += 1
    out = []
    for r in rows:
        ik, mw, smi, hc = mol_info(r["inchi"])
        if ik is None:
            stats["no usable InChI"] += 1
            continue
        if not hc:
            stats["not a hydrocarbon"] += 1
            continue
        out.append({**{k: r[k] for k in ("file", "doi", "prop", "value", "P_kPa", "phases", "name")}, "inchikey": ik, "smiles": smi,
                    "source": TIER})
    df = pd.DataFrame(out)
    if len(df):
        n0 = len(df)
        df = df.drop_duplicates(["file", "inchikey", "prop", "value"]).reset_index(drop=True)
        stats["exact duplicates removed"] = n0 - len(df)
    return df, stats


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """Per compound: tm_K (median of tm and ttp rows), dhfus_kJ_mol (median), ranges, counts, spread flags."""
    recs = []
    for ik, g in df.groupby("inchikey"):
        t = g[g.prop.isin(["tm", "ttp"])].value
        h = g[g.prop == "dhfus"].value
        h = h[h > 0]
        rec = dict(inchikey=ik, smiles=g.smiles.iloc[0], name=g.name.dropna().iloc[0] if g.name.notna().any() else None,
                   tm_K=float(t.median()) if len(t) else np.nan, tm_min=float(t.min()) if len(t) else np.nan,
                   tm_max=float(t.max()) if len(t) else np.nan, n_tm=int(len(t)), n_tm_files=int(g[g.prop.isin(["tm", "ttp"])].file.nunique()),
                   dhfus_kJ_mol=float(h.median()) if len(h) else np.nan, dh_min=float(h.min()) if len(h) else np.nan,
                   dh_max=float(h.max()) if len(h) else np.nan, n_dh=int(len(h)), source=TIER)
        rec["spread_tm"] = bool(len(t) and rec["tm_max"] - rec["tm_min"] > 5.0)
        rec["spread_dh"] = bool(len(h) and (rec["dh_max"] - rec["dh_min"]) > 0.15 * rec["dhfus_kJ_mol"])
        recs.append(rec)
    return pd.DataFrame(recs)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("archives", nargs="*")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    archives = a.archives or [p for p in DEFAULT_ARCHIVES if Path(p).exists()]
    if not archives:
        sys.exit("no ThermoML archive found")
    df, stats = extract(archives)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    s = summarize(df)
    s.to_csv(str(a.out).replace(".csv", "_by_compound.csv"), index=False)
    print(f"wrote {a.out}: {len(df)} rows, {len(s)} compounds")
    for k, v in stats.most_common():
        print(f"  {v:>8}  {k}")


if __name__ == "__main__":
    main()
