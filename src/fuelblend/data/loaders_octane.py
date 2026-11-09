"""Gasoline-blend octane loaders (v4 Task 5): FACE set I/II and Co-Optima bRON.

Both loaders return one row per measured blend with columns
``blend_id, base_fuel, components, ron, mon, provenance, composition_provenance,
unresolved_members, composition_resolved, target_ok, target_ok_reason`` (plus source-specific
extras).  ``components`` is a JSON list of ``{name, inchikey, frac, basis}`` where ``inchikey`` is
a real InChIKey or ``LUMP:<id>`` for a pseudo-component whose molecular composition is unknown.

FACE CSV layout (``face_gasoline_set_I/II.csv``) -- decoded from the real files
-------------------------------------------------------------------------------
Header: ``Name,n-N,C_1,n-P,C_6,iso-P,C_2,C_5,C_4,C_3,FACE A gasoline,...,FACE J gasoline,C_7,RON,MON``.
It does NOT alternate component-name / fraction columns.  Every column between ``Name`` and ``RON``
is a *slot*; the header text of a slot is an opaque label (``C_1`` ...) and the CELL holds the
component NAME present in that blend (``n-Heptane``, ``iso-Octane``, ``Toluene``, ``Ethanol``,
``1,2,4-Trimethylbenzene``, ``Cyclopentane``, ``1-Hexene``, ``FACE A gasoline`` ...) or is empty.
The file contains NO fractions of any component.  Consequences:

* hydrocarbon members (PRF / TPRF / MC mixtures) have unknown ratios -> one unresolved lump
  ``LUMP:<blend>`` whose ``unresolved_members`` are the listed names;
* FACE base gasolines are the pseudo-components ``LUMP:FACE_A`` ... (no GC composition in the file);
* only for names ``PRF+N``, ``TPRF {1,2,3}+N``, ``TMB+N``, ``CP+N``, ``Hex+N`` and ``FACE X+N`` is N
  read as the ethanol VOLUME percentage (an inference from the naming convention, consistent with
  the Ethanol cell being filled iff N>0; it is not stated in the file, hence
  ``composition_provenance = "unknown"``).  ``TPRF 4+N`` is a series index, not a percentage, and
  is left unparsed.
* RON/MON are read directly from the file (``provenance = "raw_file"``).

Because every FACE row contains a lump, ``target_ok`` is False for all of them (the lump's own
RON/MON, where known, is a model input only).

Co-Optima bRON (``bron_data_with_cas.xls``, sheet ``Components``)
-----------------------------------------------------------------
One row per pure blendstock (CAS, NAME) with up to 5 data points per property, columns
``base_RON_bRON_i, vol_frac_bRON_i, blend_RON_bRON_i`` (same for MON): a base gasoline of RON
``base_RON`` is blended with ``vol_frac`` (a 0-1 VOLUME fraction) of the blendstock and the blend
RON/MON is measured.  The base gasoline has no composition -> ``LUMP:BRON_BASE_RON<base_RON>``,
carried as an input (``base_ron``/``base_mon``); ``target_ok`` is False under the lump rule while
``base_lump_input_only`` marks that the base RON/MON are inputs.  One output row per distinct
(blendstock, vol_frac) over the RON and MON data points.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

COLUMNS = ["blend_id", "base_fuel", "components", "ron", "mon", "provenance",
           "composition_provenance", "unresolved_members", "composition_resolved",
           "target_ok", "target_ok_reason", "source", "base_ron", "base_mon",
           "base_lump_input_only"]

# CAS of the pure components that appear by NAME in the FACE files (typed here; the InChIKey
# itself comes from the `chemicals` identifier database, never typed)
_ETHANOL_CAS = "64-17-5"
_ETHANOL_NAME = re.compile(r"^(PRF|TPRF [123]|TMB|CP|Hex|FACE [A-J])\+(\d+)$")


def _inchikey(cas: str) -> str | None:
    from chemicals.identifiers import search_chemical

    try:
        return search_chemical(cas).InChI_key
    except Exception:
        return None


def load_face_sets(paths: Iterable[str | Path]) -> pd.DataFrame:
    """Parse FACE set I/II CSVs (see module docstring). Rows are kept in file order."""
    rows: list[dict] = []
    for path in paths:
        raw = pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
        slot_cols = [c for c in raw.columns if c not in ("Name", "RON", "MON")]
        for _, r in raw.iterrows():
            name = r["Name"].strip()
            if not name:
                continue
            members = [r[c].strip() for c in slot_cols if r[c].strip()]
            face = [m for m in members if m.lower().startswith("face ")]
            eth_cell = any(m.lower() == "ethanol" for m in members)
            hc = [m for m in members if m.lower() != "ethanol" and m not in face]
            m = _ETHANOL_NAME.match(name)
            eth = int(m.group(2)) if m else None
            if eth is not None and (eth > 0) != eth_cell:
                eth = None                                 # name and Ethanol cell disagree
            if face:
                base = face[0].replace(" gasoline", "").replace(" ", "_")   # FACE_A
                lump_key, base_fuel, lump_members = f"LUMP:{base}", base, face
            else:
                lump_key, base_fuel, lump_members = f"LUMP:{name}", name, hc
            eth_frac = eth / 100.0 if eth is not None else None
            comps: list[dict] = []
            if eth_frac:
                comps.append({"name": "Ethanol", "inchikey": _inchikey(_ETHANOL_CAS),
                              "frac": eth_frac, "basis": "volume"})
            comps.append({"name": lump_key, "inchikey": lump_key,
                          "frac": 1.0 - (eth_frac or 0.0),
                          "basis": "volume" if eth_frac is not None else "unresolved"})
            rows.append({
                "blend_id": name, "base_fuel": base_fuel, "components": json.dumps(comps),
                "ron": pd.to_numeric(r["RON"], errors="coerce"),
                "mon": pd.to_numeric(r["MON"], errors="coerce"),
                "provenance": "raw_file", "composition_provenance": "unknown",
                "unresolved_members": json.dumps(lump_members),
                "composition_resolved": False, "target_ok": False,
                "target_ok_reason": "contains lump pseudo-component (no fractions in file)",
                "source": f"face:{Path(path).name}", "base_ron": float("nan"),
                "base_mon": float("nan"), "base_lump_input_only": True,
            })
    return pd.DataFrame(rows, columns=COLUMNS)


def load_cooptima_bron(path: str | Path) -> pd.DataFrame:
    """Parse the Co-Optima bRON workbook (see module docstring)."""
    raw = pd.read_excel(path)
    rows: list[dict] = []
    for _, r in raw.iterrows():
        cas = r["CAS"] if isinstance(r["CAS"], str) and r["CAS"].strip() else None
        name = str(r["NAME"]).strip()
        key = _inchikey(cas) if cas else None
        stock_key = key or "LUMP:" + re.sub(r"\W+", "_", name).strip("_").upper()
        pts: dict[float, dict] = {}
        for prop in ("RON", "MON"):
            for i in range(5):
                f = r.get(f"vol_frac_b{prop}_{i}")
                b = r.get(f"blend_{prop}_b{prop}_{i}")
                if pd.isna(f) or pd.isna(b):
                    continue
                d = pts.setdefault(round(float(f), 6), {"frac": float(f)})
                d[prop.lower()] = float(b)
                d[f"base_{prop.lower()}"] = float(r[f"base_{prop}_b{prop}_{i}"])
        for _, d in sorted(pts.items()):
            if "base_ron" in d:
                base_id = f"BRON_BASE_RON{d['base_ron']:g}"
            else:
                base_id = f"BRON_BASE_MON{d['base_mon']:g}"
            comps = [
                {"name": name, "inchikey": stock_key, "frac": d["frac"], "basis": "volume"},
                {"name": f"LUMP:{base_id}", "inchikey": f"LUMP:{base_id}",
                 "frac": 1.0 - d["frac"], "basis": "volume"},
            ]
            rows.append({
                "blend_id": f"BRON:{cas or stock_key}:{d['frac']:.4f}", "base_fuel": base_id,
                "components": json.dumps(comps), "ron": d.get("ron", float("nan")),
                "mon": d.get("mon", float("nan")), "provenance": "raw_file",
                "composition_provenance": "raw_file", "unresolved_members": json.dumps([]),
                "composition_resolved": False, "target_ok": False,
                "target_ok_reason": "base gasoline is a lump (RON/MON are inputs only)",
                "source": f"cooptima_bron:{Path(path).name}",
                "base_ron": d.get("base_ron", float("nan")),
                "base_mon": d.get("base_mon", float("nan")), "base_lump_input_only": True,
            })
    return pd.DataFrame(rows, columns=COLUMNS)
