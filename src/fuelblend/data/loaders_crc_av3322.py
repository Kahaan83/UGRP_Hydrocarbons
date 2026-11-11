"""CRC Report AV-33-22 (World Jet Fuel Survey Part 1, 2023-2024; UDRI, UDR-TR-2025-39) supplementary-workbook loader (v5).

Reads ``data/raw/crc_av3322/CRC-AV-33-22_Supplementary-Data.xlsx`` (sheet ``Data`` only; ``README`` and ``Bulk Modulus`` are
ignored; see that folder's PROVENANCE.md for the legal notice) with openpyxl and writes, to ``data/real_fuels/v5/``:

* ``crc_av3322_composition.csv``  long: one row per (fuel, UDRI FC-M-101 hydrocarbon-type bin), mass %, GCxGC tier A.
  Bins are class x carbon number (no molecule list); ``plus`` marks open-ended bins (C21+, C16+ ...), ``lower`` the
  'n-alkane-C03 and lower' bin.
* ``crc_av3322_properties.csv``   long: UDRI-measured rows (``source_tag = 'udri_measured'``) and supplier certificate-of-analysis
  rows (``'supplier_coa'``, the workbook's ``Specification`` rows, stored as ``<name>_method`` / ``<name>_value`` pairs).
  Censored CoA entries such as ``'<-60'`` keep ``qualifier = '<'`` and the numeric bound in ``value``; they are never coerced.
* ``crc_av3322_bulk.csv``         wide, one row per fuel.
* ``crc_av3322_build_log.json``.

The sheet is transposed: rows = properties, columns = CRC sample numbers 1..38.  Property names follow the AJFTD / FuelLib
loaders (``flash_point_C``, ``viscosity_kin_cSt``, ``density_kg_m3``, ``vapor_pressure_kPa``, ``freezing_point_C``,
``heat_of_combustion_net_MJ_kg``, ``distillation_T{pct}_C`` ...).  D2887 simulated distillation is stored as
``distillation_simdis_T{pct}_C`` where the sheet label '0.05% recovered_SimDis' means 5 % (fraction x 100); IBP -> T0, FBP -> T100.
The D86-*correlation* rows (derived from SimDis) are NOT loaded.

Double counting (ASSUMPTION, not verified): these 38 fuels are new 2023-2024 survey samples, independent of FuelLib/NJFCP
POSF fuels, AJFTD (2015-2016 tests) and CRC AV-27-18 / Report 647.  ``double_count_group`` is ``'CRC_AV3322'`` for all rows.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "data" / "raw" / "crc_av3322" / "CRC-AV-33-22_Supplementary-Data.xlsx"
OUT = REPO / "data" / "real_fuels" / "v5"
SHEET = "Data"
PROVENANCE = "crc_av3322_gcxgc_tierA"
DOUBLE_COUNT_GROUP = "CRC_AV3322"
PSI_TO_KPA = 6.894757293168
N_META = 9                      # metadata columns before the first sample column
NOT_A_VALUE = {"", "none", "nan", "n/a", "na", "-", "not specified", "no test"}

COMP_COLS = ["fuel_id", "crc_sample_no", "fuel_type", "group", "carbon_number", "plus", "lower", "mass_pct", "bin_label",
             "method", "provenance", "data_tier", "double_count_group"]
PROP_COLS = ["fuel_id", "crc_sample_no", "fuel_type", "property", "value", "unit_std", "T_K", "qualifier", "method",
             "source_tag", "row_label", "raw_value", "raw_unit", "provenance", "data_tier", "double_count_group"]
GROUPS = {"n-alkane": "n-alkane", "iso-alkane": "iso-alkane", "alkene": "alkene", "monocycloalkane": "monocycloalkane",
          "dicycloalkane": "dicycloalkane", "tricycloalkane": "tricycloalkane", "alkylbenzene": "alkylbenzene",
          "cycloaromatic": "cycloaromatic", "diaromatic": "diaromatic"}
NON_BIN_COMPOSITION = {"Hydrogen content", "Average Molecular Wt", "Average Molecular Formula - Carbon",
                       "Average Molecular Formula - Hydrogen"}
BIN_RE = re.compile(r"^(?P<g>[A-Za-z-]+)-C(?P<n>\d+)(?P<s>\+| and lower)?$")

# supplier CoA rows: sheet name (before _method/_value) -> (property, unit_std)
COA = {"Flash point": ("flash_point_C", "C"), "Freezing point": ("freezing_point_C", "C"),
       "Viscosity": ("viscosity_kin_cSt", "cSt"), "Net Heat of Combustion": ("heat_of_combustion_net_MJ_kg", "MJ/kg"),
       "Smoke point": ("smoke_point_mm", "mm"), "Density": ("density_kg_m3", "kg/m3"),
       "Aromatics": ("aromatics_vol_pct", "vol%"), "Naphthalenes": ("naphthalenes_vol_pct", "vol%")}
# UDRI measured, keyed by (subcategory, unit-in-sheet) -> (property, unit_std, multiplier)
MEASURED = {("Viscosity", "mm^2/s"): ("viscosity_kin_cSt", "cSt", 1.0),
            ("Density", "kg/m^3"): ("density_kg_m3", "kg/m3", 1.0),
            ("Vapor Pressure", "psi"): ("vapor_pressure_kPa", "kPa", PSI_TO_KPA),
            ("Surface Tension", "mN/m"): ("surface_tension_N_m", "N/m", 1e-3),
            ("Specific Heat", "kJ/kg*K"): ("heat_capacity_kJ_kg_K", "kJ/kg/K", 1.0)}


def parse_value(raw) -> tuple[float, str]:
    """-> (number or NaN, qualifier '<' / '>' / '').  Strings like '<-60' keep the bound and the qualifier."""
    if raw is None:
        return math.nan, ""
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        return float(raw), ""
    s = str(raw).strip()
    if s.lower() in NOT_A_VALUE:
        return math.nan, ""
    q = s[0] if s[:1] in "<>" else ""
    try:
        return float(s.lstrip("<>").strip()), q
    except ValueError:
        return math.nan, q


def _num(x) -> float:
    return float(x) if isinstance(x, (int, float)) and not isinstance(x, bool) else math.nan


def read_sheet(path: Path = RAW) -> tuple[list[dict], list[dict]]:
    """-> (fuels, rows).  fuels = [{col, sample_no, fuel_id, fuel_type}], rows = one dict per property row."""
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    raw = [list(r) for r in wb[SHEET].iter_rows(values_only=True)]
    wb.close()
    hdr = raw[0]
    if str(hdr[0]).strip() != "Heading \\ CRC Sample No." or str(hdr[2]).strip() != "Category":
        raise ValueError("unexpected header row in CRC AV-33-22 workbook")
    cols = [j for j in range(N_META, len(hdr)) if isinstance(hdr[j], (int, float))]
    fuels = [{"col": j, "sample_no": int(hdr[j]), "fuel_id": f"CRC_AV3322_{int(hdr[j]):02d}",
              "fuel_type": str(raw[1][j]).strip()} for j in cols]
    if str(raw[1][0]).startswith("Property") is False:
        raise ValueError("fuel-type row not found at row 1")
    rows = []
    for r in raw[2:]:
        r = r + [None] * (len(hdr) - len(r))
        if r[0] is None:
            continue
        rows.append({"name": str(r[0]).strip(), "unit": str(r[1]).strip(), "category": r[2], "sub": r[3], "method": r[4],
                     "temp": r[8], "values": {f["fuel_id"]: r[f["col"]] for f in fuels}})
    return fuels, rows


def parse_bin(label: str) -> dict:
    m = BIN_RE.match(label.strip())
    if not m or m["g"].lower() not in GROUPS:
        raise ValueError(f"unrecognised composition bin label {label!r}")
    s = m["s"] or ""
    return {"group": GROUPS[m["g"].lower()], "carbon_number": int(m["n"]), "plus": s == "+", "lower": s.strip() == "and lower"}


def build_compositions(fuels, rows, log: dict | None = None) -> pd.DataFrame:
    log = {} if log is None else log
    bins = [r for r in rows if r["sub"] == "Hydrocarbon Type" and r["unit"] == "mass%" and r["name"] not in NON_BIN_COMPOSITION]
    out = []
    for r in bins:
        b = parse_bin(r["name"])
        for f in fuels:
            v = r["values"][f["fuel_id"]]
            if not isinstance(v, (int, float)) or isinstance(v, bool) or math.isnan(v):
                raise ValueError(f"{f['fuel_id']} {r['name']}: non-numeric composition value {v!r}")
            out.append({"fuel_id": f["fuel_id"], "crc_sample_no": f["sample_no"], "fuel_type": f["fuel_type"], **b,
                        "mass_pct": float(v), "bin_label": r["name"], "method": r["method"], "provenance": PROVENANCE,
                        "data_tier": "A", "double_count_group": DOUBLE_COUNT_GROUP})
    df = pd.DataFrame(out, columns=COMP_COLS)
    sums = df.groupby("fuel_id")["mass_pct"].sum()
    log.update(n_fuels=len(fuels), n_bins=len(bins), composition_sum_min=round(float(sums.min()), 3),
               composition_sum_max=round(float(sums.max()), 3),
               excluded_non_bin_rows=sorted(NON_BIN_COMPOSITION))
    return df


def _row(f, prop, val, unit, T_K, q, method, tag, label, raw_val, raw_unit):
    return {"fuel_id": f["fuel_id"], "crc_sample_no": f["sample_no"], "fuel_type": f["fuel_type"], "property": prop,
            "value": val, "unit_std": unit, "T_K": T_K, "qualifier": q, "method": method, "source_tag": tag,
            "row_label": label, "raw_value": raw_val, "raw_unit": raw_unit, "provenance": PROVENANCE, "data_tier": "A",
            "double_count_group": DOUBLE_COUNT_GROUP}


def _tk(temp) -> float:
    return float(temp) + 273.15 if isinstance(temp, (int, float)) and not isinstance(temp, bool) else math.nan


def _distillation_name(label: str, tag: str, prefix: str) -> str | None:
    """'10% recovered_D86' / '0.05% recovered_SimDis' / 'Initial boiling point_X' -> distillation_T{pct}_C name."""
    base = label.rsplit("_", 1)[0]
    if base.startswith("Initial boiling point"):
        pct = 0
    elif base.startswith("Final boiling point"):
        pct = 100
    else:
        m = re.match(r"^([\d.]+)% recovered$", base)
        if not m:
            return None
        pct = float(m[1]) * (100.0 if tag == "simdis" else 1.0)
        pct = int(round(pct))
    return f"{prefix}_T{pct}_C"


def build_properties(fuels, rows, log: dict | None = None) -> pd.DataFrame:
    log = {} if log is None else log
    out, dropped = [], []
    for r in rows:
        sub, nm, unit = r["sub"], r["name"], r["unit"]
        if r["category"] == "Specification":
            continue
        pname = None
        mult, ustd = 1.0, None
        if (sub, unit) in MEASURED:
            pname, ustd, mult = MEASURED[(sub, unit)]
            T_K = _tk(r["temp"])
            if math.isnan(T_K):
                raise ValueError(f"{nm}: missing measurement temperature")
        elif sub == "Distillation" and r["method"] == "D86":
            T_K = math.nan
            if unit == "Celsius":
                pname, ustd = _distillation_name(nm, "d86", "distillation"), "C"
            elif nm.startswith("Residue"):
                pname, ustd = "distillation_residue_pct", "%"
            elif nm.startswith("Loss"):
                pname, ustd = "distillation_loss_pct", "%"
        elif sub == "Simulated Distillation" and nm.endswith("_SimDis"):
            T_K, ustd = math.nan, "C"
            pname = _distillation_name(nm, "simdis", "distillation_simdis")
        elif sub == "Hydrogen" and nm.lower().startswith("hydrogen"):
            T_K, pname, ustd = math.nan, "hydrogen_content_D3701_mass_pct", "mass%"
        if pname is None:
            continue
        for f in fuels:
            raw_v = r["values"][f["fuel_id"]]
            v, q = parse_value(raw_v)
            if math.isnan(v):
                dropped.append((f["fuel_id"], nm))
                continue
            out.append(_row(f, pname, v * mult, ustd, T_K, q, r["method"], "udri_measured", nm, v, unit))

    # supplier CoA: '<name>_method' row followed by '<name>_value' row (matched by name, not position)
    by_name = {r["name"]: r for r in rows if r["category"] == "Specification"}
    for base, (pname, ustd) in COA.items():
        rv, rm = by_name.get(f"{base}_value"), by_name.get(f"{base}_method")
        if rv is None or rm is None:
            raise ValueError(f"CoA pair missing for {base!r}")
        T_K = _tk(rv["temp"])        # 'N/A' -> NaN (viscosity -20 C and density 15 C carry a temperature)
        for f in fuels:
            raw_v = rv["values"][f["fuel_id"]]
            v, q = parse_value(raw_v)
            if math.isnan(v):
                dropped.append((f["fuel_id"], rv["name"]))
                continue
            meth = rm["values"][f["fuel_id"]]
            out.append(_row(f, pname, v, ustd, T_K, q, None if meth is None else str(meth).strip(), "supplier_coa",
                            rv["name"], str(raw_v), rv["unit"]))
    log["n_property_rows"] = len(out)
    log["properties_dropped_non_numeric"] = len(dropped)
    log["censored_rows"] = [(o["fuel_id"], o["property"], o["raw_value"]) for o in out if o["qualifier"]]
    return pd.DataFrame(out, columns=PROP_COLS)


def build_bulk(fuels, rows, comp: pd.DataFrame, prop: pd.DataFrame) -> pd.DataFrame:
    by = {r["name"]: r for r in rows if r["sub"] == "Hydrocarbon Type"}
    out = pd.DataFrame({"fuel_id": [f["fuel_id"] for f in fuels], "crc_sample_no": [f["sample_no"] for f in fuels],
                        "fuel_type": [f["fuel_type"] for f in fuels]})
    for col, nm in (("hydrogen_content_mass_pct", "Hydrogen content"), ("avg_mw_g_mol", "Average Molecular Wt"),
                    ("avg_formula_C", "Average Molecular Formula - Carbon"),
                    ("avg_formula_H", "Average Molecular Formula - Hydrogen")):
        out[col] = [_num(by[nm]["values"][f["fuel_id"]]) for f in fuels]
    out["composition_sum_pct"] = out["fuel_id"].map(comp.groupby("fuel_id")["mass_pct"].sum())
    coa = prop[prop.source_tag == "supplier_coa"]
    for col, p, T in (("coa_flash_point_C", "flash_point_C", None), ("coa_freezing_point_C", "freezing_point_C", None),
                      ("coa_viscosity_m20C_cSt", "viscosity_kin_cSt", None),
                      ("coa_net_heat_MJ_kg", "heat_of_combustion_net_MJ_kg", None), ("coa_smoke_point_mm", "smoke_point_mm", None),
                      ("coa_density_15C_kg_m3", "density_kg_m3", None), ("coa_aromatics_vol_pct", "aromatics_vol_pct", None)):
        s = coa[coa.property == p].set_index("fuel_id")
        out[col] = out["fuel_id"].map(s["value"])
        if (s["qualifier"] != "").any():
            out[col.replace("coa_", "coa_", 1) + "_qualifier"] = out["fuel_id"].map(s["qualifier"]).fillna("")
    out["provenance"], out["data_tier"], out["double_count_group"] = PROVENANCE, "A", DOUBLE_COUNT_GROUP
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log: dict = {}
    fuels, rows = read_sheet(RAW)
    comp = build_compositions(fuels, rows, log)
    prop = build_properties(fuels, rows, log)
    bulk = build_bulk(fuels, rows, comp, prop)
    log["fuel_types"] = {k: int(v) for k, v in bulk["fuel_type"].value_counts().items()}
    log["n_composition_rows"] = int(len(comp))
    comp.to_csv(OUT / "crc_av3322_composition.csv", index=False)
    prop.to_csv(OUT / "crc_av3322_properties.csv", index=False)
    bulk.to_csv(OUT / "crc_av3322_bulk.csv", index=False)
    (OUT / "crc_av3322_build_log.json").write_text(json.dumps(log, indent=1, default=str), encoding="utf-8")
    print(f"composition: {len(comp)} rows, {comp.fuel_id.nunique()} fuels; properties: {len(prop)} rows; bulk: {len(bulk)} rows")


if __name__ == "__main__":
    main()
