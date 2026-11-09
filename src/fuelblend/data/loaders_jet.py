"""Jet-fuel class x carbon-number corpus loaders (v4 Task 6).

``load_gc_fuels``  : ``data/raw/gc_composition/fuel_gc_compositions.csv`` -> one row per fuel with
                     class fractions and class x carbon-number pseudo-components.
``load_jet_properties`` : property targets / descriptors, long format, ONLY from raw files
                     (``av2718_fuel_properties.csv`` for the POSF fuel, ``crc647_fuel_properties.csv`` for CRC 647).

Raw GC file layout (verified against the real file): ``fuel_code, fuel_name, reference, class_or_component,
frac_basis, pct, source_page_or_table, note``.  A fuel's rows are NOT a partition: the file lists class totals
AND their sub-totals/sub-bins (POSF-10325 rows sum to 152.28 %).  ``PARTITION`` below names, per fuel, the
non-overlapping rows that make up the classes (checked against the sum), and the remaining carbon-number rows
(``..._C<n>``) are the carbon-number bins of a class.

Carbon-number rule (task brief): a class's carbon-number split is
  * ``source``          - bins given by the cited table (scaled to the class total if they cover >= 98 % of it);
  * ``D86/D2887_same_fuel`` - otherwise, the measured boiling curve of THE SAME fuel in a raw file
                          (class-independent n-alkane-Tb mapping = ``fuel_hindcast.carbon_distribution``);
  NOTE: the distillation rows of such a fuel (POSF-10325) are role=descriptor in the property table AND drive its
  carbon-number split; do not feed them to a model as features alongside the derived composition (leakage).
  * ``assumed_cn_split`` - otherwise (cn = None, never invented).  Such a fuel is ``headline_ok = False``.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.mixture.fuel_hindcast import carbon_distribution

# non-overlapping rows that make each fuel's classes (names in column class_or_component)
_CYC3 = ["monocycloparaffins_total", "dicycloparaffins_total", "tricycloparaffins_total"]
PARTITION: dict[str, dict[str, list[str]] | str] = {
    "JetA_API14avg": {"paraffin_unsplit": ["paraffins_total"], "cycloalkane": _CYC3,
                      "aromatic": ["aromatics_total"], "alkene_alkyne": ["olefins_total"]},
    "POSF-10325": {"aromatic": ["aromatics_total"], "iso_paraffin": ["iso_paraffins_total"],
                   "n_paraffin": ["n_paraffins_total"], "cycloalkane": _CYC3},
    "POSF-13511": {"aromatic": ["aromatics_alkylbenzenes_total", "aromatics_cycloaromatics_total"],
                   "iso_paraffin": ["iso_paraffins_total"], "n_paraffin": ["n_paraffins_total"],
                   "cycloalkane": _CYC3},
    "POSF-13516": {"iso_paraffin": ["iso_paraffins_total"], "n_paraffin": ["n_paraffins_total"],
                   "cycloalkane": _CYC3},
    "UDRI_HTA_repfuel": {"iso_paraffin": ["iso_paraffins_total"], "n_paraffin": ["n_paraffins_total"],
                         "cycloalkane": ["cycloparaffins_total"],
                         "aromatic": ["aromatics_alkylbenzenes_total", "aromatics_alkylnaphthalenes_total",
                                      "aromatics_indanes_tetralins_total"]},
    # ATSDR Tables 3/4: a list of identified compounds only (~41 % of the fuel); every row is a CN bin
    "JP5_PotterSimmons": "leaf",
    "JP8_PotterSimmons": "leaf",
}
FUEL_NAMES_AV2718 = {"POSF-10325": 10325}        # GC fuel -> fuel_code in av2718 (POSF number = CRC AV-27-18 code)
_CN = re.compile(r"_C(\d+)(?:_(plus|and_lower))?$")
_TOL = 0.02

COLUMNS = ["fuel_id", "fuel_name", "source_ref", "basis", "class_fractions", "components", "n_pseudo_components",
           "class_sum_raw", "unresolved_fraction", "composition_complete", "corpus_status", "cn_split_method",
           "cn_assumed_classes", "cn_split_class_independent", "headline_ok", "headline_reason", "composition_provenance",
           "property_provenance", "inventory_fuel_id", "n_target_props"]


def _class_of(name: str) -> str | None:
    if name.startswith("n_paraffins"):
        return "n_paraffin"
    if name.startswith("iso_paraffins"):
        return "iso_paraffin"
    if name.startswith("paraffins"):
        return "paraffin_unsplit"
    if name.startswith(("cycloparaffins", "monocyclo", "dicyclo", "tricyclo")):
        return "cycloalkane"
    if name.startswith("aromatics"):
        return "aromatic"
    if name.startswith("olefins"):
        return "alkene_alkyne"
    return None


def _distillation_cn(av: pd.DataFrame, code: int) -> tuple[dict[int, float], str] | None:
    """Carbon-number mass distribution from the measured boiling curve of fuel ``code`` (D86 preferred)."""
    sub = av[(av.fuel_code == code) & av.property.str.fullmatch(r"distillation_T\d+_C")]
    if sub.empty:
        return None
    meth = "D86" if sub.method.astype(str).str.contains("D86").any() else "D2887"
    sub = sub[sub.method.astype(str).str.contains(meth)]
    keymap = {0: "ibp", 100: "fbp"}
    row = {}
    for _, r in sub.iterrows():
        pct = int(re.search(r"T(\d+)_C", r.property).group(1))
        row[f"d86_{keymap.get(pct, f't{pct}')}"] = float(r.value_std)
    cd = carbon_distribution(pd.Series(row))
    if not cd:
        return None
    cd = {int(n): f for n, f in cd.items() if f >= 1e-3}          # drop numerical dust
    s = sum(cd.values())
    return {n: f / s for n, f in cd.items()}, f"{meth}_same_fuel"


def load_gc_fuels(gc_csv: str | Path, av2718_csv: str | Path) -> pd.DataFrame:
    gc = pd.read_csv(gc_csv)
    av = pd.read_csv(av2718_csv)
    rows = []
    for code, g in gc.groupby("fuel_code", sort=False):
        if code not in PARTITION:
            continue
        spec = PARTITION[code]
        assert (g.frac_basis == "mass").all(), code
        name, ref = g.fuel_name.iloc[0], g.reference.iloc[0]
        pct_of = g.groupby("class_or_component").pct.sum()
        cls_pct: dict[str, float] = {}
        bins: dict[str, dict[int, tuple[float, str]]] = {}
        if spec == "leaf":
            for comp, p in pct_of.items():
                c, m = _class_of(comp), _CN.search(comp)
                cls_pct[c] = cls_pct.get(c, 0.0) + p
                bins.setdefault(c, {})
                n = int(m.group(1))
                bins[c][n] = (bins[c].get(n, (0.0, ""))[0] + p, comp)
            part_rows = set(pct_of.index)
        else:
            part_rows = {r for rs in spec.values() for r in rs}
            for c, rs in spec.items():
                miss = [r for r in rs if r not in pct_of.index]
                assert not miss, (code, miss)
                cls_pct[c] = float(sum(pct_of[r] for r in rs))
            for comp, p in pct_of.items():
                m = _CN.search(comp)
                c = _class_of(comp)
                if comp in part_rows or not m or c not in spec:
                    continue
                n = int(m.group(1))
                bins.setdefault(c, {})
                bins[c][n] = (bins[c].get(n, (0.0, ""))[0] + p, comp)
        total = sum(cls_pct.values()) / 100.0
        complete = abs(total - 1.0) <= _TOL
        scale = 1.0 / total if complete else 1.0       # normalise only a complete composition
        dist = _distillation_cn(av, FUEL_NAMES_AV2718[code]) if code in FUEL_NAMES_AV2718 else None
        comps, assumed = [], []
        for c, tp in cls_pct.items():
            K = sum(p for p, _ in bins.get(c, {}).values())
            parts: list[tuple[int | None, float, str, str, str | None]] = []   # cn, pct, cn_source, prov, label
            if K > 0:
                f = tp / K if K >= 0.98 * tp else 1.0
                for n, (p, lab) in sorted(bins[c].items()):
                    parts.append((n, p * f, "source", "cited_typed", lab))
            rest = tp - sum(p for _, p, *_ in parts)
            if rest > 1e-9 and tp > 0:
                if dist is not None:
                    for n, w in sorted(dist[0].items()):
                        parts.append((n, rest * w, dist[1], "formula_derived", None))
                else:
                    parts.append((None, rest, "assumed_cn_split", "cited_typed", None))
                    assumed.append(c)
            for n, p, src, prov, lab in parts:
                comps.append({"class": c, "cn": n, "cn_label": lab, "frac": round(p / 100.0 * scale, 6),
                              "basis": "mass", "cn_source": src, "provenance": prov})
        cfrac = {c: round(p / 100.0 * scale, 6) for c, p in cls_pct.items()}
        srcs = {c["cn_source"] for c in comps}
        method = ("assumed_cn_split" if assumed else "+".join(sorted(srcs - {"source"})) or "source")
        reasons = []
        indep = any(c["cn_source"].endswith("_same_fuel") for c in comps)
        if indep:
            reasons.append("cn_split_class_independent: one n-alkane-Tb-mapped distillation distribution applied to every "
                           "class (aromatics/cycloalkanes get the n-paraffin shape); vector partly formula_derived; "
                           "descriptor distillation rows also drive it (leakage if used as features)")
        if not complete:
            reasons.append(f"partial_composition (listed {total:.3f} of the fuel; unresolved {1 - total:.3f})")
        if assumed:
            reasons.append(f"assumed_cn_split: {sorted(set(assumed))}")
        rows.append({
            "fuel_id": code, "fuel_name": name, "source_ref": ref, "basis": "mass",
            "class_fractions": json.dumps(cfrac), "components": json.dumps(comps),
            "n_pseudo_components": len(comps), "class_sum_raw": round(total, 6),
            "unresolved_fraction": 0.0 if complete else round(1.0 - total, 6),
            "composition_complete": complete,
            "corpus_status": "included" if complete else "excluded_partial_composition",
            "cn_split_method": method, "cn_assumed_classes": json.dumps(sorted(set(assumed))),
            "cn_split_class_independent": indep,
            "headline_ok": complete and not assumed, "headline_reason": "; ".join(reasons),
            "composition_provenance": "cited_typed", "property_provenance": "raw_file",
            "inventory_fuel_id": np.nan, "n_target_props": 0})
    return pd.DataFrame(rows, columns=COLUMNS)


def included_jet_fuel_ids(path: str | Path | pd.DataFrame) -> list[str]:
    """Fuel ids with corpus_status == "included".  Task 7 splits MUST use this, not every row of the CSV:
    the ATSDR JP-5/JP-8 rows (excluded_partial_composition) are not compositions."""
    df = path if isinstance(path, pd.DataFrame) else pd.read_csv(path)
    return df.loc[df.corpus_status == "included", "fuel_id"].tolist()


# --------------------------------------------------------------------------- properties
BTU_LB_TO_MJ_KG = 0.002326          # 1 Btu/lb = 2.326 J/g exactly
TARGETS = ("density_kg_m3", "viscosity_kin_cSt", "flash_point_C", "freezing_point_C", "pour_point_C",
           "net_heat_content_btu_lb")
DESCRIPTORS = ("aromatics_vol_pct", "olefins_vol_pct", "saturates_vol_pct", "naphthalenes_vol_pct",
               "hydrogen_wt_pct", "d86_residue_pct")
PROP_COLS = ["fuel_id", "property", "value", "unit", "T_C", "method", "role", "source_file", "source_page",
             "provenance", "raw_value", "raw_unit"]


def _long(df: pd.DataFrame, fuel_id, src: str) -> pd.DataFrame:
    keep = df[df.property.isin(TARGETS + DESCRIPTORS) | df.property.str.fullmatch(r"distillation_T\d+_C")].copy()
    keep = keep[pd.to_numeric(keep.value_std, errors="coerce").notna()]
    out = []
    for _, r in keep.iterrows():
        prop, val, unit = r.property, float(r.value_std), r.unit_std
        if prop == "net_heat_content_btu_lb":
            prop, val, unit = "net_heat_of_combustion_MJ_kg", val * BTU_LB_TO_MJ_KG, "MJ/kg"
        T = r.T_C
        if (pd.isna(T) and prop == "density_kg_m3" and isinstance(r.unit_raw, str) and "60F" in r.unit_raw):
            T = 15.56                                        # unit_raw says "g/cm3 (60F)"
        out.append({"fuel_id": fuel_id(r), "property": prop, "value": val, "unit": unit, "T_C": T,
                    "method": r.method if isinstance(r.method, str) else "",
                    "role": "target" if r.property in TARGETS else "descriptor", "source_file": src,
                    "source_page": r.get("page", np.nan), "provenance": "raw_file",
                    "raw_value": r.value_raw, "raw_unit": r.unit_raw})
    return pd.DataFrame(out, columns=PROP_COLS)


def load_jet_properties(av2718_csv: str | Path, crc647_csv: str | Path) -> pd.DataFrame:
    """Long property table; every value is a row of a raw file (nothing derived except unit conversion)."""
    av = pd.read_csv(av2718_csv)
    crc = pd.read_csv(crc647_csv)
    parts = []
    for fid, code in FUEL_NAMES_AV2718.items():
        parts.append(_long(av[av.fuel_code == code], lambda r, f=fid: f, "av2718_fuel_properties.csv"))
    parts.append(_long(crc, lambda r: f"CRC647_{int(r.fuel_code)}", "crc647_fuel_properties.csv"))
    return pd.concat(parts, ignore_index=True)
