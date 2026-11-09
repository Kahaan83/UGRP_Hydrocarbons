"""v4 Task 6: jet fuels with cited class-level compositions.

Outputs
  data/real_fuels/v4/jet_class_compositions.csv  one row per GC-composition fuel (class x carbon-number pseudo-components)
  data/real_fuels/v4/jet_properties_v4.csv       long table of property targets/descriptors, raw files only
                                                 (POSF-10325 from av2718; the 57 inventory_v4 CRC647_* fuels from crc647)
  results/v4/P2/jet_corpus.md                    count summary
  data/MISSING.md                                block 'v4-task6' (idempotent): properties with no local source (-> Task 8)
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.loaders_jet import TARGETS, load_gc_fuels, load_jet_properties  # noqa: E402

GC = ROOT / "data" / "raw" / "gc_composition" / "fuel_gc_compositions.csv"
AV = ROOT / "data" / "raw" / "pour_point" / "av2718_fuel_properties.csv"
CRC = ROOT / "data" / "raw" / "pour_point" / "crc647_fuel_properties.csv"
V4 = ROOT / "data" / "real_fuels" / "v4"
OUT_C, OUT_P = V4 / "jet_class_compositions.csv", V4 / "jet_properties_v4.csv"
REPORT = ROOT / "results" / "v4" / "P2" / "jet_corpus.md"
MISSING = ROOT / "data" / "MISSING.md"
TARGET_NAMES = [("net_heat_of_combustion_MJ_kg" if t == "net_heat_content_btu_lb" else t) for t in TARGETS]
BEGIN, END = "<!-- v4-task6:begin -->", "<!-- v4-task6:end -->"

DOC = {
    "POSF-10325": "CRC AV-27-18 lists POSF 10325 but not freezing/pour/heat of combustion; fetch the AFRL/CRC property sheet "
                  "for POSF 10325 (e.g. Edwards, 'Reference Jet Fuels for Combustion Testing', AIAA 2017-0146, or the "
                  "ASCENT 31a deck's property slide). AJFTD fuel_id 2 (label '10325': freezing -50 C, D86) is NOT linked: "
                  "its identity with POSF 10325 is unconfirmed and its flash (48 C) disagrees with AV-27-18 (55 C)",
    "POSF-13511": "no row in av2718/crc647/AJFTD; fetch POSF 13511 property sheet (Edwards AIAA 2017-0146 / AFRL "
                  "'Alternative Jet Fuels' tables) and a D86 or D2887 curve (needed for the carbon-number split)",
    "POSF-13516": "no row in av2718/crc647/AJFTD; fetch POSF 13516 property sheet (Edwards AIAA 2017-0146 / AFRL) "
                  "and a D86 or D2887 curve (needed for the carbon-number split)",
    "UDRI_HTA_repfuel": "UDRI 'HTA representative fuel' (ASCENT 31a deck slide 12) has no property rows locally; fetch the "
                        "fuel's specification-property table from the same ASCENT 31a report and a D86/D2887 curve",
    "JetA_API14avg": "API 2010 'Jet A' 14-sample average: fetch API (2010b) Table of fuel properties cited by "
                     "ATSDR TP-121 Table 4-2, and a mean D86 curve",
}


def _targets_present(props: pd.DataFrame, fid: str) -> set[str]:
    return set(props[(props.fuel_id == fid) & (props.role == "target")].property)


def write_missing(comp: pd.DataFrame, props: pd.DataFrame) -> int:
    rows = []
    for _, r in comp.iterrows():
        fid = r.fuel_id
        if fid in DOC:
            miss = [t for t in TARGET_NAMES if t not in _targets_present(props, fid)]
            cn = "; carbon-number split assumed (no same-fuel distillation)" if "assumed_cn_split" in r.headline_reason else ""
            rows.append(f"| 2026-11-10 | v4 Task 6: {fid} property targets | missing: {', '.join(miss)}{cn} | {DOC[fid]} |")
    rows.append("| 2026-11-10 | v4 Task 6: JP5_PotterSimmons / JP8_PotterSimmons composition | ATSDR TP-121 Tables 3/4 list identified "
                "compounds only (JP-5 40.9 %, JP-8 42.3 % of mass); rows kept as `excluded_partial_composition`, never normalised | "
                "fetch Potter & Simmons (1998) full GC tables (or a DHA/GCxGC type analysis) for the same two fuels |")
    rows.append("| 2026-11-10 | v4 Task 6: inventory_v4 RP-1 fuels (NIST_RP1, RP1_Binary_2comp, RP1_Edwards_3comp, RP1_NIST_14comp) property targets | "
                "no linked measured properties; `rp1_properties.csv` (NISTIR 6646) describes 'RP-1 original sample', whose identity with "
                "SRM 4572 (`NIST_RP1`) is unconfirmed; surrogates have no measurements | fetch the SRM 4572 certificate / Fortin et al. 2012 property table |")
    block = [BEGIN, "", "### v4 Task 6: jet properties with no local source (2026-11-10)", "",
             "| Date Logged | Item / Source | Expected Content | Current Status & Mitigation |", "|---|---|---|---|",
             *rows, "", END]
    text = MISSING.read_text(encoding="utf-8") if MISSING.exists() else ""
    text = re.sub(rf"\n*{re.escape(BEGIN)}.*?{re.escape(END)}\n*", "\n", text, flags=re.S).rstrip() + "\n\n" + "\n".join(block) + "\n"
    MISSING.write_text(text, encoding="utf-8")
    return len(rows)


def main() -> None:
    comp = load_gc_fuels(GC, AV)
    props = load_jet_properties(AV, CRC)
    inv = pd.read_csv(V4 / "inventory_v4.csv")
    inv_ids = set(inv.fuel_id)
    # keep only fuels that exist in the v4 corpus: the GC fuels and the inventory_v4 CRC647 fuels
    props = props[props.fuel_id.isin(set(comp.fuel_id) | inv_ids)].reset_index(drop=True)
    assert not set(comp.fuel_id) & inv_ids, "GC fuel ids must not duplicate inventory_v4 ids"
    comp["n_target_props"] = comp.fuel_id.map(lambda f: len(_targets_present(props, f)))
    comp.to_csv(OUT_C, index=False)
    props.to_csv(OUT_P, index=False)
    n_missing = write_missing(comp, props)

    crc_ids = sorted(f for f in inv_ids if f.startswith("CRC647_"))
    have = sorted(set(props.fuel_id) & set(crc_ids))
    L = ["# Jet corpus (v4 Task 6)", "",
         f"Built by `scripts/v4/build_jet_corpus.py` (loader `src/fuelblend/data/loaders_jet.py`).  Outputs: "
         f"`{OUT_C.relative_to(ROOT).as_posix()}` ({len(comp)} fuels), `{OUT_P.relative_to(ROOT).as_posix()}` ({len(props)} rows).", "",
         "## Fuels (class x carbon-number compositions)", "",
         f"* fuels: {len(comp)}; composition_complete: {int(comp.composition_complete.sum())}; "
         f"excluded_partial_composition: {int((~comp.composition_complete).sum())}",
         f"* **headline_ok: {int(comp.headline_ok.sum())}**; not headline: {int((~comp.headline_ok).sum())} "
         f"(assumed_cn_split: {int(comp.headline_reason.str.contains('assumed_cn_split').sum())}, "
         f"partial composition: {int((~comp.composition_complete).sum())})",
         f"* link to inventory_v4: none of the {len(comp)} GC fuels is in inventory_v4 (those POSF ids were quarantined in Task 3); "
         f"ids are disjoint (asserted), `inventory_fuel_id` is empty.  The 57 `CRC647_*` inventory fuels are linked on the property side only.", "",
         "| fuel_id | complete | cn_split_method | headline_ok | pseudo-components | target props | reason |", "|---|---|---|---|---|---|---|"]
    for _, r in comp.iterrows():
        L.append(f"| {r.fuel_id} | {r.composition_complete} | {r.cn_split_method} | {r.headline_ok} | "
                 f"{r.n_pseudo_components} | {r.n_target_props} | {r.headline_reason} |")
    L += ["", "## Properties present per source (raw files only)", "",
          "| source | fuels | rows | target rows | targets present (fuels with it) |", "|---|---|---|---|---|"]
    for src, g in props.groupby("source_file"):
        t = g[g.role == "target"]
        per = ", ".join(f"{p} ({t[t.property == p].fuel_id.nunique()})" for p in sorted(t.property.unique()))
        L.append(f"| {src} | {g.fuel_id.nunique()} | {len(g)} | {len(t)} | {per} |")
    L += ["", f"* inventory_v4 CRC647 fuels with properties: {len(have)} of {len(crc_ids)}",
          f"* `data/MISSING.md` rows added (block `v4-task6`): {n_missing}", "",
          "## Notes / limits", "",
          "* The raw GC file lists class totals AND sub-totals (POSF-10325 raw rows sum to 152.28 %); the loader sums only a "
          "named non-overlapping partition per fuel (POSF-10325: 0.9999).",
          "* POSF-10325's carbon-number split comes from the same fuel's av2718 ASTM **D2887** simulated distillation (no D86 "
          "exists for it in a raw file) through `fuel_hindcast.carbon_distribution`; that function applies ONE carbon-number "
          "distribution to every class (class-independent), which is an approximation.  Column `cn_split_method` records it.",
          "* Property values were never interpolated; nothing is NaN-filled (absence = no row).  Btu/lb -> MJ/kg is the exact "
          "factor 0.002326; the density `(60F)` unit label gives T = 15.56 C."]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L[:14]))


if __name__ == "__main__":
    main()
