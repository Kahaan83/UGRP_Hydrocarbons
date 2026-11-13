"""Build data/raw/smoke_lit/*.csv from numbers transcribed out of the literature PDFs (see data/raw/smoke_lit/SOURCES.md).

Every number below was read from the cited table/page (checked against the page image where the text layer was doubtful).
Nothing is estimated.  Run: python scripts/v5/ingest/build_smoke_lit_raw.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "raw" / "smoke_lit"

SMILES = {
    "iso-cetane": "CC(C)CC(C)(C)CC(C)(C)CC(C)C",
    "1,2,4-trimethylbenzene": "Cc1ccc(C)c(C)c1",
    "1,3,5-trimethylbenzene": "Cc1cc(C)cc(C)c1",
    "iso-octane": "CC(C)CC(C)(C)C",
    "toluene": "Cc1ccccc1",
    "m-xylene": "Cc1cccc(C)c1",
    "n-dodecane": "CCCCCCCCCCCC",
    "1-methylnaphthalene": "Cc1cccc2ccccc12",
    "decalin": "C1CCC2CCCCC2C1",
    "methylcyclohexane": "CC1CCCCC1",
    "benzene": "c1ccccc1",
    "tetralin": "c1ccc2c(c1)CCCC2",
    "1,2-dihydronaphthalene": "C1Cc2ccccc2C=C1",
    "cyclohexylbenzene": "c1ccc(cc1)C1CCCCC1",
    "n-propylbenzene": "CCCc1ccccc1",
    "n-butylcyclohexane": "CCCCC1CCCCC1",
    "1,4-dimethylcyclooctane": "CC1CCC(C)CCCC1",
    "n-hexylbenzene": "CCCCCCc1ccccc1",
    "n-octane": "CCCCCCCC",
    "n-undecane": "CCCCCCCCCCC",
    "n-decane": "CCCCCCCCCC",
    "n-hexadecane": "CCCCCCCCCCCCCCCC",
    "n-heptane": "CCCCCCC",
    # isomer mix in the source; m-xylene used as the representative structure (see SOURCES.md)
    "xylenes": "Cc1cccc(C)c1",
}

COLS = ["source_id", "doi", "table_or_page", "sample_id", "kind", "components", "smiles", "fractions", "fraction_unit",
        "fraction_basis", "property", "value", "unit", "T_C", "method", "qualifier", "note"]
rows: list[dict] = []


def add(source_id, doi, table, sample, kind, comps, fracs, unit_f, basis, prop, value, unit, method="", qual="", note="", T_C=""):
    comps = list(comps)
    rows.append(dict(source_id=source_id, doi=doi, table_or_page=table, sample_id=sample, kind=kind,
                     components=json.dumps(comps), smiles=json.dumps([SMILES[c] for c in comps]),
                     fractions=json.dumps(list(fracs)), fraction_unit=unit_f if fracs else "",
                     fraction_basis=basis if fracs else "none",
                     property=prop, value=value, unit=unit, T_C=T_C, method=method, qualifier=qual, note=note))


def sp_tsi(source_id, doi, table, sample, kind, comps, fracs, unit_f, basis, sp, tsi, method, qual="", note_sp="", note_tsi=""):
    if sp is not None:
        add(source_id, doi, table, sample, kind, comps, fracs, unit_f, basis, "smoke_point_mm", sp, "mm", method, qual, note_sp)
    if tsi is not None:
        add(source_id, doi, table, sample, kind, comps, fracs, unit_f, basis, "tsi", tsi, "-", method, qual,
            note_tsi or note_sp)


# --------------------------------------------------------------------------------------------- S3 Mensch 2009 thesis
S3, S3D = "S3_mensch2009", ""
M = "smoke point lamp (Mensch thesis apparatus, ASTM D1322 type)"
# Table 4 (PDF p45, printed p34): pure compounds, MW (g/mol), SP (mm), TSI.
for name, mw, sp, tsi in [("cyclohexylbenzene", 160.26, 8.9, 69), ("1-methylnaphthalene", 142.20, 5.5, 100),
                          ("decalin", 138.26, 22.7, 20), ("tetralin", 132.21, 7.1, 71),
                          ("1,2-dihydronaphthalene", 130.19, 5.8, 87), ("1,2,4-trimethylbenzene", 120.20, 7.2, 63),
                          ("1,3,5-trimethylbenzene", 120.20, 7.3, 62), ("n-propylbenzene", 120.20, 8.5, 53),
                          ("iso-octane", 114.23, 40.0, 6.8), ("m-xylene", 106.17, 7.8, 51),
                          ("methylcyclohexane", 98.19, 40.8, 5.0), ("toluene", 92.14, 8.4, 40), ("benzene", 78.12, 9.1, 30)]:
    q = ""
    nsp = f"MW {mw}"
    if name == "iso-octane":
        q = "text_conflict"
        nsp += "; thesis text (p45) says pure iso-octane could not be measured, yet Table 4 lists 40.0 mm: kept as tabulated"
    if name in ("methylcyclohexane", "iso-octane"):
        nsp += "; >35 mm, near the 50 mm range limit of the lamp"
    sp_tsi(S3, S3D, "Table 4 (PDF p45)", f"pure_{name}", "pure", [name], [], "", "none", sp, tsi, M, q, nsp,
           f"MW {mw}; TSI as printed; a=4.07, b=-4.8 set by the 1-MN=100 and MCH=5 anchors (those two are assigned, not independent)")
# Appendix B (PDF p72-73): binary mixtures, mole fractions, MW, SP, TSI measured.
APPB = {
    ("iso-cetane", "1,2,4-trimethylbenzene"): [(0.11, 0.89, 131.48, 8.1, 61), (0.24, 0.76, 145.76, 10.2, 53),
                                               (0.42, 0.58, 164.36, 11.9, 52), (0.53, 0.47, 176.39, 14.2, 46),
                                               (0.66, 0.34, 190.00, 18.5, 37), (0.81, 0.19, 206.30, 24.5, 30)],
    ("iso-octane", "toluene"): [(0.10, 0.90, 94.39, 8.4, 41), (0.24, 0.76, 97.35, 10.2, 34), (0.39, 0.61, 100.79, 12.3, 29),
                                (0.49, 0.51, 102.99, 14.5, 24), (0.66, 0.34, 106.69, 20.0, 17), (0.72, 0.28, 108.05, 21.7, 15),
                                (0.78, 0.22, 109.48, 24.6, 13), (0.85, 0.15, 110.98, 29.4, 11), (0.92, 0.08, 112.56, 35.0, 8.3)],
    ("m-xylene", "n-dodecane"): [(0.17, 0.83, 159.35, 36.9, 13), (0.33, 0.67, 149.43, 26.8, 18), (0.34, 0.66, 148.43, 26.4, 18),
                                 (0.36, 0.64, 147.49, 23.9, 20), (0.37, 0.63, 146.52, 24.0, 20), (0.39, 0.61, 145.61, 23.0, 21),
                                 (0.50, 0.50, 138.31, 19.5, 24), (0.65, 0.35, 128.57, 14.1, 32), (0.78, 0.22, 120.58, 11.9, 36),
                                 (0.88, 0.12, 113.68, 9.8, 43)],
    ("1,2,4-trimethylbenzene", "1-methylnaphthalene"): [(0.75, 0.25, 125.64, 6.9, 69), (0.60, 0.40, 128.92, 6.2, 80),
                                                        (0.40, 0.60, 133.32, 6.2, 83), (0.25, 0.75, 136.64, 6.2, 85)],
    ("decalin", "1-methylnaphthalene"): [(0.95, 0.05, 138.47, 20.0, 23), (0.89, 0.11, 138.69, 17.4, 28),
                                         (0.78, 0.22, 139.12, 13.3, 38), (0.68, 0.32, 139.53, 11.3, 46),
                                         (0.52, 0.48, 140.15, 9.2, 57), (0.32, 0.68, 140.94, 7.2, 75), (0.13, 0.87, 141.67, 6.1, 89)],
    ("methylcyclohexane", "1-methylnaphthalene"): [(0.98, 0.02, 99.20, 29.9, 8.7), (0.95, 0.05, 100.24, 23.7, 12),
                                                   (0.91, 0.09, 102.31, 18.9, 17), (0.81, 0.19, 106.39, 14.7, 25),
                                                   (0.72, 0.28, 110.54, 11.5, 34), (0.62, 0.38, 114.82, 10.5, 40),
                                                   (0.42, 0.58, 123.83, 7.3, 64), (0.21, 0.79, 132.85, 6.3, 81)],
}
for (a, b), lst in APPB.items():
    for i, (xa, xb, mw, sp, tsi) in enumerate(lst, 1):
        sp_tsi(S3, S3D, "Appendix B (PDF p72-73)", f"Mensch_bin_{a}__{b}_{i}", "binary", [a, b], [xa, xb], "fraction", "mole",
               sp, tsi, M, "", f"MW {mw}; fraction order = column order of the appendix block ({a}, {b}), checked against printed MW")
# Appendix A (PDF p70-71): repeated SP of iso-octane/toluene standards (liquid-volume %). date SP [T_amb_C] [P_inHg] MW TSI
APPA = {
    ("Std1", 60, 40): """11/27/2007 14.3 103.0 24
12/12/2007 15.3 103.0 23
1/3/2008 14.3 103.0 24
1/25/2008 15.0 103.0 23
1/28/2008 14.7 103.0 24
2/6/2008 14.3 103.0 24
3/7/2008 14.7 19.5 28.86 103.0 24
3/18/2008 14.0 20.0 28.56 103.0 25
3/25/2008 15.0 18.0 29.21 103.0 23
5/19/2008 14.7 20.3 28.49 103.0 24
6/18/2008 14.0 20.1 28.73 103.0 25
7/22/2008 14.0 23.3 28.88 103.0 25
8/25/2008 14.3 20.1 28.91 103.0 24
9/17/2008 14.7 24.0 29.15 103.0 24
9/23/2008 14.3 24.4 29.37 103.0 24
10/7/2008 15.0 21.6 29.25 103.0 23
10/16/2008 14.7 103.0 24
10/28/2008 15.0 21.5 28.68 103.0 23
11/26/2008 14.5 21.5 28.85 103.0 24""",
    ("Std2", 75, 25): """11/27/2007 19.7 106.7 17
12/12/2007 22.0 106.7 15
12/18/2007 20.0 106.7 17
1/3/2008 19.0 106.7 18
1/14/2008 20.0 106.7 17
1/25/2008 20.0 106.7 17
2/6/2008 20.0 106.7 17
2/26/2008 20.0 106.7 17
3/18/2008 20.0 20.0 28.56 106.7 17
4/30/2008 20.0 21.0 28.98 106.7 17
5/9/2008 20.0 19.6 28.68 106.7 17
5/29/2008 20.0 20.6 28.99 106.7 17
7/22/2008 19.3 23.3 28.88 106.7 18
9/17/2008 19.3 24.0 29.15 106.7 18
10/16/2008 19.3 106.7 18
11/26/2008 19.8 21.5 28.85 106.7 17""",
    ("Std3", 80, 20): "1/14/2008 21.7 108.1 16",
    ("Std4", 85, 15): """2/26/2008 25.0 109.5 13
3/25/2008 24.0 19 29.21 109.5 14
4/30/2008 24.3 21 28.98 109.5 14
5/9/2008 24.3 19.6 28.68 109.5 14
5/29/2008 24.7 20.6 28.99 109.5 13
9/3/2008 25.3 24.7 109.5 13""",
    ("Std5", 90, 10): """8/25/2008 29.3 20.1 28.91 111.0 11
9/3/2008 29.0 24.7 111.0 11
9/23/2008 30.0 24.4 29.37 111.0 10""",
    ("Std6", 95, 5): """4/22/2008 32.7 23.3 29.20 112.6 9.2
6/18/2008 31.7 20.1 28.73 112.6 9.7
6/23/2008 32.0 19.7 28.91 112.6 9.5
10/7/2008 30.0 21.6 29.25 112.6 10.5""",
}
for (std, vi, vt), block in APPA.items():
    for ln in block.splitlines():
        tk = ln.split()
        date, sp, mw, tsi = tk[0], float(tk[1]), float(tk[-2]), float(tk[-1])
        mid = tk[2:-2]
        extra = f"date {date}; MW {mw}"
        if len(mid) == 2:
            extra += f"; ambient {mid[0]} C, {mid[1]} inHg"
        elif len(mid) == 1:
            extra += f"; ambient {mid[0]} C (pressure not printed)"
        sp_tsi(S3, S3D, "Appendix A (PDF p70-71)", f"Mensch_{std}_iC8_tol_{vi}_{vt}v", "binary", ["iso-octane", "toluene"],
               [vi, vt], "percent", "volume", sp, tsi, M, "replicate", extra + "; the thesis gives no mean for this standard")
# Table 7 (PDF p62): MURI JP-8 surrogates, MOLE fractions (Table 3 gives volume recipes; not used).
for sid, comps, fr, mw, sp, tsi in [
        ("Mensch_MURI_3comp_A", ["iso-cetane", "n-dodecane", "1,3,5-trimethylbenzene"], [52, 23, 25], 186.8, 22.5, 29.0),
        ("Mensch_MURI_3comp_B", ["iso-cetane", "n-dodecane", "1-methylnaphthalene"], [64, 23, 13], 202.5, 24.2, 29.3),
        ("Mensch_MURI_4comp_A", ["iso-octane", "n-dodecane", "1,3,5-trimethylbenzene", "1-methylnaphthalene"],
         [57.6, 22.8, 14.7, 5], 193.9, 24.0, 28.1)]:
    sp_tsi(S3, S3D, "Table 7 (PDF p62)", sid, "surrogate", comps, fr, "percent", "mole", sp, tsi, M, "",
           f"MW {mw}; printed mole % sum to {sum(fr):g}")
# Table 9 (PDF p64): JP-8 surrogate (22), mole %.
sp_tsi(S3, S3D, "Table 9 (PDF p64)", "Mensch_JP8_surrogate_22", "surrogate",
       ["iso-octane", "n-dodecane", "1,3,5-trimethylbenzene"], [54, 18, 28], "percent", "mole", 17.8, 24.0, M, "", "MW 125.9")
# Table 8 (PDF p63): drum of JP-8, three occasions (SP, TSI) + reported average.
for i, (sp, tsi) in enumerate([(22.8, 22.6), (22.2, 23.3), (23.8, 21.4)], 1):
    sp_tsi(S3, S3D, "Table 8 (PDF p63)", "Mensch_JP8_drum", "real_fuel", [], [], "", "none", sp, tsi, M, "replicate",
           f"occasion {i}; MW 153.3 assumed from average formula C11H21; reported mean is SP 23.0 mm, TSI 22")
sp_tsi(S3, S3D, "Table 8 (PDF p63)", "Mensch_JP8_drum", "real_fuel", [], [], "", "none", 23.0, 22, M, "reported_mean",
       "mean of the three occasions as printed")

# --------------------------------------------------------------------------------------------- S4 Rainville/Lemaire 2025
S4, S4D = "S4_rainville2025", "10.1021/acs.energyfuels.5c00385"
# Table 5 (PDF p7): iso-octane/toluene calibration fuels; reference standard SP, manual, FURTI, VBA (mm, mean of 3 tests, +/-).
T5 = [(60, 40, 14.7, (16.6, 1.90), (14.34, 0.36), (14.97, 0.27)), (75, 25, 20.2, (19.1, 1.10), (19.45, 0.75), (20.65, 0.45)),
      (80, 20, 22.7, (20.2, 2.50), (22.04, 0.66), (22.37, 0.33)), (85, 15, 25.8, (23.5, 2.30), (24.93, 0.87), (26.14, 0.34)),
      (90, 10, 30.2, (27.4, 2.80), (29.23, 0.97), (30.87, 0.67)), (95, 5, 35.4, (38.3, 2.90), (34.39, 1.01), (36.11, 0.71)),
      (100, 0, 42.8, (41.4, 1.40), (42.21, 0.59), (42.46, 0.34))]
for vi, vt, ref, man, fur, vba in T5:
    comps, fr = (["iso-octane", "toluene"], [vi, vt]) if vt else (["iso-octane"], [])
    kind = "binary" if vt else "pure"
    sid = f"Rainville_iC8_tol_{vi}_{vt}v"
    for meth, (v, u) in [("manual ASTM D1322", man), ("FURTI", fur), ("VBA", vba)]:
        add(S4, S4D, "Table 5 (PDF p7)", sid, kind, comps, fr, "percent", "volume", "smoke_point_mm", v, "mm", meth, "",
            f"+/- {u}; mean of 3 tests")
    add(S4, S4D, "Table 5 (PDF p7)", sid, kind, comps, fr, "percent", "volume", "smoke_point_mm", ref, "mm",
        "ASTM D1322 standard value", "reference", "standard value quoted by the paper for the calibration fuel, not measured there")
# Table 1 (PDF p2) recipes; Table 6 (PDF p9) SP measured by VBA (the method the paper selects in section 2.3.3 for later work).
S = {
    "Eddings2005_S1": ("S1", "mole", ["n-octane", "n-dodecane", "n-hexadecane", "xylenes", "tetralin", "decalin"],
                       [3.5, 40.0, 5.0, 8.5, 8.0, 35.0], 28.1, 0.34),
    "Eddings2005_S2": ("S2", "mole", ["n-octane", "n-dodecane", "n-hexadecane", "xylenes", "tetralin", "decalin"],
                       [3.0, 30.0, 12.0, 15.0, 13.0, 27.0], 24.8, 0.28),
    "AksitMoss_S3": ("S3", "mass", ["n-decane", "1,3,5-trimethylbenzene"], [70, 30], 21.6, 0.28),
    "Honnet2009_Aachen": ("S4", "mass", ["n-decane", "1,2,4-trimethylbenzene"], [80, 20], 32.9, 0.36),
    "Dooley2012_2ndGen": ("S5", "mole", ["n-dodecane", "iso-octane", "1,3,5-trimethylbenzene", "n-propylbenzene"],
                          [40.41, 29.48, 7.28, 22.83], 25.5, 0.39),
    "KimVioli2014_UM2": ("S6", "mole", ["n-dodecane", "decalin", "iso-cetane", "toluene"], [28.97, 31.88, 14.24, 24.91], 24.8, 0.19),
    "Yu_S7": ("S7", "mole", ["n-dodecane", "decalin", "iso-cetane", "iso-octane", "toluene"], [40.11, 0.5, 12.49, 9.8, 37.1], 25.4, 0.21),
    "Rainville2022_S8": ("S8", "mole", ["n-dodecane", "1,3,5-trimethylbenzene", "iso-cetane", "n-propylbenzene"],
                         [40.41, 7.28, 29.48, 22.83], 24.1, 0.24),
}
for sid, (lab, basis, comps, fr, sp, u) in S.items():
    note = f"paper label {lab}; +/- {u}"
    if sid in ("Dooley2012_2ndGen", "Honnet2009_Aachen", "KimVioli2014_UM2"):
        note += f"; same recipe as {sid} in data/raw/surrogates/fuel_surrogates.csv"
    if "xylenes" in comps:
        note += "; xylenes are an isomer mix in the source, m-xylene used as representative structure"
    add(S4, S4D, "Table 6 (PDF p9) + Table 1 (PDF p2)", sid, "surrogate", comps, fr, "percent", basis, "smoke_point_mm", sp,
        "mm", "VBA", "", note)
add(S4, S4D, "Table 6 (PDF p9)", "JetA_Rainville2025", "real_fuel", [], [], "", "none", "smoke_point_mm", 22.8, "mm", "VBA", "",
    "+/- 0.29; Jet-A, no recipe (approximate formula C11.6H21.7)")

# --------------------------------------------------------------------------------------------- S5 Botero preprint 153
S5, S5D = "S5_botero_preprint153", ""
# Table 1 (PDF p8): n-heptane/toluene, volume %, SP (mm) and TSI.
for sid, comps, fr, sp, spu, tsi, tsiu, q in [
        ("Botero_heptane", ["n-heptane"], [], 73.0, 18, -2.6, None, "literature_value"),
        ("Botero_H90T10", ["n-heptane", "toluene"], [90, 10], 39.0, 2.8, 2.1, 3.30, ""),
        ("Botero_H70T30", ["n-heptane", "toluene"], [70, 30], 20.5, 1.7, 11.0, 4.81, ""),
        ("Botero_H50T50", ["n-heptane", "toluene"], [50, 50], 12.0, 1.1, 19.9, 6.53, ""),
        ("Botero_H30T70", ["n-heptane", "toluene"], [30, 70], 9.8, 0.8, 30.5, 8.55, ""),
        ("Botero_toluene", ["toluene"], [], 7.8, 0.5, 39.4, 10.1, "")]:
    kind = "pure" if len(comps) == 1 else "binary"
    nsp = f"+/- {spu}" + ("; heptane value is quoted from the preprint's ref [65] (table footnote), not measured there" if q else "")
    ntsi = (f"+/- {tsiu}" if tsiu else "no uncertainty printed") + ("; derived from the literature SP" if q else "")
    add(S5, S5D, "Table 1 (PDF p8)", sid, kind, comps, fr, "percent", "volume", "smoke_point_mm", sp, "mm",
        "wick-fed smoke point lamp", q, nsp)
    add(S5, S5D, "Table 1 (PDF p8)", sid, kind, comps, fr, "percent", "volume", "tsi", tsi, "-",
        "wick-fed smoke point lamp", q, ntsi)
add(S5, S5D, "Table 1 (PDF p8)", "Botero_gasoline", "real_fuel", [], [], "", "none", "smoke_point_mm", 11.8, "mm",
    "wick-fed smoke point lamp", "", "+/- 0.4; commercial gasoline, class composition only (no recipe); no TSI printed")

# --------------------------------------------------------------------------------------------- S6 Boehm 2022
S6, S6D = "S6_boehm2022", "10.1021/acs.energyfuels.1c03794"
LAMP = "wick-fed smoke point lamp"
add(S6, S6D, "Table 5 (PDF p10)", "Boehm_n-butylcyclohexane", "pure", ["n-butylcyclohexane"], [], "", "none",
    "smoke_point_mm", 45.1, "mm", LAMP, "", "+/- 1.3 (95 % repeatability)")
add(S6, S6D, "Table 5 (PDF p10)", "Boehm_n-butylcyclohexane", "pure", ["n-butylcyclohexane"], [], "", "none", "tsi", 8.2, "-", LAMP)
add(S6, S6D, "Table 5 (PDF p10)", "Boehm_1,4-dimethylcyclooctane", "pure", ["1,4-dimethylcyclooctane"], [], "", "none",
    "smoke_point_mm", 36.8, "mm", LAMP, "", "+/- 1.2; not in the paper's training set")
add(S6, S6D, "Table 5 (PDF p10)", "Boehm_1,4-dimethylcyclooctane", "pure", ["1,4-dimethylcyclooctane"], [], "", "none",
    "tsi", 10.8, "-", LAMP)
sc = ["n-hexylbenzene", "m-xylene", "n-octane", "iso-cetane", "n-undecane", "decalin"]
sf = [0.076, 0.117, 0.118, 0.166, 0.273, 0.249]
add(S6, S6D, "Table 5 (PDF p10)", "Boehm_surrogate1", "surrogate", sc, sf, "fraction", "mole", "smoke_point_mm", 26.0, "mm", LAMP,
    "", "+/- 1.0; decalin is cis-decalin in the paper (stereo dropped); printed mole fractions sum to 0.999")
add(S6, S6D, "Table 5 (PDF p10)", "Boehm_surrogate1", "surrogate", sc, sf, "fraction", "mole", "tsi", 18.5, "-", LAMP)
# Table 2 (PDF p5, image-only table, read from a 200 dpi render): 'This work' column.
for lab, name, comps, fr, sp, tsi in [
        ("A", "1-methylnaphthalene", ["1-methylnaphthalene"], [], 5.12, 100),
        ("B", "methylcyclohexane", ["methylcyclohexane"], [], 43.56, 5),
        ("C", "toluene", ["toluene"], [], 7.59, 41.8),
        ("E", "1,3,5-trimethylbenzene", ["1,3,5-trimethylbenzene"], [], 7.19, 58.9),
        ("F", "toluene+iso-octane_40_60v", ["toluene", "iso-octane"], [40, 60], 14.88, 22.4)]:
    kind = "binary" if fr else "pure"
    sid = f"Boehm_T2_{name}"
    anchor = lab in ("A", "B")
    add(S6, S6D, "Table 2 (PDF p5)", sid, kind, comps, fr, "percent", "volume", "smoke_point_mm", sp, "mm", LAMP, "",
        f"fuel {lab}, 'This work' column")
    add(S6, S6D, "Table 2 (PDF p5)", sid, kind, comps, fr, "percent", "volume", "tsi", tsi, "-", LAMP,
        "anchor_assigned" if anchor else "",
        f"fuel {lab}" + ("; assigned scale anchor (1-MN=100, MCH=5), not independent" if anchor else ""))

# --------------------------------------------------------------------------------------------- S7 YSI unified (Dataverse)
S7, S7D = "S7_ysi_unified", "10.7910/DVN/7HGFT8"
x = pd.read_excel(OUT / "ysi_unified" / "ysi_database_volume2.xlsx")
x.columns = ["species", "formula", "cas", "ysi", "ysi_err", "smiles"]
ysi_rows = []
for r in x.itertuples():
    if not isinstance(r.smiles, str):
        continue
    m = Chem.MolFromSmiles(r.smiles)
    if m is None or any(a.GetSymbol() != "C" for a in m.GetAtoms()):
        continue
    ysi_rows.append(dict(source_id=S7, doi=S7D, table_or_page="YSI Database Volume 2.xlsx (sheet Database)",
                         sample_id=f"ysi_{r.cas}", kind="pure", components=json.dumps([r.species]),
                         smiles=json.dumps([r.smiles]), fractions="[]", fraction_unit="", fraction_basis="none",
                         property="ysi", value=float(r.ysi), unit="YSI", T_C="",
                         method="1000 ppm doped into methane/N2 coflow flame, unified scale", qualifier="",
                         note=f"CAS {r.cas}; unified-scale error {r.ysi_err}"))

for fn, src in [("S3_mensch2009_smoke_tsi.csv", S3), ("S4_rainville2025_smoke.csv", S4),
                ("S5_botero_smoke_tsi.csv", S5), ("S6_boehm2022_smoke_tsi.csv", S6)]:
    df = pd.DataFrame([r for r in rows if r["source_id"] == src], columns=COLS)
    df.to_csv(OUT / fn, index=False)
    print(fn, len(df))
pd.DataFrame(ysi_rows, columns=COLS).to_csv(OUT / "S7_ysi_unified_hc.csv", index=False)
print("S7_ysi_unified_hc.csv", len(ysi_rows), "of", len(x))
