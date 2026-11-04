#!/usr/bin/env python
"""Extract Fortin et al. 2012 (RP-1/RP-2 orthogonality, NIST copy) Tables 1-4.

Output: data/raw/rp1/fortin2012_tables.csv (canonical long rows)
  fuel_id, fuel_class, property, value_std, unit_std, vol_pct, source_pdf, page

Properties: initial_boiling_sustained_C, initial_boiling_vapor_rise_C,
            distillation_T<fac>_C  (fac = 0.025 | 5 | 10 | ... 90-95 whatever printed)
Values are Sydney-Young-corrected to standard atmospheric pressure (as printed).
"""
import re
from pathlib import Path

import pdfplumber
import pandas as pd

PDF = Path("data/raw/rp1/fortin_2012_rp1_blends.pdf")
OUT = Path("data/raw/rp1/fortin2012_tables.csv")

NUM = re.compile(r"^-?\d+(?:\.\d+)?$")
BLENDS_RP1 = ["Blend 1", "Blend 2", "Blend 3", "Blend 4", "Blend 5", "Blend 6",
              "Blend 7", "Blend 8", "Blend 9", "Blend 10", "Blend 11", "RP-1 4572"]
BLENDS_RP2 = ["RP-2 EAFB", "RP-2 5433"]


def _cells(line_text):
    return [c for c in re.split(r"\s+", line_text.strip()) if c]


def parse_boiling(pdf, pno, blends, colname_first):
    page = pdf.pages[pno]
    rows = []
    for ln in page.extract_text_lines():
        t = ln["text"].replace("\u2103", "C").strip()
        m = re.match(rf"^(?:{colname_first})?\s*({r'|'.join(re.escape(b) for b in blends)})", t, re.I)
        # simpler: line starts with the blend id after strip
        head = t.split()[0] if t.split() else ""
        hits = [b for b in blends if t.lower().startswith(b.lower())]
        if not hits:
            continue
        blend = hits[0]
        rest = t[len(blend):]
        # pair numbers in (pressure) style
        nums = re.findall(r"-?\d+(?:\.\d+)?", rest)
        if len(nums) < 3:
            continue
        p_kpa = float(nums[0])
        T_sb = float(nums[1])
        T_vr = float(nums[2])
        rows.append((blend, p_kpa, T_sb, T_vr, pno))
    return rows


def parse_distillation(pdf, pno, blends, expected_cells):
    page = pdf.pages[pno]
    rows = []
    started = set()
    for ln in page.extract_text_lines():
        t = ln["text"].replace("\u2103", "C").strip()
        cells = _cells(t)
        if not cells or not NUM.match(cells[0]):
            continue
        vf = float(cells[0])
        vals = [c for c in cells[1:] if NUM.match(c)]
        # value columns may carry per-blend pressure trailing numbers (e.g. '83.0', 'kPa')
        if len(vals) < expected_cells:
            continue
        # last expected_cells numbers correspond to the blends
        vals = vals[len(vals) - expected_cells:] if len(vals) > expected_cells else vals
        for b, v in zip(blends, vals):
            rows.append((b, vf, float(v), pno))
        started.add(vf)
    return rows


def main():
    pdf = pdfplumber.open(PDF)
    recs = []
    for blend, p_kpa, T_sb, T_vr, pno in parse_boiling(pdf, 18, BLENDS_RP1, "Blend"):
        recs += [
            dict(fuel_id=blend, fuel_class="RP-1", property="initial_boiling_sustained_C",
                 value_std=T_sb, unit_std="C", vol_pct="", source_pdf=PDF.name, page=pno + 1,
                 T_C="", P_kPa=101.325, note=f"sydney-young corrected; measured at {p_kpa} kPa"),
            dict(fuel_id=blend, fuel_class="RP-1", property="initial_boiling_vapor_rise_C",
                 value_std=T_vr, unit_std="C", vol_pct="", source_pdf=PDF.name, page=pno + 1,
                 T_C="", P_kPa=101.325, note=f"sydney-young corrected; measured at {p_kpa} kPa"),
        ]
    for blend, p_kpa, T_sb, T_vr, pno in parse_boiling(pdf, 18, BLENDS_RP2, "Blend"):
        recs += [
            dict(fuel_id=blend, fuel_class="RP-2", property="initial_boiling_sustained_C",
                 value_std=T_sb, unit_std="C", vol_pct="", source_pdf=PDF.name, page=pno + 1,
                 T_C="", P_kPa=101.325, note=f"sydney-young corrected; measured at {p_kpa} kPa"),
            dict(fuel_id=blend, fuel_class="RP-2", property="initial_boiling_vapor_rise_C",
                 value_std=T_vr, unit_std="C", vol_pct="", source_pdf=PDF.name, page=pno + 1,
                 T_C="", P_kPa=101.325, note=f"sydney-young corrected; measured at {p_kpa} kPa"),
        ]

    # Table 3 (page 19 -> index 18): 12 RP-1 blends
    d3 = parse_distillation(pdf, 19, BLENDS_RP1, 12)
    for blend, vf, T, pno in d3:
        recs.append(dict(fuel_id=blend, fuel_class="RP-1",
                         property=f"distillation_T{g(vf)}_C", value_std=T, unit_std="C",
                         vol_pct=vf, source_pdf=PDF.name, page=pno + 1, T_C="", P_kPa=101.325,
                         note="sydney-young corrected"))
    # Table 4 (page 20 -> index 19): 2 RP-2 blends
    d4 = parse_distillation(pdf, 20, BLENDS_RP2, 2)
    for blend, vf, T, pno in d4:
        recs.append(dict(fuel_id=blend, fuel_class="RP-2",
                         property=f"distillation_T{g(vf)}_C", value_std=T, unit_std="C",
                         vol_pct=vf, source_pdf=PDF.name, page=pno + 1, T_C="", P_kPa=101.325,
                         note="sydney-young corrected"))
    pdf.close()

    df = pd.DataFrame(recs)
    if df.empty:
        raise SystemExit("nothing parsed!")
    # sanity: T range and per-blend counts
    print(df.groupby(["fuel_class", "property"]).size().to_string())
    print("T range:", df.value_std.min(), "-", df.value_std.max(), "C")
    assert df.value_std.between(150, 350).all(), "distillation/bp values outside 150-350 C"
    counts = df[df.property.str.startswith("distillation")].groupby("fuel_id").size()
    print(counts.to_string())
    fv = df[df.property.str.startswith("distillation")].groupby("fuel_id").vol_pct.nunique()
    assert fv.nunique() <= 1, "inconsistent volume-fraction grid across blends"
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}: {len(df)} rows")


def g(v):
    v = float(v)
    return f"{v:g}"


if __name__ == "__main__":
    main()
