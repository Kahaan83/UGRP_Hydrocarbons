#!/usr/bin/env python
"""Build data/processed/fortin2012_rp1_blends.csv from Fortin et al. 2012.

Sources:
  - data/raw/rp1/fortin2012_tables.csv  (distillation curves, previously extracted)
  - Hard-coded Table 1/2 initial boiling temperatures (Tables 1 & 2 of the paper)
  - Hard-coded T90 rows for RP-1 Blends 1-11 (Table 3, missing from prior extraction)
  - Hard-coded RP-2 EAFB T0.025 (Table 4, missing from prior extraction)

Paper: Fortin T.J. et al., "Assessment of Compositional Variability of RP-1 and RP-2
       with the Advanced Distillation Curve Approach", Energy & Fuels 26:6711 (2012),
       doi:10.1021/ef3004009. PDF: data/raw/rp1/fortin_2012_rp1_blends.pdf.

Output schema:
  fuel_id          : string ('Blend 1'..'Blend 11', 'RP-1 4572', 'RP-2 EAFB', 'RP-2 5433')
  fuel_class       : 'RP-1' or 'RP-2'
  property         : canonical property name
  value_std        : numeric value in unit_std
  unit_std         : 'C' (all temperatures)
  vol_pct          : distillate volume fraction (%) — blank for initial boiling rows
  P_kPa            : ambient pressure at measurement (kPa) — from each blend column header
  source_pdf       : 'fortin_2012_rp1_blends.pdf'
  doi              : '10.1021/ef3004009'
  page             : PDF page number (1-based)
  note             : measurement/correction details
  tier             : 'A' (all — directly measured, published in peer-reviewed journal)

Properties:
  distillation_T<frac>_C   : ADC temperature at volume fraction <frac> %, Sydney-Young corrected
  initial_boiling_sustained_C  : sustained boiling temperature, Sydney-Young corrected
  initial_boiling_vapor_rise_C : vapor rise temperature (= IBT), Sydney-Young corrected

Composition notes (no molecular composition available — see README.md):
  The 11 orthogonal blends were prepared to cover RP-1 spec space; exact recipes are
  proprietary. GC-MS analysis confirms primarily linear/branched paraffins, cycloparaffins,
  alkenes, and minor aromatics — same class composition as RP-1. RP-2 has fewer aromatics.
  No quantitative molecular composition is reported in this paper.

Run: python build_fortin2012_processed.py
"""
import pandas as pd
from pathlib import Path

SRC_CSV = Path("data/raw/rp1/fortin2012_tables.csv")
OUT = Path("data/processed/fortin2012_rp1_blends.csv")
DOI = "10.1021/ef3004009"
PDF = "fortin_2012_rp1_blends.pdf"

# Ambient pressures at which each blend was measured (from Table 3 column headers)
BLEND_PRESSURE_RP1 = {
    "Blend 1": 83.0, "Blend 2": 83.5, "Blend 3": 84.0, "Blend 4": 83.6,
    "Blend 5": 84.1, "Blend 6": 83.5, "Blend 7": 82.9, "Blend 8": 81.9,
    "Blend 9": 81.7, "Blend 10": 81.8, "Blend 11": 80.8, "RP-1 4572": 83.7,
}
BLEND_PRESSURE_RP2 = {
    "RP-2 EAFB": 82.7, "RP-2 5433": 83.1,
}

# RP-1 T90 data from Table 3 (page 20); RP-1 4572 column blank at T90 row
RP1_T90 = {
    "Blend 1": 239.8, "Blend 2": 247.1, "Blend 3": 246.9, "Blend 4": 235.9,
    "Blend 5": 249.6, "Blend 6": 235.1, "Blend 7": 255.4, "Blend 8": 255.7,
    "Blend 9": 250.9, "Blend 10": 258.1, "Blend 11": 254.6,
}

# Table 1 (page 19): RP-1 initial boiling, Sydney-Young corrected
# (fuel_id, measured_pressure_kPa, sustained_C, vapor_rise_C)
TABLE1 = [
    ("Blend 1",    83.0, 208.0, 211.2),
    ("Blend 2",    83.5, 207.6, 210.3),
    ("Blend 3",    84.0, 204.5, 207.3),
    ("Blend 4",    83.6, 206.5, 208.6),
    ("Blend 5",    84.1, 210.2, 212.8),
    ("Blend 6",    83.5, 205.5, 207.7),
    ("Blend 7",    82.9, 212.1, 214.0),
    ("Blend 8",    81.9, 206.1, 208.9),
    ("Blend 9",    81.7, 209.3, 211.0),
    ("Blend 10",   81.8, 211.1, 213.6),
    ("Blend 11",   80.8, 206.3, 208.9),
    ("RP-1 4572",  83.7, 201.5, 203.9),
]

# Table 2 (page 19): RP-2 initial boiling, Sydney-Young corrected
TABLE2 = [
    ("RP-2 EAFB", 82.7, 202.9, 206.4),
    ("RP-2 5433", 83.1, 214.3, 216.5),
]


def main():
    rows = []

    # 1. Load existing distillation data (already extracted by fortin_extract.py)
    df_src = pd.read_csv(SRC_CSV)
    # Drop any T90 rows — we add them fresh (existing ones are RP-2 only, which is fine)
    # Actually keep RP-2 T90 from the existing CSV; only add RP-1 T90 which is missing
    rp2_t90_existing = set(
        df_src[df_src.property == "distillation_T90_C"]["fuel_id"].tolist()
    )

    for _, r in df_src.iterrows():
        fid = str(r["fuel_id"])
        p_meas = BLEND_PRESSURE_RP1.get(fid) or BLEND_PRESSURE_RP2.get(fid)
        rows.append({
            "fuel_id": fid,
            "fuel_class": str(r["fuel_class"]),
            "property": str(r["property"]),
            "value_std": float(r["value_std"]),
            "unit_std": "C",
            "vol_pct": r["vol_pct"] if pd.notna(r.get("vol_pct")) else "",
            "P_kPa": p_meas if p_meas else "",
            "source_pdf": PDF,
            "doi": DOI,
            "page": int(r["page"]),
            "note": str(r.get("note", "sydney-young corrected")),
            "tier": "A",
        })

    # 2. Add RP-1 T90 rows (missing from fortin2012_tables.csv)
    for fid, val in RP1_T90.items():
        p_meas = BLEND_PRESSURE_RP1[fid]
        rows.append({
            "fuel_id": fid,
            "fuel_class": "RP-1",
            "property": "distillation_T90_C",
            "value_std": val,
            "unit_std": "C",
            "vol_pct": 90.0,
            "P_kPa": p_meas,
            "source_pdf": PDF,
            "doi": DOI,
            "page": 20,
            "note": "sydney-young corrected; added from Table 3 row 90%",
            "tier": "A",
        })

    # 3. Add RP-2 EAFB T0.025 = 206.1 (marked '-' for 5433; not missing from EAFB)
    if "RP-2 EAFB" not in rp2_t90_existing or True:
        # Check if EAFB T0.025 is already present
        existing_rp2_eafb_t0025 = df_src[
            (df_src.fuel_id == "RP-2 EAFB") &
            (df_src.property == "distillation_T0.025_C")
        ]
        if len(existing_rp2_eafb_t0025) == 0:
            rows.append({
                "fuel_id": "RP-2 EAFB",
                "fuel_class": "RP-2",
                "property": "distillation_T0.025_C",
                "value_std": 206.1,
                "unit_std": "C",
                "vol_pct": 0.025,
                "P_kPa": 82.7,
                "source_pdf": PDF,
                "doi": DOI,
                "page": 21,
                "note": "sydney-young corrected; added from Table 4 row 0.025%",
                "tier": "A",
            })

    # 4. Add initial boiling temperatures (Tables 1 & 2, page 19)
    for tbl, fuel_class, data in [
        (1, "RP-1", TABLE1),
        (2, "RP-2", TABLE2),
    ]:
        for fid, p_meas, sustained, vapor_rise in data:
            note_base = f"sydney-young corrected; measured at {p_meas} kPa (Table {tbl}, page 19)"
            rows.append({
                "fuel_id": fid,
                "fuel_class": fuel_class,
                "property": "initial_boiling_sustained_C",
                "value_std": sustained,
                "unit_std": "C",
                "vol_pct": "",
                "P_kPa": p_meas,
                "source_pdf": PDF,
                "doi": DOI,
                "page": 19,
                "note": note_base + "; uncertainty ~2 °C (k=2)",
                "tier": "A",
            })
            rows.append({
                "fuel_id": fid,
                "fuel_class": fuel_class,
                "property": "initial_boiling_vapor_rise_C",
                "value_std": vapor_rise,
                "unit_std": "C",
                "vol_pct": "",
                "P_kPa": p_meas,
                "source_pdf": PDF,
                "doi": DOI,
                "page": 19,
                "note": note_base + "; uncertainty ~0.2 °C (k=2); thermodynamically consistent IBT",
                "tier": "A",
            })

    df = pd.DataFrame(rows)

    # Sanity checks
    dist = df[df.property.str.startswith("distillation_")]
    assert dist.value_std.between(195, 270).all(), \
        f"distillation temps outside 195-270 C: {dist[~dist.value_std.between(195,270)][['fuel_id','property','value_std']]}"
    ib = df[df.property.str.startswith("initial_boiling")]
    assert ib.value_std.between(195, 220).all(), \
        f"initial boiling temps outside 195-220 C: {ib[~ib.value_std.between(195,220)][['fuel_id','property','value_std']]}"

    # Row counts summary
    print("=== fortin2012_rp1_blends.csv summary ===")
    print(f"Total rows: {len(df)}")
    print("\nRows by fuel_class and fuel_id:")
    print(df.groupby(["fuel_class", "fuel_id"]).size().to_string())
    print("\nRows by property:")
    print(df.groupby("property").size().sort_index().to_string())
    print("\nValue range by property:")
    for prop, g in df.groupby("property"):
        print(f"  {prop}: {g.value_std.min():.1f} – {g.value_std.max():.1f} C (n={len(g)})")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT}: {len(df)} rows")


if __name__ == "__main__":
    main()
