#!/usr/bin/env python3
"""parse_kaust_jameel.py - Parse octane data from Abdul Jameel et al. (2018).

Source: Abdul Gani Abdul Jameel, Vincent Van Oudenhoven, Abdul-Hamid Emwas, S. Mani Sarathy.
"Predicting Octane Number Using Nuclear Magnetic Resonance Spectroscopy and Artificial Neural Networks"
Energy & Fuels 2018, 32, 6309-6329. DOI: 10.1021/acs.energyfuels.8b00556.
KAUST Repository: https://repository.kaust.edu.sa/handle/10754/627695

Outputs:
  - data/interim/external_staging/kaust_jameel_table3_pure.csv (128 pure hydrocarbons)
  - data/interim/external_staging/kaust_jameel_table4_blends.csv (123 multi-component blends)
  - data/interim/external_staging/kaust_jameel_table5_face.csv (30 FACE gasoline + ethanol blends)

NOTE: This script ONLY outputs staged files in data/interim/external_staging/.
It does NOT modify active canonical datasets or the production pipeline.
"""

from pathlib import Path
import re
import sys
import pandas as pd
import pdfplumber

RAW_PDF = Path("data/raw/kaust_jameel2018/acs.energyfuels.8b00556.pdf")
STAGING_DIR = Path("data/interim/external_staging")


def parse_table3_pure(pdf: pdfplumber.PDF) -> pd.DataFrame:
    """Parse Table 3: RON and MON of pure hydrocarbons (Pages 25-27, index 24-26)."""
    rows = []
    for p_idx in [24, 25, 26]:
        txt = pdf.pages[p_idx].extract_text() or ""
        for line in txt.splitlines():
            line = line.strip()
            # Match optional manuscript line number, then s_no, then name, then ron, mon (optional letters like 108a)
            m = re.match(r"^(?:\d+\s+)?(\d{1,3})\s+([a-zA-Z0-9\-\,\(\)\[\]\s\']+?)\s+([\-\d\.]+[a-zA-Z]?)\s+([\-\d\.]+[a-zA-Z]?)$", line)
            if m:
                s_no_str, name, ron_str, mon_str = m.groups()
                try:
                    s_no = int(s_no_str)
                    # strip any footnote characters
                    ron = float(re.sub(r"[^\d\.\-]", "", ron_str))
                    mon = float(re.sub(r"[^\d\.\-]", "", mon_str))
                    if 1 <= s_no <= 128:
                        rows.append({
                            "s_no": s_no,
                            "name": name.strip(),
                            "ron": ron,
                            "mon": mon,
                            "source_ref": "Flynn & Wall [70]",
                            "provenance": "Abdul Jameel (2018) Table 3",
                        })
                except ValueError:
                    continue

    df = pd.DataFrame(rows).drop_duplicates(subset=["s_no"]).sort_values("s_no").reset_index(drop=True)
    return df


def parse_table4_blends(pdf: pdfplumber.PDF) -> pd.DataFrame:
    """Parse Table 4: RON and MON of hydrocarbon blends (Pages 28-30, index 27-29)."""
    rows = []
    # In Table 4, rows 1-123 are listed across pages 28-30
    for p_idx in [27, 28, 29]:
        txt = pdf.pages[p_idx].extract_text() or ""
        for line in txt.splitlines():
            line = line.strip()
            # Match S.no and blend name, numbers, ron, mon, [ref]
            # Pattern: optional line num, s_no (1-123), Name, series of numbers, ron, mon, [ref]
            m = re.match(r"^(?:\d+\s+)?(\d{1,3})\s+([A-Za-z0-9\s\+\-]+?)\s+([\d\.\s]+?)\s+\[(\d+)\]$", line)
            if m:
                s_no_str, blend_name, nums_str, ref = m.groups()
                if "TPRF 3+20" in blend_name:
                    s_no_str = "22"
                try:
                    s_no = int(s_no_str)
                    if not (1 <= s_no <= 123):
                        continue
                    tokens = [float(x) for x in nums_str.split()]
                    if len(tokens) < 3:
                        continue
                    ron = tokens[-2]
                    mon = tokens[-1]
                    comp_vals = tokens[:-2]

                    # Assign components based on blend class & number of values
                    comp_dict = {}
                    name_clean = blend_name.strip()
                    if 1 <= s_no <= 9:  # PRF+ethanol: n-heptane, i-octane, ethanol
                        comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "ethanol": comp_vals[2]}
                    elif 10 <= s_no <= 39:  # TPRF 1-4: n-heptane, i-octane, toluene, ethanol
                        if len(comp_vals) == 4:
                            comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2], "ethanol": comp_vals[3]}
                        elif len(comp_vals) == 3:
                            comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2]}
                    elif 40 <= s_no <= 44:  # TMB+ethanol: 1,2,4-trimethylbenzene, ethanol
                        comp_dict = {"1,2,4-trimethylbenzene": comp_vals[0], "ethanol": comp_vals[1]}
                    elif 45 <= s_no <= 49:  # CP+ethanol: cyclopentane, ethanol
                        comp_dict = {"cyclopentane": comp_vals[0], "ethanol": comp_vals[1]}
                    elif 50 <= s_no <= 54:  # Hex+ethanol: 1-hexene, ethanol
                        comp_dict = {"1-hexene": comp_vals[0], "ethanol": comp_vals[1]}
                    elif 55 <= s_no <= 97:  # TPRF 1-43: n-heptane, i-octane, toluene
                        if len(comp_vals) == 3:
                            comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2]}
                    elif 98 <= s_no <= 123:  # MC 1-26: n-heptane, i-octane, toluene, cyclopentane
                        if len(comp_vals) == 4:
                            comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2], "cyclopentane": comp_vals[3]}

                    rows.append({
                        "s_no": s_no,
                        "blend_name": name_clean,
                        "composition_vol_pct": comp_dict,
                        "ron": ron,
                        "mon": mon,
                        "source_ref": f"[{ref}]",
                        "provenance": "Abdul Jameel (2018) Table 4",
                    })
                except Exception:
                    continue

    df = pd.DataFrame(rows).drop_duplicates(subset=["s_no"]).sort_values("s_no").reset_index(drop=True)
    return df


def parse_table5_face(pdf: pdfplumber.PDF) -> pd.DataFrame:
    """Parse Table 5: RON and MON of gasoline-ethanol blends (Page 31, index 30)."""
    rows = []
    txt = pdf.pages[30].extract_text() or ""
    for line in txt.splitlines():
        line = line.strip()
        # Pattern: optional line num, s_no (1-30), name e.g. FACE A+10, base_vol, ethanol_vol, ron, mon
        m = re.match(r"^(?:\d+\s+)?(\d{1,2})\s+(FACE\s+[A-Z]\+\d+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)$", line)
        if m:
            s_no, name, base_vol, eth_vol, ron, mon = m.groups()
            base_fuel = name.split("+")[0].strip()
            rows.append({
                "s_no": int(s_no),
                "blend_name": name.strip(),
                "base_fuel": base_fuel,
                "base_vol_pct": float(base_vol),
                "ethanol_vol_pct": float(eth_vol),
                "ron": float(ron),
                "mon": float(mon),
                "source_ref": "Foong et al. [26]",
                "provenance": "Abdul Jameel (2018) Table 5",
            })

    df = pd.DataFrame(rows).drop_duplicates(subset=["s_no"]).sort_values("s_no").reset_index(drop=True)
    return df


def main():
    if not RAW_PDF.exists():
        print(f"ERROR: Raw PDF not found at {RAW_PDF}")
        sys.exit(1)

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Opening {RAW_PDF}...")

    with pdfplumber.open(RAW_PDF) as pdf:
        # Table 3
        df_t3 = parse_table3_pure(pdf)
        out_t3 = STAGING_DIR / "kaust_jameel_table3_pure.csv"
        df_t3.to_csv(out_t3, index=False)
        print(f"Saved Table 3 (Pure): {len(df_t3)} rows -> {out_t3}")

        # Table 4
        df_t4 = parse_table4_blends(pdf)
        out_t4 = STAGING_DIR / "kaust_jameel_table4_blends.csv"
        df_t4.to_csv(out_t4, index=False)
        print(f"Saved Table 4 (Blends): {len(df_t4)} rows -> {out_t4}")

        # Table 5
        df_t5 = parse_table5_face(pdf)
        out_t5 = STAGING_DIR / "kaust_jameel_table5_face.csv"
        df_t5.to_csv(out_t5, index=False)
        print(f"Saved Table 5 (FACE Blends): {len(df_t5)} rows -> {out_t5}")

    print("All KAUST tables parsed and staged successfully!")


if __name__ == "__main__":
    main()
