#!/usr/bin/env python3
"""parse_kaust_jameel.py - Comprehensive parser for KAUST Octane Benchmark Tables 3, 4, and 5.

Source: Abdul Jameel et al. Energy & Fuels 2018, DOI: 10.1021/acs.energyfuels.8b00556.
Tables:
  - Table 3: Pure hydrocarbons (128 compounds with RON and MON)
  - Table 4: Hydrocarbon blends (123 blends with exact volume fractions: PRF, TPRF 1-43, MC 1-26, TMB, CP, Hex)
  - Table 5: Real fuels / FACE gasoline + ethanol blends (30 blends)

NOTE: Standalone parser script. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import re
import json
import pandas as pd
import pdfplumber

RAW_PDF = Path("data/raw/kaust_jameel2018/acs.energyfuels.8b00556.pdf")
STAGING_DIR = Path("data/interim/external_staging")

def parse_table3_pure(pdf: pdfplumber.PDF) -> pd.DataFrame:
    rows = []
    for p_idx in [24, 25, 26]:
        txt = pdf.pages[p_idx].extract_text() or ""
        for line in txt.splitlines():
            line = line.strip()
            m = re.match(r"^(?:\d+\s+)?(\d{1,3})\s+([a-zA-Z0-9\-\,\(\)\[\]\s\']+?)\s+([\-\d\.]+[a-zA-Z]?)\s+([\-\d\.]+[a-zA-Z]?)$", line)
            if m:
                s_no_str, name, ron_str, mon_str = m.groups()
                try:
                    s_no = int(s_no_str)
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
    rows = []
    for p_idx in [27, 28, 29]:
        txt = pdf.pages[p_idx].extract_text() or ""
        for line in txt.splitlines():
            line = line.strip()
            # Match S.no and blend name, numbers, ron, mon, [ref]
            m = re.match(r"^(?:\d+\s+)?(\d{1,3})\s+([A-Za-z0-9\s\+\-]+?)\s+([\d\.\s]+?)\s+\[(\d+)\]$", line)
            if not m:
                continue

            s_no_str, blend_name, nums_str, ref = m.groups()
            blend_name = blend_name.strip()
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

                comp_dict = {}

                # 1-9: PRF+ethanol: n-heptane, i-octane, ethanol
                if 1 <= s_no <= 9:
                    comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "ethanol": comp_vals[2]}

                # 10-24: TPRF 1-3 + ethanol: n-heptane, i-octane, toluene, ethanol
                elif 10 <= s_no <= 24:
                    comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2], "ethanol": comp_vals[3]}

                # 25-39: TPRF 4 series (TPRF 4+1 to 4+15): n-heptane, i-octane, toluene, ethanol
                elif 25 <= s_no <= 39:
                    comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2], "ethanol": comp_vals[3]}

                # 40-44: TMB+ethanol: 1,2,4-trimethylbenzene, ethanol
                elif 40 <= s_no <= 44:
                    comp_dict = {"1,2,4-trimethylbenzene": comp_vals[0], "ethanol": comp_vals[1]}

                # 45-49: CP+ethanol: cyclopentane, ethanol
                elif 45 <= s_no <= 49:
                    comp_dict = {"cyclopentane": comp_vals[0], "ethanol": comp_vals[1]}

                # 50-54: Hex+ethanol: 1-hexene, ethanol
                elif 50 <= s_no <= 54:
                    comp_dict = {"1-hexene": comp_vals[0], "ethanol": comp_vals[1]}

                # 55-97: TPRF 1-43: n-heptane, i-octane, toluene
                elif 55 <= s_no <= 97:
                    if len(comp_vals) == 4:
                        comp_dict = {"n-heptane": comp_vals[1], "i-octane": comp_vals[2], "toluene": comp_vals[3]}
                        blend_name = f"TPRF {int(comp_vals[0])}"
                    elif len(comp_vals) == 3:
                        comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2]}

                # 98-123: MC 1-26: n-heptane, i-octane, toluene, diisobutylene
                elif 98 <= s_no <= 123:
                    if len(comp_vals) == 5:
                        comp_dict = {"n-heptane": comp_vals[1], "i-octane": comp_vals[2], "toluene": comp_vals[3], "diisobutylene": comp_vals[4]}
                        blend_name = f"MC {int(comp_vals[0])}"
                    elif len(comp_vals) == 4:
                        comp_dict = {"n-heptane": comp_vals[0], "i-octane": comp_vals[1], "toluene": comp_vals[2], "diisobutylene": comp_vals[3]}

                # Filter out zero-percentage components
                comp_dict = {k: v for k, v in comp_dict.items() if v > 0.001}

                rows.append({
                    "s_no": s_no,
                    "blend_name": blend_name,
                    "composition_vol_pct": json.dumps(comp_dict),
                    "n_comp": len(comp_dict),
                    "ron": ron,
                    "mon": mon,
                    "source_ref": f"[{ref}]",
                    "provenance": "Abdul Jameel (2018) Table 4"
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows).drop_duplicates(subset=["s_no"]).sort_values("s_no").reset_index(drop=True)
    return df

def parse_table5_face(pdf: pdfplumber.PDF) -> pd.DataFrame:
    rows = []
    txt = pdf.pages[30].extract_text() or ""
    for line in txt.splitlines():
        line = line.strip()
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
        print(f"ERROR: {RAW_PDF} not found.")
        return

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    with pdfplumber.open(RAW_PDF) as pdf:
        df3 = parse_table3_pure(pdf)
        df4 = parse_table4_blends(pdf)
        df5 = parse_table5_face(pdf)

    df3.to_csv(STAGING_DIR / "kaust_jameel_table3_pure.csv", index=False)
    df4.to_csv(STAGING_DIR / "kaust_jameel_table4_blends.csv", index=False)
    df5.to_csv(STAGING_DIR / "kaust_jameel_table5_face.csv", index=False)

    print(f"Table 3 (Pure): {len(df3)}/128 rows parsed")
    print(f"Table 4 (Blends): {len(df4)}/123 rows parsed (resolved compositions: {(df4['n_comp'] > 0).sum()})")
    print(f"Table 5 (FACE Blends): {len(df5)}/30 rows parsed")

if __name__ == "__main__":
    main()
