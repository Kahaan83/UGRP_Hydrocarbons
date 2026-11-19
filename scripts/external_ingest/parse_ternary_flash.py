#!/usr/bin/env python3
"""parse_ternary_flash.py - Parse multi-component alkane flash points from literature.

Sources:
1. Liaw et al. "Measurement of Lower Flash Point for Ternary Mixture, n-Nonane+n-Decane+n-Dodecane System"
   Korean J. Chem. Eng.
   File: data/raw/ternary_flash/liaw_ternary_nonane_decane_dodecane.pdf
2. Li et al. "Prediction of the Flash Point of Binary and Ternary Alkane Mixtures"
   Advances in Materials Science and Engineering 2014. DOI: 10.1155/2014/239023.
   File: data/raw/ternary_flash/li2014_ternary_alkanes.pdf

Outputs:
  - data/interim/external_staging/ternary_flash_liaw_korean.csv
  - data/interim/external_staging/ternary_flash_li2014.csv

NOTE: Staged output ONLY. DOES NOT alter active canonical datasets.
"""

from pathlib import Path
import re
import pandas as pd
import pdfplumber

RAW_DIR = Path("data/raw/ternary_flash")
STAGING_DIR = Path("data/interim/external_staging")


def parse_liaw_korean():
    pdf_path = RAW_DIR / "liaw_ternary_nonane_decane_dodecane.pdf"
    if not pdf_path.exists():
        print(f"Skipping {pdf_path}: not found")
        return pd.DataFrame()

    rows = []
    # Table 2 on page 4
    with pdfplumber.open(pdf_path) as pdf:
        p4 = pdf.pages[3]
        txt = p4.extract_text() or ""
        # Find lines with 3 mole fractions and an exp flash point
        # Format: x1 x2 x3 Exp ...
        # e.g.: 0.20 0.70 0.10 43.0 43.15 43.31
        for match in re.finditer(r"(\d\.\d{2})\s+(\d\.\d{2})\s+(\d\.\d{2})\s+(\d{2}\.\d)", txt):
            x1, x2, x3, exp_fp = match.groups()
            fx1, fx2, fx3, f_exp = float(x1), float(x2), float(x3), float(exp_fp)
            if abs((fx1 + fx2 + fx3) - 1.0) < 0.02:
                rows.append({
                    "system": "n-nonane + n-decane + n-dodecane",
                    "n_nonane_x": fx1,
                    "n_decane_x": fx2,
                    "n_dodecane_x": fx3,
                    "flash_point_C": f_exp,
                    "method": "ASTM D93 (closed cup)",
                    "source": "Liaw et al. (Korean J Chem Eng)",
                })

    df = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    return df


def parse_li2014():
    pdf_path = RAW_DIR / "li2014_ternary_alkanes.pdf"
    if not pdf_path.exists():
        print(f"Skipping {pdf_path}: not found")
        return pd.DataFrame()

    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        # Table 9 & 10 on page 4
        p4 = pdf.pages[3]
        txt = p4.extract_text() or ""
        lines = txt.splitlines()

        curr_table = None
        for line in lines:
            if "Table9:" in line or "Table 9:" in line or "nonane+decane+undecane" in line:
                curr_table = "nonane_decane_undecane"
            elif "Table10:" in line or "Table 10:" in line or "nonane+decane+dodecane" in line:
                curr_table = "nonane_decane_dodecane"

            # Match 3 fractions and exp flash point: 0.19 0.17 0.63 54.50
            m = re.match(r"^(\d\.\d{2})\s+(\d\.\d{2})\s+(\d\.\d{2})\s+(\d{2}\.\d{2})\s+[\d\.]+", line.strip())
            if m and curr_table:
                x1, x2, x3, exp_fp = [float(v) for v in m.groups()]
                if abs((x1 + x2 + x3) - 1.0) < 0.02:
                    if curr_table == "nonane_decane_undecane":
                        rows.append({
                            "system": "n-nonane + n-decane + n-undecane",
                            "n_nonane_x": x1,
                            "n_decane_x": x2,
                            "n_undecane_x": x3,
                            "n_dodecane_x": 0.0,
                            "flash_point_C": exp_fp,
                            "method": "Pensky-Martens closed cup",
                            "source": "Li et al. (AMSE 2014)",
                        })
                    elif curr_table == "nonane_decane_dodecane":
                        rows.append({
                            "system": "n-nonane + n-decane + n-dodecane",
                            "n_nonane_x": x1,
                            "n_decane_x": x2,
                            "n_undecane_x": 0.0,
                            "n_dodecane_x": x3,
                            "flash_point_C": exp_fp,
                            "method": "Pensky-Martens closed cup",
                            "source": "Li et al. (AMSE 2014)",
                        })

    df = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True)
    return df


def main():
    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    df_liaw = parse_liaw_korean()
    if not df_liaw.empty:
        out_liaw = STAGING_DIR / "ternary_flash_liaw_korean.csv"
        df_liaw.to_csv(out_liaw, index=False)
        print(f"Saved Liaw Korean: {len(df_liaw)} rows -> {out_liaw}")

    df_li = parse_li2014()
    if not df_li.empty:
        out_li = STAGING_DIR / "ternary_flash_li2014.csv"
        df_li.to_csv(out_li, index=False)
        print(f"Saved Li 2014: {len(df_li)} rows -> {out_li}")

    print("All ternary flash points parsed successfully!")


if __name__ == "__main__":
    main()
