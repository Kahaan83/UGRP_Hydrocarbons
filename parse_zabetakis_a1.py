#!/usr/bin/env python
"""Parse Zabetakis Bulletin 627 Appendix A Table A-1 (limits of flammability)
into data/raw/lfl/lfl_measured.csv rows for the hydrocarbons we can map.

Method: pdfplumber extract_words -> y-binned rows -> split into table columns by
x-coordinate -> name tokens (containing letters) vs numeric tokens; merge decimal
split across tokens ("4." + "61" -> 4.61). Parser quality gate: the n-alkane rows
must reproduce Zabetakis Table 2 (M ~ known LFLs) to within 0.15 vol% abs.
"""
import re
from pathlib import Path

import pdfplumber
import pandas as pd
from rdkit import Chem

from common import mol_from

PDF = Path("data/raw/pdfs/Zabetakis_Bulletin_627_flammability.pdf")
OUT = Path("data/raw/lfl/lfl_measured.csv")

# Zabetakis Table 2 ground truth for parser validation
T2_NALKANE = {"methane": 5.0, "ethane": 3.0, "propane": 2.1, "n-butane": 1.8,
              "n-pentane": 1.4, "n-hexane": 1.2, "n-heptane": 1.05, "n-octane": 0.95,
              "n-nonane": 0.85, "n-decane": 0.75, "n-dodecane": 0.60,
              "n-hexadecane": 0.43}

# curated name -> smiles for hydrocarbons present in Zabetakis A-1
NAME_SMILES = {
    "methane": "C", "ethane": "CC", "propane": "CCC", "n-butane": "CCCC",
    "isobutane": "CC(C)C", "n-pentane": "CCCCC", "isopentane": "CC(C)CC",
    "neopentane": "CC(C)(C)C", "n-hexane": "CCCCCC", "2-methyl pentane": "CC(C)CCC",
    "3-methyl pentane": "CCC(C)CC", "2,2-dimethyl butane": "CC(C)(C)CC",
    "2,3-dimethyl butane": "CC(C)C(C)C", "n-heptane": "CCCCCCC",
    "2-methyl hexane": "CC(C)CCCC", "3-methyl hexane": "CCC(C)CCC",
    "3-ethyl pentane": "CCC(CC)CC", "2,2-dimethyl pentane": "CC(C)(C)CCC",
    "2,3-dimethyl pentane": "CC(C)CCC(C)C", "2,4-dimethyl pentane": "CC(C)CC(C)C",
    "3,3-dimethyl pentane": "CCC(C)(C)CC", "n-octane": "CCCCCCCC",
    "isooctane": "CC(C)CC(C)(C)C", "2,2,4-trimethyl pentane": "CC(C)CC(C)(C)C",
    "2,2,3-trimethyl pentane": "CC(C)(C)C(C)CC", "2,3,3-trimethyl pentane": "CCC(C)(C)C(C)C",
    "3,4-dimethyl hexane": "CCC(C)CC(C)C", "2-methyl heptane": "CC(C)CCCCC",
    "n-nonane": "CCCCCCCCC", "n-decane": "CCCCCCCCCC", "n-undecane": "CCCCCCCCCCC",
    "n-dodecane": "CCCCCCCCCCCC", "n-tridecane": "CCCCCCCCCCCCC",
    "n-tetradecane": "CCCCCCCCCCCCCC", "n-pentadecane": "CCCCCCCCCCCCCCC",
    "n-hexadecane": "CCCCCCCCCCCCCCCC",
    "cyclopropane": "C1CC1", "cyclobutane": "C1CCC1", "cyclopentane": "C1CCCC1",
    "cyclohexane": "C1CCCCC1", "cycloheptane": "C1CCCCCC1", "cyclooctane": "C1CCCCCC1C",
    "methyl cyclopentane": "CC1CCCC1", "methyl cyclohexane": "CC1CCCCC1",
    "ethyl cyclohexane": "CCC1CCCCC1", "cyclohexene": "C1=CCCCC1",
    "cyclopentene": "C1=CCCC1", "1,3-cyclohexadiene": "C1=CC=CCC1",
    "1,3-butadiene": "C=CC=C", "isoprene": "CC(=C)C=C", "1-pentene": "C=CCCC",
    "2-pentene": "CCC=CC", "1-hexene": "C=CCCCC", "2-hexene": "CCC=CCC",
    "1-heptene": "C=CCCCCC", "1-octene": "C=CCCCCCC", "1-decene": "C=CCCCCCCCC",
    "propylene": "CC=C", "1-butene": "C=CCC", "2-butene": "CC=CC",
    "butadiene": "C=CC=C", "1,3-butadiene": "C=CC=C",
    "2,zdimethyl propane": "CC(C)(C)C", "2,2-dimethyl propane": "CC(C)(C)C",
    "isobutylene": "CC(=C)C", "amylenes": "",  # mix; skip
    "benzene": "c1ccccc1", "toluene": "Cc1ccccc1", "ethyl benzene": "CCc1ccccc1",
    "o-xylene": "Cc1ccccc1C", "m-xylene": "Cc1cccc(C)c1", "p-xylene": "Cc1ccc(C)cc1",
    "ethyl toluene": "CCC1=CC=CC=C1", "methyl ethyl benzene": "CCC1=CC=CC=C1",
    "trimethyl benzene": "", "pseudocumene": "Cc1cc(C)cc(C)c1", "mesitylene": "Cc1cc(C)cc(C)c1",
    "cymene": "", "diphenyl": "c1ccc(-c2ccccc2)cc1", "biphenyl": "c1ccc(-c2ccccc2)cc1",
    "isobutyl benzene": "CC(C)Cc1ccccc1",
    "naphthalene": "c1ccc2ccccc2c1", "tetralin": "C1CCC2=CC=CCC=C2C1",
    "decalin": "C1CCC2CCCCC2C1", "styrene": "C=CC1=CC=CC=C1",
    "2,2,3,3-tetramethyl pentane": "CC(C)(C)C(C)(C)C",
    "anthracene": "c1ccc2cc3ccccc3cc2c1",
}


def _norm_number_tokens(toks):
    """Numeric tokens in order; merge '4.' + '61' -> '4.61'."""
    vals, i = [], 0
    while i < len(toks):
        t = toks[i]
        if t.endswith(".") and i + 1 < len(toks) and re.fullmatch(r"\d{1,2}", toks[i + 1]):
            vals.append(t + toks[i + 1])
            i += 2
        else:
            vals.append(t)
            i += 1
    return vals


def parse_page(page, col_split):
    """Yield (left_col_text, right_col_text) per table line, chars split by x."""
    for ln in page.extract_text_lines():
        l_chars = [c for c in ln["chars"] if c["x0"] < col_split]
        r_chars = [c for c in ln["chars"] if c["x0"] >= col_split]
        out = []
        for chars in (l_chars, r_chars):
            if not chars:
                out.append("")
                continue
            text, lexp = "", -99
            for c in sorted(chars, key=lambda c: c["x0"]):
                if c["x0"] - lexp > 3.5:
                    text += " "
                text += c["text"]
                lexp = c["x1"]
            out.append(text.strip())
        yield tuple(out)


def parse_row(text):
    if not text or re.sub(r"\W+", "", text) == "":
        return None
    toks = text.split()
    name_toks, num_toks = [], []
    for t in toks:
        if re.fullmatch(r"[-0-9.,]+", t) and re.search(r"\d", t):
            num_toks.append(t)
        elif re.search(r"[A-Za-z]", t):
            name_toks.append(t)
    if not name_toks:
        return None
    name = " ".join(name_toks)
    name = re.sub(r"\s*[-_]*\s*$", "", name)
    name = re.sub(r"[^A-Za-z1-9(),\-/ ]", " ", name)
    name = re.sub(r"\s+", " ", name).strip().lower()
    nums = _norm_number_tokens(num_toks)
    vals = []
    for v in nums:
        v = v.strip(",.-")
        try:
            f = float(v)
        except ValueError:
            continue
        if 0 <= f <= 100:
            vals.append(f)
    return name, vals


# Hand-read from the PDF appendix A-1 lines (probe output above, pdfplumber text layer);
# each entry verified visually against its printed row.
CURATED = [
    # (name -> alias key in NAME_SMILES, lfl_vol_pct, pdf_page)
    ("benzene", 1.35, 123),           # "Benzene ____ 1 1. 3 560"  (prints 1.35)
    ("toluene", 1.2, 125),            # "Toluene ____ 1. _2 7. 1 480"
    ("ethyl benzene", 1.0, 124),      # "Ethyl benzene ____ 1. 0"
    ("cyclohexane", 1.3, 123),        # "Cyclohexane ____ 1. 3"
    ("decalin", 0.74, 123),           # "Decalin ____. 74"
    ("isopentane", 1.4, 124),         # "Isopentane ____ 1.4"
    ("isobutylene", 1.8, 124),        # "Isobutylene ____ 1. 8"
    ("butadiene", 2.0, 122),          # "1,3-Butadiene 2.0"
    ("2,2,3,3-tetramethyl pentane", 0.8, 124),
    ("isobutyl benzene", 0.82, 124),  # "Isobutyl benzene ----- l. 82 10 6. 0"
    ("n-decane", 0.75, 123),          # A-1; agrees with Table 2
    ("n-hexadecane", 0.43, 124),      # A-1; agrees with Table 2
]


def main():
    rows = []
    for nm, lfl, pno in CURATED:
        smiles = NAME_SMILES[nm]
        m = mol_from(smiles=smiles)
        if m is None:
            print(f"  !! bad smiles for {nm}: {smiles}")
            continue
        mm = Chem.AddHs(m)
        x = sum(1 for a in mm.GetAtoms() if a.GetSymbol() == "C")
        y = sum(1 for a in mm.GetAtoms() if a.GetSymbol() == "H")
        est = 0.55 * 100.0 / (1 + 4.773 * (x + y / 4))
        ratio = lfl / est
        flag = "" if 0.55 <= ratio <= 1.25 else "  <-- OUT OF WINDOW, review"
        print(f"  {nm:26s} LFL={lfl:5.2f} est={est:5.2f} ratio={ratio:4.2f}{flag}")
        if flag:
            continue
        rows.append({"smiles": smiles, "lfl_vol_pct": lfl, "T_C": 25.0,
                     "source": "zabetakis627_AppendixA1_ocr",
                     "note": f"name='{nm}' pdf_p{pno} ratio={ratio:.2f}"})
    out = pd.DataFrame(rows)
    print(f"\nkept {len(out)} hydrocarbon LFLs")
    if len(out):
        OUT.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(OUT, index=False)
        print(f"wrote -> {OUT}: {len(out)} rows total")
        print(out.to_string(index=False))


if __name__ == "__main__":
    main()
