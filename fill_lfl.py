#!/usr/bin/env python
"""Build data/raw/lfl/lfl.csv (lower flammable limit, vol % in air) for every hydrocarbon in the unified registry.

    python fill_lfl.py            # reads data/unified/compounds.csv, writes data/raw/lfl/lfl.csv
    python build_unified.py       # picks lfl.csv up (adds lfl_vol_pct, lfl_T_C, lfl_source to compounds.csv)

Priority for each compound (first hit wins):
  1. YOUR measured values in data/raw/lfl/lfl_measured.csv   (columns: smiles OR inchikey, lfl_vol_pct, T_C, source)
  2. Zabetakis, Bureau of Mines Bulletin 627 (1965), Table 2, n-alkanes C1-C16, L25 (vol %, 25 C, 1 atm).
     Values were read from the OCR text layer of the scanned bulletin; two cells (nonane, decane) were garbled
     ('1.85', '2.75') and read as 0.85 and 0.75 after cross-checking against the bulletin's own rule L25 ~ 107/M
     (eq. 28) and against 0.55*Cst.  -> please confirm against page 21 of the PDF.
  3. ESTIMATE  LFL = 0.55 * Cst  (Jones rule, quoted in the bulletin as eq. 21 for paraffins), with
     Cst = 100 / (1 + 4.773 * (x + y/4)) vol % for CxHy.  For the n-alkanes C6-C16 this reproduces the measured
     values to within ~1-8 %; for aromatics/cycloalkanes/alkenes expect ~10-20 %.
     Effect on flash point: a 15 % LFL error moves the predicted flash point by roughly 2 K.
Every row records its `source` so the baselines can separate measured from estimated.
"""
from __future__ import annotations
from pathlib import Path

import pandas as pd
from rdkit import Chem

from common import mol_from

UNI, LFL_DIR = Path("data/unified"), Path("data/raw/lfl")

# n-alkane: (carbon number) -> L25 vol %   [Zabetakis 627, Table 2; MW 16.04 ... 226.44]
ZABETAKIS_T2 = {1: 5.0, 2: 3.0, 3: 2.1, 4: 1.8, 5: 1.4, 6: 1.2, 7: 1.05, 8: 0.95, 9: 0.85, 10: 0.75,
                11: 0.68, 12: 0.60, 13: 0.55, 14: 0.50, 15: 0.46, 16: 0.43}


def n_alkane_key(n: int) -> str:
    return Chem.MolToInchiKey(Chem.MolFromSmiles("C" * n))


def cst_vol_pct(mol) -> float:
    """Stoichiometric fuel concentration in air (vol %) for a hydrocarbon CxHy."""
    m = Chem.AddHs(mol)
    x = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
    y = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "H")
    return 100.0 / (1.0 + 4.773 * (x + y / 4.0))


def main():
    cp = UNI / "compounds.csv"
    if not cp.exists():
        raise SystemExit("run build_unified.py first (it writes data/unified/compounds.csv)")
    comp = pd.read_csv(cp)
    comp = comp[comp.is_hc.astype(bool)]

    measured: dict[str, tuple[float, float, str]] = {}
    mp = LFL_DIR / "lfl_measured.csv"
    if mp.exists():
        for r in pd.read_csv(mp).to_dict("records"):
            key = r.get("inchikey") if isinstance(r.get("inchikey"), str) else None
            if key is None and isinstance(r.get("smiles"), str):
                m = mol_from(smiles=r["smiles"])
                key = Chem.MolToInchiKey(m) if m is not None else None
            if key and pd.notna(r.get("lfl_vol_pct")):
                measured[key] = (float(r["lfl_vol_pct"]), float(r["T_C"]) if pd.notna(r.get("T_C")) else 25.0,
                                 str(r.get("source", "user")))
    table2 = {n_alkane_key(n): (v, n) for n, v in ZABETAKIS_T2.items()}

    rows = []
    for r in comp.itertuples():
        if r.inchikey in measured:
            v, T, src = measured[r.inchikey]
            rows.append((r.inchikey, r.smiles, v, T, f"measured: {src}", ""))
        elif r.inchikey in table2:
            v, n = table2[r.inchikey]
            rows.append((r.inchikey, r.smiles, v, 25.0, "measured: Zabetakis627_Table2 (OCR-read; verify p.21)", f"n-C{n}"))
        else:
            m = mol_from(smiles=r.smiles)
            v = 0.55 * cst_vol_pct(m)
            rows.append((r.inchikey, r.smiles, round(v, 4), 25.0, "estimated: 0.55*Cst (Jones rule)", r.cls))
    out = pd.DataFrame(rows, columns=["inchikey", "smiles", "lfl_vol_pct", "T_C", "source", "note"])
    LFL_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(LFL_DIR / "lfl.csv", index=False)
    kinds = out.source.str.split(":").str[0].value_counts().to_dict()
    print(f"wrote {LFL_DIR/'lfl.csv'}: {len(out)} compounds  {kinds}")
    print("To replace an estimate with a measured value add a row to data/raw/lfl/lfl_measured.csv and rerun.")


if __name__ == "__main__":
    main()
