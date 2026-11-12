"""v5 data audit: LFL overlay that uses the tabulated values already on disk before falling back to an estimate.

v4 state: extras/fuelblend_data/data/raw/lfl/lfl.csv holds 695 compounds, 25 measured (Zabetakis Bulletin 627) and 670 from
the Jones rule (0.55 x stoichiometric concentration).  Two tabulated sources sit unused next to it:
  * lfl_chemicals.csv  : IEC 60079-20-1 values exported from the `chemicals` library (tabulated lookups only);
  * Co-Optima database : `Pure_LFL_LEL` in raw/cooptima/fuel_engine_co_optimization.csv (compiled; origin per row not stated).

Priority per compound (first that exists): Zabetakis measured > IEC 60079-20-1 table > Co-Optima compiled > Jones estimate.
Output: data/processed/v5/lfl_v5.csv (an overlay; lfl.csv and every v4 result are left untouched).
The Co-Optima CAS is resolved to a structure with `chemicals`; a row is used only if its InChIKey is already in lfl.csv.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "extras" / "fuelblend_data" / "data" / "raw"
OUT = ROOT / "data" / "processed" / "v5" / "lfl_v5.csv"


def cooptima_lfl() -> dict[str, float]:
    from chemicals import CAS_from_any
    from chemicals.identifiers import search_chemical

    c = pd.read_csv(RAW / "cooptima" / "fuel_engine_co_optimization.csv", low_memory=False, encoding="latin-1")
    c = c[c.Pure_LFL_LEL.notna() & c.Pure_CAS.notna()]
    out: dict[str, list[float]] = {}
    for r in c.itertuples():
        try:
            smi = search_chemical(CAS_from_any(str(r.Pure_CAS).strip())).smiles
            mol = Chem.MolFromSmiles(smi)
            val = float(r.Pure_LFL_LEL)
        except Exception:
            continue
        if mol is None or not (0.1 < val < 20):
            continue
        out.setdefault(Chem.inchi.MolToInchiKey(mol), []).append(val)
    return {k: float(pd.Series(v).median()) for k, v in out.items()}


def build() -> pd.DataFrame:
    base = pd.read_csv(RAW / "lfl" / "lfl.csv")
    iec = pd.read_csv(RAW / "lfl" / "lfl_chemicals.csv")
    iec = dict(iec[iec.lfl_vol_pct.notna()][["inchikey", "lfl_vol_pct"]].itertuples(index=False, name=None))
    co = cooptima_lfl()
    rows = []
    for r in base.itertuples():
        measured = str(r.source).startswith("measured")
        if measured:
            val, src = r.lfl_vol_pct, "zabetakis627_measured"
        elif r.inchikey in iec:
            val, src = iec[r.inchikey], "iec60079_20_1_table"
        elif r.inchikey in co:
            val, src = co[r.inchikey], "cooptima_compiled"
        else:
            val, src = r.lfl_vol_pct, "jones_estimate"
        rows.append(dict(inchikey=r.inchikey, smiles=r.smiles, lfl_vol_pct=float(val), lfl_T_C=r.T_C, lfl_source_v5=src,
                         lfl_v4=float(r.lfl_vol_pct), lfl_source_v4=r.source, rel_change=float(val) / float(r.lfl_vol_pct) - 1.0,
                         iec_value=iec.get(r.inchikey), cooptima_value=co.get(r.inchikey)))
    return pd.DataFrame(rows)


def main() -> None:
    d = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT, index=False)
    print(d.lfl_source_v5.value_counts().to_string())
    ch = d[d.lfl_source_v5.isin(["iec60079_20_1_table", "cooptima_compiled"])]
    print(f"replaced estimates: {len(ch)}; mean rel change {ch.rel_change.mean():+.3f}; mean |rel change| {ch.rel_change.abs().mean():.3f}; "
          f"max |rel change| {ch.rel_change.abs().max():.3f}")
    both = d[d.iec_value.notna() & d.cooptima_value.notna()]
    print(f"IEC vs Co-Optima where both exist: n={len(both)}, mean |rel diff| {((both.cooptima_value - both.iec_value) / both.iec_value).abs().mean():.3f}")


if __name__ == "__main__":
    main()
