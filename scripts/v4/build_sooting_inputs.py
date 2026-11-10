"""Task 14: assemble tracked sooting inputs from the git-ignored extras/ raw files.

Provenance classes written to the `provenance` column:
  raw_file     number read from a downloaded raw file (YSI .tab, AJFTD export)
  cited_typed  hand-typed in extras/fuelblend_data/build_tsi_table.py with a
               literature citation, NOT checkable against any file in this repo
Rows of extras tsi_table.csv with derived=True (166+9 rows, smoke point computed from
YSI by a hand-set interpolation + TSI inversion) are CIRCULAR and are excluded.
"""
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
EX = ROOT / "extras" / "fuelblend_data" / "data"
OUT = ROOT / "data" / "processed" / "v4" / "sooting"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tab = pd.read_csv(EX / "raw/sooting/YSI_Database_Volume_2.tab", sep="\t")
    ysi = pd.read_csv(EX / "raw/sooting/ysi_compounds.csv")
    ysi = ysi[ysi.is_hc].copy()
    # traceability: every stored YSI must equal UnifiedYSI in the raw Dataverse table
    raw = tab.groupby(["CAS#", "SMILES"]).UnifiedYSI.apply(list).to_dict()
    chk = ysi.apply(lambda r: bool(np.isclose(raw.get((r.cas, r.smiles_raw), [np.nan]), r.ysi).any()), axis=1)
    assert chk.mean() > 0.99, f"YSI mismatch vs raw tab: {1-chk.mean():.3f}"
    ysi = ysi[chk]
    ysi["mw"] = [Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in ysi.smiles]
    ysi["provenance"] = "raw_file"
    ysi["source_file"] = "Harvard Dataverse YSI Database Volume 2 (UnifiedYSI), Das et al. 2017/2018"
    ysi = ysi.groupby("inchikey", as_index=False).agg(
        smiles=("smiles", "first"), name=("name", "first"), ysi=("ysi", "median"),
        ysi_unc=("ysi_unc", "median"), mw=("mw", "first"), provenance=("provenance", "first"),
        source_file=("source_file", "first"))
    ysi.to_csv(OUT / "ysi_hc.csv", index=False)

    tsi = pd.read_csv(EX / "processed/tsi_table.csv")
    tier = tsi[~tsi.derived].copy()
    tier = tier.merge(ysi[["inchikey", "ysi"]].rename(columns={"ysi": "ysi_raw"}), on="inchikey", how="left")
    tier["mw"] = [Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in tier.smiles]
    tier["provenance"] = "cited_typed"
    tier = tier[["inchikey", "smiles", "compound_name", "mw", "smoke_point_mm", "censored", "ysi_raw",
                 "source", "provenance", "note"]].rename(columns={"ysi_raw": "ysi"})
    tier.to_csv(OUT / "smoke_point_pure_citedtyped.csv", index=False)

    a = pd.read_csv(EX / "processed/smoke_point_ajftd.csv")
    a = a.assign(provenance="raw_file", censored=a.value > 50.0).drop(columns=["ysi", "derived"])
    a.to_csv(OUT / "smoke_point_ajftd.csv", index=False)

    # broader sensitivity set: all AJFTD fuels with D1322 smoke point + total aromatics
    f = pd.read_csv(EX / "processed/ajftd_fuels.csv")
    c = pd.read_csv(EX / "processed/ajftd_composition.csv")
    sp = f[(f.property == "smoke_point_mm") & (f.method == "D1322")]
    sp = sp.groupby("fuel_id", as_index=False).value.median()
    c = c[["fuel_id", "comp::aromatics_total[% vol]", "comp::diaromatics[% vol]"]].drop_duplicates("fuel_id")
    c.columns = ["fuel_id", "aromatics_pct_vol", "diaromatics_pct_vol"]
    b = sp.merge(c, on="fuel_id")
    b = b[(b.value > 5) & (b.value < 50) & b.aromatics_pct_vol.between(0, 60)]
    b["provenance"] = "raw_file"
    b.to_csv(OUT / "smoke_point_ajftd_broad.csv", index=False)
    print(len(ysi), len(tier), len(a), len(b))


if __name__ == "__main__":
    main()
