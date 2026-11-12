"""Transcription tests for Mueller et al. 2016 diesel surrogates (numbers only; doi:10.1021/acs.energyfuels.5b02879)."""
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import MolToInchiKey

ROOT = Path(__file__).resolve().parents[1]
SURR = ROOT / "data/raw/surrogates/fuel_surrogates.csv"
PROPS = ROOT / "data/raw/surrogates/mueller2016_properties.csv"
COMPOUNDS = ROOT / "data/unified/compounds.csv"

# Table 3 wt % columns as printed (page 1449); mol % is in fuel_surrogates.csv. Palette abbreviation -> wt %.
TABLE3_WT = {
    "V0a": {"NHXD": 32.2, "HMN": 42.0, "TDEC": 10.5, "1MN": 15.3},
    "V0b": {"NOD": 32.1, "HMN": 32.8, "TMB": 8.1, "TET": 14.8, "1MN": 12.3},
    "V1": {"NHXD": 3.2, "NOD": 27.3, "HMN": 35.1, "NBCX": 3.8, "TDEC": 4.0, "TMB": 4.8, "TET": 10.8, "1MN": 10.9},
    "V2": {"NOD": 15.2, "NEI": 1.2, "2MHPD": 10.2, "NBCX": 14.8, "TIPCX": 12.8, "PHP": 6.4, "TIPB": 16.6, "TET": 12.0, "1MN": 10.9},
}
ABBREV = {  # name in fuel_surrogates.csv -> Table 1 abbreviation
    "n-hexadecane": "NHXD", "n-octadecane": "NOD", "n-eicosane": "NEI", "2,2,4,4,6,8,8-heptamethylnonane": "HMN",
    "2-methylheptadecane": "2MHPD", "n-butylcyclohexane": "NBCX", "1,3,5-triisopropylcyclohexane": "TIPCX",
    "trans-decalin": "TDEC", "perhydrophenanthrene": "PHP", "1,2,4-trimethylbenzene": "TMB",
    "1,3,5-triisopropylbenzene": "TIPB", "tetralin": "TET", "1-methylnaphthalene": "1MN",
}
# Palette compounds whose InChIKey is not in data/unified/compounds.csv (listed, never guessed)
EXPECTED_MISSING = {"PHP"}


@pytest.fixture(scope="module")
def recipes():
    d = pd.read_csv(SURR)
    d = d[d.fuel_code.str.startswith("Mueller2016_")].copy()
    d["abbr"] = d.component_name.map(ABBREV)
    return d


def test_four_surrogates_with_expected_component_counts(recipes):
    counts = recipes.groupby("fuel_code").size().to_dict()
    assert counts == {"Mueller2016_V0a": 4, "Mueller2016_V0b": 5, "Mueller2016_V1": 8, "Mueller2016_V2": 9}
    assert recipes.abbr.notna().all()
    assert (recipes.frac_type == "mole").all()
    assert recipes.reference.str.contains("10.1021/acs.energyfuels.5b02879").all()
    assert recipes.reference.str.contains("Table 3").all()


def test_mole_and_mass_fractions_sum_to_100(recipes):
    for code, g in recipes.groupby("fuel_code"):
        assert abs(g.frac_value.sum() - 100) <= 0.5, code
    for v, d in TABLE3_WT.items():
        assert abs(sum(d.values()) - 100) <= 0.5, v


def test_mass_fractions_consistent_with_mole_fractions(recipes):
    for code, g in recipes.groupby("fuel_code"):
        v = code.split("_")[1]
        mw = {r.abbr: Descriptors.MolWt(Chem.MolFromSmiles(r.component_smiles)) for r in g.itertuples()}
        n = {r.abbr: r.frac_value for r in g.itertuples()}
        tot = sum(n[a] * mw[a] for a in n)
        assert set(n) == set(TABLE3_WT[v])
        for a in n:
            assert abs(100 * n[a] * mw[a] / tot - TABLE3_WT[v][a]) <= 1.0, (v, a)  # paper rounding; tolerance 1.0 wt %


def test_smiles_parse_and_inchikeys_in_registry(recipes):
    reg = set(pd.read_csv(COMPOUNDS).inchikey)
    missing = set()
    for r in recipes.drop_duplicates("abbr").itertuples():
        m = Chem.MolFromSmiles(r.component_smiles)
        assert m is not None, r.abbr
        if MolToInchiKey(m) not in reg:
            missing.add(r.abbr)
    assert missing == EXPECTED_MISSING


def test_three_rows_fixture_as_printed(recipes):
    g = recipes[recipes.fuel_code == "Mueller2016_V0a"].set_index("abbr").frac_value
    assert (g["NHXD"], g["HMN"], g["TDEC"], g["1MN"]) == (27.8, 36.3, 14.8, 21.1)
    p = pd.read_csv(PROPS)
    f = p[(p.fuel_code == "Mueller2016_V2") & (p.property == "flash_point_C")].iloc[0]
    assert (f.value, f.table, f.page) == (80, "Table 5", 1457)
    c = p[(p.fuel_code == "Mueller2016_CF") & (p.property == "cloud_point_C")].iloc[0]
    assert (c.value, c.table, c.page) == (-19.2, "Figure 14", 1456)


def test_properties_measured_only_columns():
    p = pd.read_csv(PROPS)
    assert list(p.columns) == ["fuel_code", "property", "value", "unit", "T_C", "method", "table", "page"]
    assert set(p.property) == {"flash_point_C", "density_kg_m3", "cloud_point_C"}
    assert p.method.str.contains("10.1021/acs.energyfuels.5b02879").all()
    assert (p[p.property == "density_kg_m3"].T_C == 20).all()
