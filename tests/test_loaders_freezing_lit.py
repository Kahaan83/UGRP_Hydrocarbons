"""Tests for the literature freezing-point / composition ingest (data/raw/freezing_lit)."""
import json

import pandas as pd
from rdkit import Chem

from fuelblend.data import loaders_freezing_lit as L


def test_row_counts():
    p = L.load_properties()
    assert p.groupby("source_id").size().to_dict() == {"vozka2018_si": 25, "ascent065A_2025_f2": 11, "ascent065A_slides2024_f3": 8,
                                                      "kittel2023": 79}
    assert (p[p.source_id == "vozka2018_si"].property == "freezing_point_C").sum() == 22
    assert (p[p.source_id == "kittel2023"].property == "freezing_point_C").sum() == 7
    assert (p.property == "freezing_point_C").sum() == 22 + 19 + 7
    assert len(L.load_vozka_gcxgc()) == 272


def test_spot_values():
    p = L.load_properties().set_index(["source_id", "sample_id", "property"]).value
    assert p[("vozka2018_si", "Jet A", "freezing_point_C")] == -51.0           # SI p9, Table S4 (rendered page checked)
    assert p[("vozka2018_si", "C-30", "freezing_point_C")] == -52.0            # p9
    assert p[("vozka2018_si", "T-60", "freezing_point_C")] == -58.0            # p9
    assert p[("vozka2018_si", "MFAT", "freezing_point_C")] == -48.5            # p10
    assert p[("ascent065A_2025_f2", "nC12", "freezing_point_C")] == -9.5       # F2 p6 Table A-1
    assert p[("ascent065A_2025_f2", "100% n-alkane fuel (C9-C16)", "freezing_point_C")] == -40.7   # F2 p10
    assert p[("ascent065A_slides2024_f3", "0.17% nC16 / 0.52% nC13 / 99.31% iC8", "freezing_point_C")] == -46.1   # F3 slide 21
    assert p[("kittel2023", "HEFA 215", "freezing_point_C")] == -49.1          # F6 p4 Table 1
    assert p[("kittel2023", "Jet HC", "smoke_point_mm")] == 18.9
    assert p[("kittel2023", "HEFA Cam", "flash_point_C")] == 43.5


def test_censoring_flags():
    p = L.load_properties()
    k = p[p.source_id == "kittel2023"]
    fz = k[k.property == "freezing_point_C"].set_index("sample_id").qualifier
    assert fz["Jet HC"] == "<" and fz["FCC HN"] == "<" and fz["Jet A-1"] == ""
    sp = k[k.property == "smoke_point_mm"].set_index("sample_id").qualifier
    assert sp["HEFA Cam"] == ">" and sp["HEFA 215"] == ">"
    assert "FCC HN" not in sp.index
    assert (p[p.source_id == "vozka2018_si"].query("kind == 'real_fuel' and components != '[]'").qualifier == "blend_basis_unconfirmed").all()


def test_smiles_parse_and_fractions():
    p = L.load_properties()
    for s in p.smiles:
        for smi in json.loads(s):
            assert Chem.MolFromSmiles(smi) is not None
    for r in p.itertuples():
        f = json.loads(r.fractions_norm)
        if r.fraction_basis != "none" and f:
            assert abs(sum(f) - 1) < 1e-3, r.sample_id
        if json.loads(r.smiles):
            assert len(json.loads(r.comp_inchikeys)) == len(f)
    named = p[p.smiles != "[]"]
    assert named.is_hc.all()
    assert named.system_id.str.len().eq(12).all()
    assert p[p.smiles == "[]"].system_id.str.startswith("fuel:").all()


def test_gcxgc_closure_and_mapping():
    g = L.load_vozka_gcxgc()
    tot = g.groupby("fuel_id").value.sum()
    assert ((tot - 100).abs() < 0.1).all()             # sums 99.95 .. 100.03 as printed
    assert set(g.hc_class) <= {"n_paraffin", "iso_paraffin", "monocycloparaffin", "dicycloparaffin", "alkylbenzene", "cycloaromatic", "diaromatic"}
    assert g.common_class.notna().all() and g.subkind.notna().all()
    from fuelblend.mixture.common import CLASSES
    assert set(g.common_class) <= set(CLASSES)
    a = pd.read_csv(L.ROOT / "data/real_fuels/v5/ajftd_gcxgc_composition.csv")
    assert set(g.hc_class) <= set(a.hc_class)           # same class vocabulary as the AJFTD GCxGC fuels
    # Jet A spot values, SI p2-3: n-C9 4.42, n-C10 4.73; CAME iC17 21.58 (p3)
    assert g[(g.fuel_id == "Jet A") & (g["class"] == "n-paraffins") & (g.n_carbon == 9)].value.iloc[0] == 4.42
    assert g[(g.fuel_id == "CAME") & (g["class"] == "isoparaffins") & (g.n_carbon == 17)].value.iloc[0] == 21.58


def test_blend_composition():
    b = L.blend_composition("C-30")
    a = L.vozka_as_ajftd_composition()
    assert abs(b.mass_pct.sum() - (0.7 * a[a.fuel_id == "Jet A"].mass_pct.sum() + 0.3 * a[a.fuel_id == "CAME"].mass_pct.sum())) < 1e-6
    assert list(L.blend_composition("Jet A").columns) == list(b.columns)
