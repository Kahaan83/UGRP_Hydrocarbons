"""v4 Task 6: jet fuels with cited class-level compositions.

Fixtures in tests/data/v4 are verbatim slices (real headers, real rows) of
data/raw/gc_composition/fuel_gc_compositions.csv, data/raw/pour_point/av2718_fuel_properties.csv and
data/raw/pour_point/crc647_fuel_properties.csv.  Real-file tests skip when the (gitignored) raw data is absent.
"""
import json
import math
from pathlib import Path

import pandas as pd
import pytest

from fuelblend.data.loaders_jet import load_gc_fuels, load_jet_properties

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "data" / "v4"
RAW = ROOT / "data" / "raw"
GC_REAL, AV_REAL, CRC_REAL = (RAW / "gc_composition" / "fuel_gc_compositions.csv",
                              RAW / "pour_point" / "av2718_fuel_properties.csv",
                              RAW / "pour_point" / "crc647_fuel_properties.csv")
GC_FIX, AV_FIX, CRC_FIX = (FIX / "gc_composition_fixture.csv", FIX / "av2718_fixture.csv",
                           FIX / "crc647_fixture.csv")
PROV = {"raw_file", "cited_typed", "formula_derived", "target_derived_composition",
        "placeholder_constant", "unknown"}
REAL = all(p.exists() for p in (GC_REAL, AV_REAL, CRC_REAL))
SOURCES = [pytest.param((GC_FIX, AV_FIX, CRC_FIX), id="fixture"),
           pytest.param((GC_REAL, AV_REAL, CRC_REAL), id="real",
                        marks=pytest.mark.skipif(not REAL, reason="raw data absent"))]
POSF = ["POSF-10325", "POSF-13511", "POSF-13516"]


@pytest.fixture(params=SOURCES)
def corpus(request):
    gc, av, crc = request.param
    return load_gc_fuels(gc, av), load_jet_properties(av, crc)


def _vec(row):
    fr = json.loads(row["class_fractions"])
    return tuple(round(fr.get(k, 0.0), 6) for k in
                 ("n_paraffin", "iso_paraffin", "paraffin_unsplit", "cycloalkane", "aromatic", "alkene_alkyne"))


def test_posf_10325_class_sums(corpus):
    comp, _ = corpus
    by = comp.set_index("fuel_id")
    # (1) class fractions of POSF-10325 sum to 1 +- 0.02.  The raw file lists the aromatics/cycloparaffin
    # totals AND their sub-totals (all rows sum to 152.28 %), so a naive sum would fail.
    r = by.loc["POSF-10325"]
    assert abs(sum(json.loads(r["class_fractions"]).values()) - 1.0) <= 0.02
    assert 0.98 <= r["class_sum_raw"] <= 1.02
    for f in POSF:
        assert abs(sum(json.loads(by.loc[f, "class_fractions"]).values()) - 1.0) <= 0.02, f
    # (2) no two POSF fuels share an identical fraction vector
    vecs = [_vec(by.loc[f]) for f in POSF]
    assert len(set(vecs)) == len(vecs), vecs


def test_cn_split_flags(corpus):
    comp, _ = corpus
    by = comp.set_index("fuel_id")
    # POSF-10325 has a measured boiling distribution for the same fuel (av2718 D2887), the others have none
    assert bool(by.loc["POSF-10325", "headline_ok"]) and by.loc["POSF-10325", "cn_split_method"].startswith("D2887")
    for f in ("POSF-13511", "POSF-13516", "JetA_API14avg", "UDRI_HTA_repfuel"):
        assert not bool(by.loc[f, "headline_ok"]), f
        assert "assumed_cn_split" in by.loc[f, "headline_reason"], f
    # ATSDR JP-5/JP-8 tables list only ~41-42 % of the fuel -> not a composition
    for f in ("JP5_PotterSimmons", "JP8_PotterSimmons"):
        if f in by.index:
            assert not bool(by.loc[f, "composition_complete"])
            assert by.loc[f, "unresolved_fraction"] > 0.5 and not bool(by.loc[f, "headline_ok"])


def test_components_consistent(corpus):
    comp, _ = corpus
    for _, r in comp.iterrows():
        comps = json.loads(r["components"])
        assert comps and all(c["basis"] == "mass" and c["provenance"] in PROV for c in comps)
        tot = sum(c["frac"] for c in comps)
        assert tot == pytest.approx(1.0 - r["unresolved_fraction"], abs=0.02) or not r["composition_complete"]
        assert r["composition_provenance"] in PROV
        for c in comps:
            assert c["cn"] is None or isinstance(c["cn"], int)
            assert (c["cn"] is None) == (c["cn_source"] == "assumed_cn_split")


def test_cn_distribution_from_same_fuel_distillation(corpus):
    comp, _ = corpus
    r = comp.set_index("fuel_id").loc["POSF-10325"]
    comps = json.loads(r["components"])
    cns = [c["cn"] for c in comps if c["cn"] is not None]
    assert cns and min(cns) >= 7 and max(cns) <= 18      # D2887 107-310 C is C7-C18
    assert {c["cn_source"] for c in comps} == {"D2887_same_fuel"}


def test_properties_only_from_raw_files(corpus):
    _, props = corpus
    assert set(props.provenance) == {"raw_file"}
    p = props[props.fuel_id == "POSF-10325"]
    assert {"density_kg_m3", "viscosity_kin_cSt", "flash_point_C"} <= set(p.property)
    assert p[p.property == "flash_point_C"].value.iloc[0] == 55.0     # av2718 (ASTM D56)
    # fuels with no local property file have NO rows (NaN by absence), never an invented value
    assert not props.fuel_id.isin(["POSF-13511", "POSF-13516", "UDRI_HTA_repfuel", "JetA_API14avg"]).any()
    assert not props.value.isna().any()


def test_crc647_values_reread_from_raw_not_placeholder(corpus):
    _, props = corpus
    crc = props[props.fuel_id.str.startswith("CRC647_")]
    pour = crc[crc.property == "pour_point_C"].value
    freeze = crc[crc.property == "freezing_point_C"].value
    assert len(pour) >= 1 and len(freeze) >= 1
    assert pour.nunique() > 1 or len(pour) < 3         # v3 stored a constant -48.0 placeholder
    assert not (freeze == -48.0).all()
    b = crc[(crc.fuel_id == "CRC647_2") & (crc.property == "pour_point_C")].value.iloc[0]
    assert b == -70.0                                  # raw crc647 row for Boston, MA
    h = crc[(crc.fuel_id == "CRC647_2") & (crc.property == "net_heat_of_combustion_MJ_kg")].value.iloc[0]
    assert h == pytest.approx(18509 * 0.0023260, rel=1e-4)       # Btu/lb -> MJ/kg
    d = crc[(crc.fuel_id == "CRC647_2") & (crc.property == "density_kg_m3") & crc.method.eq("ASTM D4052")]
    assert d.value.iloc[0] == pytest.approx(803.1) and d.T_C.iloc[0] == pytest.approx(15.56)


def test_no_jet_inventory_duplicates(corpus):
    comp, _ = corpus
    inv = pd.read_csv(ROOT / "data" / "real_fuels" / "v4" / "inventory_v4.csv") \
        if (ROOT / "data" / "real_fuels" / "v4" / "inventory_v4.csv").exists() else pd.DataFrame({"fuel_id": []})
    assert not set(comp.fuel_id) & set(inv.fuel_id)
    assert comp.fuel_id.is_unique
    assert "inventory_fuel_id" in comp.columns and comp.inventory_fuel_id.isna().all()
