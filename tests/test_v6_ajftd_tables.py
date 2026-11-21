"""v6 AJFTD tables (scripts/v6/data/build_ajftd_v6.py): re-scrape adds fuels without changing v5 rows."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
V5 = ROOT / "data" / "real_fuels" / "v5"
V6 = ROOT / "data" / "real_fuels" / "v6"

_spec = importlib.util.spec_from_file_location("build_ajftd_v6", ROOT / "scripts" / "v6" / "data" / "build_ajftd_v6.py")
b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b)

needs_v6 = pytest.mark.skipif(not (V6 / "ajftd_gcxgc_fuels.csv").exists(), reason="v6 AJFTD tables not built")


@pytest.mark.parametrize("label,expected", [
    ("POSF10264", "10264"), ("10264", "10264"), ("CRC003 POSF3962 JETA", "3962"), ("posf 4111", "4111"),
    ("(02240-070616)T115", None), ("AJ820130002", None), ("", None), (None, None),
])
def test_posf_number(label, expected):
    assert b.posf_number(label) == expected


def test_fuel_groups_pairs_same_posf_and_keeps_others_apart():
    wide = pd.DataFrame({"fuel_id": [32, 24450, 7], "fuel_label": ["10264", "POSF10264", "AJ820130002"]})
    gc = pd.DataFrame({"fuel_id": [32, 99], "fuel_label": ["10264", "(02240-070616)T115"]})
    g = b.fuel_groups(wide, gc).set_index("fuel_id")
    assert g.loc[32, "dup_group"] == g.loc[24450, "dup_group"] == "posf:10264"
    assert g.loc[7, "dup_group"] == "id:7" and g.loc[99, "dup_group"] == "id:99"
    assert g.loc[32, "n_ids_in_group"] == 2 and g.loc[7, "n_ids_in_group"] == 1


@needs_v6
def test_v5_rows_preserved_in_v6():
    for name in ["ajftd_fuels_long.csv", "ajftd_gcxgc_composition.csv"]:
        old = pd.read_csv(V5 / name, low_memory=False)
        new = pd.read_csv(V6 / name, low_memory=False)
        b._rows_preserved(old, new, name)


@needs_v6
def test_new_gcxgc_fuels_have_property_flags_from_v6_wide():
    f = pd.read_csv(V6 / "ajftd_gcxgc_fuels.csv").set_index("fuel_id")
    wide = pd.read_csv(V6 / "ajftd_fuels_wide.csv", low_memory=False).set_index("fuel_id")
    for fid in [20, 21, 22, 32, 64, 154, 155]:        # measured freezing in the 2026-11-20 re-scrape
        assert f.loc[fid, "has_freezing_point_C"]
        assert pd.notna(wide.loc[fid, "freezing_point_C"])
    for fid in [362, 363, 364]:                       # composition-only pages, leaf sum ~80 %
        assert not f.loc[fid, "composition_complete"]
        assert not f.loc[fid, "has_freezing_point_C"]
    # flags agree with the wide table for every fuel
    fr = wide["freezing_point_C"].notna()
    assert (f["has_freezing_point_C"] == f.index.map(lambda i: bool(fr.get(i, False)))).all()


@needs_v6
def test_groups_file_links_duplicate_ids():
    g = pd.read_csv(V6 / "ajftd_fuel_groups.csv", dtype={"posf": str}).set_index("fuel_id")
    assert g.loc[32, "dup_group"] == g.loc[24450, "dup_group"]
    assert g.loc[24401, "dup_group"] == g.loc[24445, "dup_group"]
    assert (g["n_ids_in_group"] <= 2).all()
