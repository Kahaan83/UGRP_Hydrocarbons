"""v5: CRC AV-33-22 (World Jet Fuel Survey Part 1) loader.  tests/data/crc_av3322_sample.xlsx is a SYNTHETIC 3-fuel
workbook mirroring the real sheet layout (transposed 'Data' sheet, same row labels/units); the real-file tests skip if the
raw workbook is absent."""
import math
from pathlib import Path

import pandas as pd
import pytest

from fuelblend.data import loaders_crc_av3322 as L

FIX = Path(__file__).parent / "data" / "crc_av3322_sample.xlsx"


@pytest.fixture(scope="module")
def fx():
    fuels, rows = L.read_sheet(FIX)
    log = {}
    comp = L.build_compositions(fuels, rows, log)
    prop = L.build_properties(fuels, rows, log)
    return fuels, rows, comp, prop, log


def get(prop, p, fid, T_C=None, tag=None):
    m = prop[(prop.property == p) & (prop.fuel_id == fid)]
    if tag:
        m = m[m.source_tag == tag]
    if T_C is not None:
        m = m[(m.T_K - (T_C + 273.15)).abs() < 1e-6]
    assert len(m) == 1, (p, fid, T_C, len(m))
    return m.iloc[0]


def test_fuels_and_types(fx):
    fuels, *_ = fx
    assert [f["fuel_id"] for f in fuels] == ["CRC_AV3322_01", "CRC_AV3322_02", "CRC_AV3322_03"]
    assert [f["fuel_type"] for f in fuels] == ["Jet A-1", "Jet A", "Jet A"]


def test_composition_bins_and_flags(fx):
    _, _, comp, _, log = fx
    assert log["n_bins"] == 6                                   # the 4 non-bin 'Hydrocarbon Type' rows are excluded
    assert not comp.isna().any().any() and len(comp) == 18
    one = comp[comp.fuel_id == "CRC_AV3322_01"].set_index("bin_label")
    assert one.loc["Alkylbenzene-C16+", ["group", "carbon_number", "plus", "lower"]].tolist() == ["alkylbenzene", 16, True, False]
    assert one.loc["n-alkane-C03 and lower", ["group", "carbon_number", "plus", "lower"]].tolist() == ["n-alkane", 3, False, True]
    assert one.loc["Monocycloalkane-C21+", "group"] == "monocycloalkane"
    assert one.loc["n-alkane-C10", "mass_pct"] == 20.0
    assert "Hydrogen content" not in set(comp.bin_label)
    assert list(comp.columns) == L.COMP_COLS


def test_viscosity_density_temperatures(fx):
    _, _, _, prop, _ = fx
    v = prop[(prop.source_tag == "udri_measured") & (prop.property == "viscosity_kin_cSt")]
    assert v.groupby("fuel_id").T_K.nunique().eq(2).all() and (v.value > 0).all()
    assert get(prop, "viscosity_kin_cSt", "CRC_AV3322_02", -40, "udri_measured").value == 12.15
    assert get(prop, "density_kg_m3", "CRC_AV3322_01", 20, "udri_measured").value == 796.0


def test_vapor_pressure_psi_to_kpa(fx):
    _, _, _, prop, _ = fx
    r = get(prop, "vapor_pressure_kPa", "CRC_AV3322_01", 100)
    assert math.isclose(r.value, 2.17 * 6.894757293168, rel_tol=1e-9) and r.unit_std == "kPa" and r.raw_unit == "psi"


def test_surface_tension_and_cp_units(fx):
    _, _, _, prop, _ = fx
    assert math.isclose(get(prop, "surface_tension_N_m", "CRC_AV3322_01", 25).value, 0.0255)
    assert get(prop, "heat_capacity_kJ_kg_K", "CRC_AV3322_01", 20).value == 1.86


def test_distillation_naming(fx):
    _, _, _, prop, _ = fx
    assert get(prop, "distillation_T0_C", "CRC_AV3322_01").value == 142
    assert get(prop, "distillation_T100_C", "CRC_AV3322_02").value == 294
    assert get(prop, "distillation_T10_C", "CRC_AV3322_01").value == 165
    assert get(prop, "distillation_residue_pct", "CRC_AV3322_01").value == 1.1
    # SimDis '0.05%' means 5 %, '0.1%' means 10 %; the D86 *correlation* row is not loaded
    assert get(prop, "distillation_simdis_T5_C", "CRC_AV3322_01").value == 135
    assert get(prop, "distillation_simdis_T10_C", "CRC_AV3322_01").value == 151
    assert get(prop, "distillation_simdis_T0_C", "CRC_AV3322_01").value == 100
    assert get(prop, "distillation_simdis_T100_C", "CRC_AV3322_02").value == 333
    assert not prop.row_label.str.contains("correlation").any()


def test_coa_rows_censoring_and_tags(fx):
    _, _, _, prop, log = fx
    f = get(prop, "freezing_point_C", "CRC_AV3322_02")
    assert f.qualifier == "<" and f.value == -60.0 and f.raw_value == "<-60" and f.source_tag == "supplier_coa"
    assert get(prop, "freezing_point_C", "CRC_AV3322_01").qualifier == ""
    assert log["censored_rows"] == [("CRC_AV3322_02", "freezing_point_C", "<-60")]
    assert get(prop, "flash_point_C", "CRC_AV3322_02").method == "D56"
    assert get(prop, "viscosity_kin_cSt", "CRC_AV3322_01", -20, "supplier_coa").value == 3.05
    assert get(prop, "density_kg_m3", "CRC_AV3322_01", 15, "supplier_coa").value == 799.8
    assert not ((prop.property == "density_kg_m3") & (prop.source_tag == "supplier_coa") & (prop.fuel_id == "CRC_AV3322_03")).any()
    assert set(prop.source_tag) == {"udri_measured", "supplier_coa"}
    assert (prop.double_count_group == "CRC_AV3322").all() and (prop.data_tier == "A").all()
    assert list(prop.columns) == L.PROP_COLS


def test_parse_value():
    assert L.parse_value("<-60") == (-60.0, "<") and L.parse_value(">5") == (5.0, ">")
    assert math.isnan(L.parse_value("Not specified")[0]) and L.parse_value(3)[1] == ""


def test_bad_bin_label_raises():
    with pytest.raises(ValueError):
        L.parse_bin("Aromatic-C09")


# ------------------------------------------------------------------ real workbook
@pytest.mark.skipif(not L.RAW.exists(), reason="raw CRC AV-33-22 workbook not present (gitignored)")
class TestRealWorkbook:
    @pytest.fixture(scope="class")
    def real(self):
        fuels, rows = L.read_sheet(L.RAW)
        log = {}
        return fuels, L.build_compositions(fuels, rows, log), L.build_properties(fuels, rows, log), log

    def test_counts_and_types(self, real):
        fuels, comp, prop, log = real
        assert len(fuels) == 38 and comp.fuel_id.nunique() == 38 and log["n_bins"] == 104
        assert comp.bin_label.nunique() == 104 and len(comp) == 38 * 104
        assert pd.Series([f["fuel_type"] for f in fuels]).value_counts().to_dict() == {"Jet A": 23, "Jet A-1": 15}

    def test_composition_clean(self, real):
        _, comp, _, _ = real
        assert not comp.isna().any().any()
        s = comp.groupby("fuel_id").mass_pct.sum()
        assert s.between(99.9, 100.1).all()
        assert not comp.duplicated(["fuel_id", "group", "carbon_number", "plus", "lower"]).any()

    def test_viscosity_nine_temps(self, real):
        _, _, prop, _ = real
        v = prop[(prop.source_tag == "udri_measured") & (prop.property == "viscosity_kin_cSt")]
        assert v.groupby("fuel_id").T_K.nunique().eq(9).all() and (v.value > 0).all()
        got = sorted(round(t, 2) for t in v.T_K.unique())
        assert got == [round(t + 273.15, 2) for t in (-40, -20, 0, 15, 20, 40, 60, 80, 100)]

    def test_vapor_pressure_kpa(self, real):
        _, _, prop, _ = real
        vp = prop[prop.property == "vapor_pressure_kPa"]
        assert (vp.unit_std == "kPa").all() and (vp.raw_unit == "psi").all()
        assert ((vp.value - vp.raw_value * 6.894757293168).abs() < 1e-9).all()
        assert vp.value.max() < 60            # 100-120 C values of a few psi are tens of kPa

    def test_censored_freezing_keeps_qualifier(self, real):
        _, _, prop, _ = real
        f = prop[(prop.property == "freezing_point_C") & (prop.source_tag == "supplier_coa")]
        c = f[f.raw_value.str.startswith("<")]
        assert len(c) >= 1 and (c.qualifier == "<").all()
        assert (c.value == c.raw_value.str.lstrip("<").astype(float)).all()
        assert (f[~f.raw_value.str.startswith("<")].qualifier == "").all()
