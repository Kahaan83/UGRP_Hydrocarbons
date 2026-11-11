"""v5: NREL FuelLib loader.  Fixture tests/data/fuellib/ = real FuelLib files (gc rows sliced, property files whole)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.data import loaders_fuellib as L

FIX = Path(__file__).parent / "data" / "fuellib"
REAL = L.RAW


@pytest.fixture(scope="module")
def comp():
    return L.build_compositions(FIX, {})


@pytest.fixture(scope="module")
def prop():
    return L.build_properties(FIX, {}, fuels=["posf10325", "posf11498"])


def _v(df, fuel, name, T=None):
    d = df[(df.fuel_id == fuel) & (df.property == name)]
    if T is not None:
        d = d[np.isclose(d.T_C, T)]
    assert len(d) == 1, (fuel, name, T, len(d))
    return float(d.value.iloc[0])


def test_schema_matches_v4_prefix(comp, prop):
    v4c = pd.read_csv(Path(__file__).parents[1] / "data/real_fuels/v4/compositions_v4.csv", nrows=1)
    assert list(comp.columns[:9]) == list(v4c.columns)
    v4p = pd.read_csv(Path(__file__).parents[1] / "data/real_fuels/v4/jet_properties_v4.csv", nrows=1)
    assert list(prop.columns[:12]) == list(v4p.columns)


def test_composition_basics(comp):
    assert set(comp.fuel_id) == {"FUELLIB_posf10325", "FUELLIB_posf11498", "FUELLIB_toluene"}
    for _, g in comp.groupby("fuel_id"):
        assert abs(g.mass_fraction_norm.sum() - 1) < 1e-9 and abs(g.mole_fraction.sum() - 1) < 1e-9
        assert g.inchikey.str.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]").all()
    t = comp[comp.fuel_id == "FUELLIB_toluene"].iloc[0]
    assert t.inchikey == "YXFVVABEGXRONW-UHFFFAOYSA-N" and not t.is_lumped
    assert set(comp[comp.fuel_id == "FUELLIB_posf11498"].hc_class) == {"iso_paraffin"}
    assert comp[comp.fuel_id == "FUELLIB_posf10325"].is_lumped.all()


def test_composition_spot(comp):
    g = comp[comp.fuel_id == "FUELLIB_posf10325"].set_index("gc_bin")
    assert g.loc["C4-Benzene", "fraction"] == pytest.approx(0.0332)
    assert g.loc["C4-Benzene", "n_carbon"] == 10 and g.loc["C4-Benzene", "hc_class"] == "aromatic"
    # mole fraction of C2-benzene: (1.1/106.167) / sum(w/M)
    w = g.weight_pct_raw / g.weight_pct_raw.sum()
    x = (w / g.mw) / (w / g.mw).sum()
    assert g.loc["C2-Benzene", "mole_fraction"] == pytest.approx(x["C2-Benzene"])
    assert g.loc["C2-Benzene", "mw"] == pytest.approx(106.168, abs=0.01)
    assert (g.v4_overlap == "POSF-10325").all()


def test_property_spot_10_values(prop):
    f, g = "FUELLIB_posf10325", "FUELLIB_posf11498"
    assert _v(prop, f, "viscosity_kin_cSt", -40) == pytest.approx(9.2)
    assert _v(prop, f, "density_kg_m3", 20) == pytest.approx(801.7)
    assert _v(prop, f, "flash_point_C") == 48
    assert _v(prop, f, "freezing_point_C") == -51
    assert _v(prop, f, "dcn") == pytest.approx(48.3)
    assert _v(prop, f, "distillation_T50_C") == pytest.approx(205.4)
    assert _v(prop, f, "vapor_pressure_kPa", 40) == pytest.approx(0.821955603)
    assert _v(prop, g, "viscosity_kin_cSt", -40) == pytest.approx(9.89)
    assert _v(prop, g, "flash_point_C") == pytest.approx(49.5)
    assert _v(prop, g, "vapor_pressure_kPa", 0) == pytest.approx(0.01 * L.PSI_TO_KPA)


def test_tier_and_double_count_flags(prop):
    fp = prop[~prop.hefa_wt_pct.notna()]
    assert (fp.data_tier == "B").all() and (fp.provenance == L.PROVENANCE).all()
    assert (fp.double_count_group == "AFRL-TR-2020-0017").all()
    assert (prop[prop.hefa_wt_pct.notna()].T_note != "").all()
    assert set(fp[fp.source_page.str.contains("Fig")].digitisation) == {"figure_digitised"}
    assert prop.value.notna().all()


def test_hefa_blends(prop):
    h = prop[(prop.hefa_wt_pct == 100) & (prop.fuel_id.str.contains("tall"))]
    assert _v(h, h.fuel_id.iloc[0], "density_kg_m3") == pytest.approx(757.3)
    assert _v(h, h.fuel_id.iloc[0], "viscosity_kin_cSt") == pytest.approx(4.973684211)


@pytest.mark.skipif(not (REAL / "gcData").exists(), reason="raw FuelLib data not present")
def test_real_full_load():
    c = L.build_compositions(REAL, {})
    assert c.fuel_id.nunique() == 16
    assert (c.groupby("fuel_id").mass_fraction_norm.sum().sub(1).abs() < 1e-9).all()
    p = L.build_properties(REAL, {})
    assert {"FUELLIB_posf10264", "FUELLIB_posf10289", "FUELLIB_posf10325", "FUELLIB_posf11498"} <= set(p.fuel_id)
