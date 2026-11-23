"""Boehm 2022 neat-fuel cold viscosity loader (values checked against the workbook cell by cell)."""
from fuelblend.data import loaders_cold_visc_v9 as L


def _v(df, sample, prop, t):
    x = df[(df.sample_id == sample) & (df.property == prop) & (df.T_K == t)]
    assert len(x) == 1
    return x.value.iloc[0]


def test_neat_fuels_only():
    df = L.load_cold_visc_v9()
    assert len(df) == 36 and df.sample_id.nunique() == 10
    assert not df.sample_id.str.startswith("boehm_r").any()


def test_spot_values_workbook_rows_207_and_274():
    df = L.load_cold_visc_v9()
    assert abs(_v(df, "A-1 (POSF 10264)", "viscosity_kin_cSt", 253.15) - 3.4392) < 1e-9   # This Work row 207
    assert abs(_v(df, "A-1 (POSF 10264)", "density_kg_m3", 253.15) - 806.1) < 1e-9
    assert abs(_v(df, "A-2 (POSF 10325)", "viscosity_kin_cSt", 233.15) - 9.1359) < 1e-9   # row 274


def test_units_and_temperatures():
    df = L.load_cold_visc_v9()
    assert set(df.T_K) == {253.15, 233.15}
    assert set(df[df.property == "viscosity_kin_cSt"].unit) == {"cSt"}
    assert df[df.property == "density_kg_m3"].value.between(600, 1100).all()
