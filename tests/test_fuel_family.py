"""Unit tests for fuel family classification (T1.4)."""
from fuelblend.features.fuel_family import classify_mixture


def test_classify_real_jet():
    res = classify_mixture(
        system_id="fuel:12",
        n_comp=0,
        comp_inchikeys=[],
        properties=["density_kg_m3", "viscosity_kin_cSt"],
        comp_c_nums=[],
        source="crc647",
    )
    assert res["fuel_family"] == "jet"
    assert res["real_fuel"] is True


def test_classify_gasoline_mon():
    res = classify_mixture(
        system_id="34822b5cc316",
        n_comp=23,
        comp_inchikeys=["dummy"] * 23,
        properties=["mon"],
        comp_c_nums=[7, 8, 8, 9] * 5,
        source="leenhouts",
    )
    assert res["fuel_family"] == "gasoline"
    assert res["real_fuel"] is False


def test_classify_diesel_dcn():
    res = classify_mixture(
        system_id="d0893f5b8daf",
        n_comp=4,
        comp_inchikeys=["k1", "k2", "k3", "k4"],
        properties=["dcn"],
        comp_c_nums=[10, 12, 14, 16],
        source="leenhouts",
    )
    assert res["fuel_family"] == "diesel"
    assert res["real_fuel"] is False


def test_classify_model_binary():
    res = classify_mixture(
        system_id="sys_binary_hep_oct",
        n_comp=2,
        comp_inchikeys=["k1", "k2"],
        properties=["density_kg_m3", "viscosity_dyn_mPas"],
        comp_c_nums=[7, 8],
        source="thermoml",
    )
    assert res["fuel_family"] == "model-mixture"
    assert res["real_fuel"] is False
