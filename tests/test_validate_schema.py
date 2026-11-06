"""Tests for schema validation of measurements and compositions (T1.1, T1.2)."""
import pandas as pd
import pytest

from fuelblend.data.validate import (
    validate_measurement_df,
    validate_composition_df,
    ValidationError,
    CANONICAL_PROPERTIES,
    CANONICAL_UNITS,
    REQUIRED_MEASUREMENT_PROVENANCE,
)


def test_measurement_validation_valid():
    df = pd.DataFrame([{
        "measurement_id": "m1",
        "source_id": "thermoml",
        "source_citation": "NIST TRC ThermoML Archive",
        "system_id": "530c1a396846",
        "property": "viscosity_kin_cSt",
        "value": 1.5,
        "unit": "cSt",
        "T_K": 298.15,
        "P_kPa": 101.325,
        "method": "ASTM D445",
        "composition_basis": "mass",
        "uncertainty": 0.02,
        "quality_flag": "ok",
        "license_note": "public domain",
        "raw_value": 1.5,
        "raw_unit": "cSt",
    }])
    report = validate_measurement_df(df, strict=True)
    assert report["valid"] is True
    assert len(report["errors"]) == 0


def test_measurement_validation_missing_provenance():
    df = pd.DataFrame([{
        "measurement_id": "m1",
        "source_id": "",  # missing source_id
        "source_citation": "Some citation",
        "system_id": "sys1",
        "property": "density_kg_m3",
        "value": 800.0,
        "unit": "kg/m3",
        "T_K": 293.15,
        "P_kPa": 101.325,
        "method": "",  # missing method
        "composition_basis": "mass",
        "uncertainty": None,
        "quality_flag": "ok",
        "license_note": "open",
    }])
    with pytest.raises(ValidationError) as excinfo:
        validate_measurement_df(df, strict=True)
    assert "source_id" in str(excinfo.value) or "method" in str(excinfo.value)


def test_measurement_validation_unknown_property_or_unit():
    df = pd.DataFrame([{
        "measurement_id": "m1",
        "source_id": "thermoml",
        "source_citation": "citation",
        "system_id": "sys1",
        "property": "invented_property",
        "value": 10.0,
        "unit": "invented_unit",
        "T_K": 298.15,
        "P_kPa": 101.325,
        "method": "unknown",
        "composition_basis": "mass",
        "uncertainty": None,
        "quality_flag": "ok",
        "license_note": "open",
    }])
    with pytest.raises(ValidationError) as excinfo:
        validate_measurement_df(df, strict=True)
    err = str(excinfo.value)
    assert "invented_property" in err or "invented_unit" in err


def test_composition_validation_valid():
    df = pd.DataFrame([
        {
            "system_id": "sys_bin",
            "inchikey": "DCAYPVUWAIABOU-UHFFFAOYSA-N",
            "smiles": "CCCCCCC",
            "fraction": 0.4,
            "basis": "mass",
            "is_lumped": False,
            "unexplained_fraction": 0.0,
        },
        {
            "system_id": "sys_bin",
            "inchikey": "CTQNGGLPUBDAKN-UHFFFAOYSA-N",
            "smiles": "CCCCCCCC",
            "fraction": 0.6,
            "basis": "mass",
            "is_lumped": False,
            "unexplained_fraction": 0.0,
        },
    ])
    report = validate_composition_df(df, strict=True)
    assert report["valid"] is True
    assert len(report["errors"]) == 0


def test_composition_validation_invalid_fraction_sum():
    df = pd.DataFrame([
        {
            "system_id": "sys_bad",
            "inchikey": "DCAYPVUWAIABOU-UHFFFAOYSA-N",
            "smiles": "CCCCCCC",
            "fraction": 0.3,
            "basis": "mass",
            "is_lumped": False,
            "unexplained_fraction": 0.0,
        },
        {
            "system_id": "sys_bad",
            "inchikey": "CTQNGGLPUBDAKN-UHFFFAOYSA-N",
            "smiles": "CCCCCCCC",
            "fraction": 0.5,
            "basis": "mass",
            "is_lumped": False,
            "unexplained_fraction": 0.0,
        },
    ])
    with pytest.raises(ValidationError) as excinfo:
        validate_composition_df(df, strict=True)
    assert "sys_bad" in str(excinfo.value)
    assert "sum" in str(excinfo.value).lower()


def test_composition_validation_lumped_unexplained():
    df = pd.DataFrame([
        {
            "system_id": "sys_lump",
            "inchikey": "DCAYPVUWAIABOU-UHFFFAOYSA-N",
            "smiles": "CCCCCCC",
            "fraction": 0.85,
            "basis": "mass",
            "is_lumped": False,
            "unexplained_fraction": 0.15,
        },
        {
            "system_id": "sys_lump",
            "inchikey": "LUMP_HEAVY",
            "smiles": "",
            "fraction": 0.15,
            "basis": "mass",
            "is_lumped": True,
            "unexplained_fraction": 0.15,
        },
    ])
    report = validate_composition_df(df, strict=True)
    assert report["valid"] is True
