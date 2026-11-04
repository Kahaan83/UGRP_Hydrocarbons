"""
Unit and regression tests for PubChem PUG-View experimental property extractor (T08).
"""
import sys
from pathlib import Path
import pandas as pd

# Support running from workspace root or inside fuelblend_data
BASE_DIR = Path(__file__).resolve().parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

try:
    from fetch_pubchem import (
        parse_flash_point,
        parse_boiling_point,
        parse_vapor_pressure,
        parse_heat_of_combustion,
        extract_properties_from_json,
        clean_text
    )
except ImportError:
    from fuelblend_data.fetch_pubchem import (
        parse_flash_point,
        parse_boiling_point,
        parse_vapor_pressure,
        parse_heat_of_combustion,
        extract_properties_from_json,
        clean_text
    )


def test_clean_text():
    raw = "Vapor pressure 4,040 kPa (-1.5 \u00b0C) & 5.21X10+4 mm Hg"
    cleaned = clean_text(raw)
    assert "4040" in cleaned
    assert "deg C" in cleaned
    assert "5.21e4" in cleaned


def test_parse_flash_point():
    # Fahrenheit
    res = parse_flash_point("115 \u00b0F (NTP, 1992)")
    assert res is not None
    val_c, unit, method, raw_val, raw_unit = res
    assert abs(val_c - 46.11) < 0.1
    assert unit == "C"
    assert method == "unspecified"

    # Celsius with closed cup
    res = parse_flash_point("46.0 \u00b0C (114.8 \u00b0F) - closed cup")
    assert res is not None
    val_c, unit, method, raw_val, raw_unit = res
    assert val_c == 46.0
    assert method == "closed_cup"

    # Negative Celsius with c.c.
    res = parse_flash_point("-18 \u00b0C c.c.")
    assert res is not None
    val_c, unit, method, raw_val, raw_unit = res
    assert val_c == -18.0
    assert method == "closed_cup"

    # Open cup
    res = parse_flash_point("55 \u00b0C (open cup)")
    assert res is not None
    assert res[2] == "open_cup"

    # Qualitative string returns None
    assert parse_flash_point("Flammable gas") is None
    assert parse_flash_point("NA (Gas)") is None


def test_parse_boiling_point():
    # Normal BP with Fahrenheit
    res = parse_boiling_point("231.1 \u00b0F at 760 mmHg (NTP, 1992)")
    assert res is not None
    val_c, unit, p_kpa, raw_val, raw_unit = res
    assert abs(val_c - 110.61) < 0.1
    assert p_kpa is not None
    assert abs(p_kpa - 101.325) < 0.1

    # Reduced pressure
    res = parse_boiling_point("174.15 \u00b0C @ 10 mmHg")
    assert res is not None
    val_c, unit, p_kpa, raw_val, raw_unit = res
    assert val_c == 174.15
    assert abs(p_kpa - 1.333) < 0.01

    # Qualitative
    assert parse_boiling_point("Sublimes") is None


def test_parse_vapor_pressure():
    # mmHg with temperature
    res = parse_vapor_pressure("28.4 mm Hg at 25 \u00b0C")
    assert len(res) == 1
    val_kpa, t_c, raw_val, raw_unit = res[0]
    assert abs(val_kpa - 3.786) < 0.01
    assert t_c == 25.0

    # kPa direct
    res = parse_vapor_pressure("Vapor pressure, kPa at 20 \u00b0C: 3.8")
    assert len(res) == 1
    val_kpa, t_c, raw_val, raw_unit = res[0]
    assert val_kpa == 3.8
    assert t_c == 20.0

    # atm
    res = parse_vapor_pressure("8.4 atm at 70 \u00b0F")
    assert len(res) == 1
    val_kpa, t_c, raw_val, raw_unit = res[0]
    assert abs(val_kpa - 851.13) < 0.1
    assert abs(t_c - 21.11) < 0.1

    # Bracketed mmHg
    res = parse_vapor_pressure("205.0 [mmHg]")
    assert len(res) == 1
    assert abs(res[0][0] - 27.33) < 0.1

    # Semicolon multi-point
    res = parse_vapor_pressure("1 mmHg at 61.7 \u00b0F ; 2.7 mmHg at 68 \u00b0F")
    assert len(res) == 2


def test_parse_heat_of_combustion():
    # Decane: MW = 142.285 g/mol, -6778.29 kJ/mol -> 47.638 MJ/kg
    res = parse_heat_of_combustion("-6778.29 kJ/mol at 25 \u00b0C", mw=142.285)
    assert res is not None
    val_mj, t_c, raw_val, raw_unit = res
    assert abs(val_mj - 47.64) < 0.1
    assert t_c == 25.0

    # BTU/lb: -17,780 BTU/lb -> 41.35 MJ/kg
    res = parse_heat_of_combustion("-17,780 BTU/lb")
    assert res is not None
    val_mj, t_c, raw_val, raw_unit = res
    assert abs(val_mj - 41.35) < 0.2

    # J/kg scientific notation
    res = parse_heat_of_combustion("-388.8X10+5 Joules/kg")
    assert res is not None
    val_mj, t_c, raw_val, raw_unit = res
    assert abs(val_mj - 38.88) < 0.1

    # Without MW for molar heat -> None
    assert parse_heat_of_combustion("-6778.29 kJ/mol", mw=None) is None


def test_output_csv_exists_and_populated():
    out_csv = Path(__file__).resolve().parent.parent / "data" / "raw" / "pubchem" / "pubchem_experimental.csv"
    assert out_csv.exists(), f"Missing {out_csv}"
    df = pd.read_csv(out_csv)
    assert len(df) >= 2000
    expected_props = {"flash_point_C", "vapor_pressure_kPa", "heat_of_combustion_MJ_kg", "boiling_point_C"}
    assert set(df["property"].unique()) == expected_props
    assert (df["tier"] == "B").all()


if __name__ == "__main__":
    test_clean_text()
    test_parse_flash_point()
    test_parse_boiling_point()
    test_parse_vapor_pressure()
    test_parse_heat_of_combustion()
    test_output_csv_exists_and_populated()
    print("All tests passed successfully!")
