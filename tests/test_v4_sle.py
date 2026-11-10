"""Task 13: ThermoML solid-liquid equilibrium (--sle) extraction.

Fixture tests/data/thermoml_sle_je5001389.xml is a REAL ThermoML record (copied unmodified from
ThermoML.v2020-09-30.tgz: J. Chem. Eng. Data 2014, biphenyl + hexacosane (+ dibenzofuran) SLE, DTA).
Expected values below were read by hand from the XML: block 14 (crystal = biphenyl, org 3) has
variable 'Mole fraction' of component 1 (hexacosane) = 0.1098 with SLE T = 336.94 K, p = 100 kPa;
block 15 (crystal = hexacosane) has x(hexacosane) = 0.5023 with T = 323.78 K.
"""
import collections
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXTRAS = ROOT / "extras" / "fuelblend_data"
sys.path.insert(0, str(EXTRAS))
import pytest  # noqa: E402

try:
    import thermoml_extract as tx  # noqa: E402
    _HAS_SLE = hasattr(tx, "parse_file_sle") and hasattr(tx, "main_sle")
except ImportError:
    tx, _HAS_SLE = None, False
_SKIP = pytest.mark.skipif(not _HAS_SLE, reason="extras/fuelblend_data/thermoml_extract.py missing or lacks --sle support "
                           "(git-ignored nested repo): apply results/v4/SLE/thermoml_extract_sle.patch inside extras/fuelblend_data")

FIX = ROOT / "tests" / "data" / "thermoml_sle_je5001389.xml"


def _parse():
    st = collections.Counter()
    comps, rows = tx.parse_file_sle(ET.parse(FIX).getroot(), "fix.xml", st)
    return comps, rows, st


def _binary_hc(comps, rows):
    return [r for r in rows if len(r["ids"]) == 2 and all(tx.hc_inchi(comps[i]["inchi"]) for i in r["ids"])]


@_SKIP
def test_sle_record_parsed():
    comps, rows, st = _parse()
    hc = _binary_hc(comps, rows)
    # biphenyl (org 3) + hexacosane (org 1): 12 liquidus points in two blocks (5 with crystal=biphenyl, 7 with crystal=hexacosane)
    assert len(hc) == 12, len(hc)
    name = {i: comps[i]["name"] for i in comps}
    pt = [r for r in hc if abs(r["value"] - 336.94) < 1e-9]
    assert len(pt) == 1
    r = pt[0]
    assert r["prop"] == "Solid-liquid equilibrium temperature, K"
    assert r["presentation"] == "Direct value, X" and r["prop_phase"] == "Liquid"
    assert r["comp_type"] == "mole fraction"
    assert r["P_kPa"] == 100.0
    assert [name[i] for i in r["ids"]] == ["biphenyl", "hexacosane"]
    assert abs(r["comp"]["3"] - 0.8902) < 1e-9 and abs(r["comp"]["1"] - 0.1098) < 1e-9   # last component implied from x(C26)
    assert [name[i] for i in r["crystal_ids"]] == ["biphenyl"]
    r2 = [r for r in hc if abs(r["value"] - 323.78) < 1e-9][0]
    assert abs(r2["comp"]["1"] - 0.5023) < 1e-9 and abs(r2["comp"]["3"] - 0.4977) < 1e-9
    assert [name[i] for i in r2["crystal_ids"]] == ["hexacosane"]
    # dibenzofuran (oxygenated) ternaries must not pass the hydrocarbon filter, and unsupported solvent-ratio scales are rejected
    assert not tx.hc_inchi(comps["2"]["inchi"])
    assert all(len(r["ids"]) == 2 or not all(tx.hc_inchi(comps[i]["inchi"]) for i in r["ids"]) for r in rows)
    assert st["sle: unsupported/ambiguous composition scale"] > 0


@_SKIP
def test_default_mode_unchanged():
    # (a) whitelist untouched
    assert tx.WANTED == {
        "Mass density, kg/m3": "density", "Viscosity, Pa*s": "viscosity",
        "Kinematic viscosity, m2/s": "viscosity", "Vapor or sublimation pressure, kPa": "vapor"}
    assert "Solid-liquid equilibrium temperature, K" not in tx.WANTED
    # (b) the default parser still ignores the SLE property entirely
    st = collections.Counter()
    _, rows = tx.parse_file(ET.parse(FIX).getroot(), "fix.xml", tx.WANTED, 1000.0, st)
    assert rows == []
    # (c) the pre-existing regression suite still passes
    p = subprocess.run([sys.executable, str(EXTRAS / "tests" / "test_thermoml.py")], cwd=str(EXTRAS),
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-800:] + p.stderr[-800:]


if __name__ == "__main__":
    assert _HAS_SLE, "extras/fuelblend_data/thermoml_extract.py lacks --sle; apply results/v4/SLE/thermoml_extract_sle.patch"
    test_sle_record_parsed()
    test_default_mode_unchanged()
    print("test_v4_sle OK")
