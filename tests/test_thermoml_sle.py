"""ThermoML solid-liquid equilibrium extractor, checked on a REAL record block.

Fixture tests/data/thermoml_sle_benzene_cyclohexane_fluid2007_01_029.xml is trimmed (header, two Compound blocks, three
PureOrMixtureData blocks) from ThermoML.v2020-09-30.tgz, file 10.1016/j.fluid.2007.01.029.xml (Tsuji et al. 2007,
benzene + cyclohexane): block 8 = 9 liquidus points (liquid mole fraction of benzene, crystal phase = benzene),
block 6 = eutectic temperature (232.11 K), block 7 = eutectic composition (x benzene = 0.2709).
"""
import collections
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from fuelblend.data import thermoml_sle as ts

FIX = Path(__file__).resolve().parent / "data" / "thermoml_sle_benzene_cyclohexane_fluid2007_01_029.xml"
BENZENE_IK, CYCLOHEXANE_IK = "UHOVQNZJYSORNB-UHFFFAOYSA-N", "XDTMQSROBMDMFD-UHFFFAOYSA-N"


@pytest.fixture(scope="module")
def parsed():
    root = ET.parse(FIX).getroot()
    st, drops = collections.Counter(), []
    comps, doi, rows = ts.parse_file(root, "fix.xml", st, 1000.0, drops)
    out = ts.to_output(rows, comps, st)
    return comps, doi, rows, out, drops, st


def test_doi_and_counts(parsed):
    _, doi, rows, out, drops, _ = parsed
    assert doi == "10.1016/j.fluid.2007.01.029"
    assert len(out) == 10 and not drops
    assert sorted(r["point_kind"] for r in out).count("liquidus") == 9


def test_liquidus_points_use_liquid_composition_and_crystal_component(parsed):
    out = parsed[3]
    liq = [r for r in out if r["point_kind"] == "liquidus"]
    assert json.loads(liq[0]["names"]) == ["benzene", "cyclohexane"]
    xs = [json.loads(r["x_liquid"]) for r in liq]
    assert xs[0] == pytest.approx([0.349, 0.651]) and xs[-1] == pytest.approx([1.0, 0.0])
    assert [r["T_K"] for r in liq][:2] == [238.42, 244.53] and liq[-1]["T_K"] == 278.59   # 278.59 K = benzene melting point
    assert all(r["P_kPa"] == 101 and r["prop_phase"] == "Liquid" and r["n_comp"] == 2 for r in liq)
    assert all(r["solid_phase_desc"] == "Crystal: benzene" and json.loads(r["crystal_inchikeys"]) == [BENZENE_IK] for r in liq)
    assert all(abs(sum(json.loads(r["x_liquid"])) - 1) < 1e-9 for r in liq)


def test_eutectic_paired_with_composition_and_flagged(parsed):
    eut = [r for r in parsed[3] if r["point_kind"] == "eutectic"]
    assert len(eut) == 1
    e = eut[0]
    assert e["T_K"] == 232.11 and e["P_kPa"] is None                 # no pressure stated
    assert json.loads(e["x_liquid"]) == pytest.approx([0.2709, 0.7291])
    assert e["solid_phase_desc"].startswith("eutectic (Crystal: benzene; cyclohexane")


def test_identity_columns(parsed):
    r = parsed[3][0]
    assert r["is_hc"] is True and json.loads(r["smiles"]) == ["c1ccccc1", "C1CCCCC1"]
    assert json.loads(r["inchikeys"]) == [BENZENE_IK, CYCLOHEXANE_IK]
    assert r["system_id"] == hashlib.md5("|".join(sorted([BENZENE_IK, CYCLOHEXANE_IK])).encode()).hexdigest()[:12]
    assert len({x["system_id"] for x in parsed[3]}) == 1


def test_pressure_filter_and_mass_fraction_conversion():
    root = ET.parse(FIX).getroot()
    st, drops = collections.Counter(), []
    _, _, rows = ts.parse_file(root, "fix.xml", st, 50.0, drops)       # P = 101 kPa > 50 -> liquidus rows dropped
    assert [r["kind"] for r in rows] == ["eutectic"] and len(drops) == 9
    comps, _, rows = ts.parse_file(root, "fix.xml", st, 1000.0, [])
    liq = [r for r in rows if r["kind"] == "liquidus"][0]
    liq["scale"] = "mass fraction"                                      # same numbers read as mass fractions
    x = json.loads(ts.to_output([liq], comps, st)[0]["x_liquid"])
    assert x[0] > 0.349 and abs(sum(x) - 1) < 1e-9                      # benzene (78 g/mol) is lighter than cyclohexane (84), so its mole fraction exceeds its mass fraction
