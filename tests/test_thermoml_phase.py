"""Regression: ThermoML phase_of() must see eVarPhase / eConstraintPhase (real records), not only ePhase.

Fixture tests/data/thermoml_vle_gas_comp_fluid2005_03_032.xml is a REAL record block (10.1016/j.fluid.2005.03.032,
methane + butane, vapour pressure at 310.95 K) trimmed from ThermoML.v2020-09-30.tgz to its one PureOrMixtureData
block.  Its only composition variable is the GAS-phase mole fraction of methane (eVarPhase = Gas), so a correct parser
must ignore it; the old phase_of() read only 'ePhase', returned None, and treated it as a liquid composition.
"""
import collections
import xml.etree.ElementTree as ET
from pathlib import Path

from fuelblend.data import thermoml_extract as tx

FIX = Path(__file__).resolve().parent / "data" / "thermoml_vle_gas_comp_fluid2005_03_032.xml"


def _blocks():
    root = ET.parse(FIX).getroot()
    return root, tx.kids(root, "PureOrMixtureData")[0]


def test_phase_of_reads_evarphase():
    _, d = _blocks()
    phases = [tx.phase_of(v, "VarPhaseID") for v in tx.kids(d, "Variable")]
    assert phases == ["Liquid", "Gas"]


def test_vapour_phase_composition_not_used_as_liquid():
    root, _ = _blocks()
    st = collections.Counter()
    _, rows = tx.parse_file(root, "fix.xml", tx.WANTED, 1000.0, st)
    # gas-phase y(methane)=0.736 must NOT be completed into a liquid binary composition
    assert rows == []
    assert st["mixture: unsupported/ambiguous composition scale"] == 1
