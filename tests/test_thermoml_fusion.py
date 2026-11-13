"""thermoml_fusion on a trimmed REAL record (10.1016/j.fluid.2017.03.012: n-eicosane, n-tetracosane, biphenyl, 1-octadecanol, ...)."""
import collections
import xml.etree.ElementTree as ET
from pathlib import Path

from fuelblend.data import thermoml_fusion as F

FIX = Path(__file__).parent / "data" / "thermoml_fusion_fluid2017_03_012_trimmed.xml"


def _rows():
    stats = collections.Counter()
    rows = F.parse_file(ET.parse(FIX).getroot(), FIX.name, stats)
    return rows, stats


def test_phase_rule():
    assert F.is_fusion_phases(["Crystal", "Liquid", "Air at 1 atmosphere"])
    assert F.is_fusion_phases(["Crystal 1", "Liquid", "Gas"])
    assert not F.is_fusion_phases(["Crystal 2", "Crystal 1", "Air at 1 atmosphere"])
    assert not F.is_fusion_phases(["Crystal 1", "Nematic liquid crystal"])
    assert not F.is_fusion_phases(["Crystal", "Gas"])


def test_solid_solid_rows_dropped_and_mixtures_skipped():
    rows, stats = _rows()
    assert all(r["phases"].count("Crystal") == 1 and "Liquid" in r["phases"] for r in rows)
    assert len(rows) == 6              # eicosane, tetracosane, 1-octadecanol x (Tm, dHfus); 4 solid-solid datasets dropped
    assert any("not exactly one Crystal" in k for k in stats)
    # tetracosane (org 2): one Tm and one fusion enthalpy kept; its 'triple point' and enthalpy between Crystal 2 and Crystal 1 are dropped
    t24 = [r for r in rows if r["org"] == "2"]
    assert sorted(r["prop"] for r in t24) == ["dhfus", "tm"]


def test_values_plausible_and_hc_filter():
    rows, _ = _rows()
    by = {(r["name"], r["prop"]): r["value"] for r in rows}
    assert 308 < by[("eicosane", "tm")] < 312
    assert 60 < by[("eicosane", "dhfus")] < 75
    df = __import__("pandas").DataFrame(rows)
    df["inchikey"] = "x"
    # summarize works on the raw rows (hydrocarbon filter is applied in extract(); 1-octadecanol is org 4 here)
    s = F.summarize(df.assign(inchikey=df.org, smiles="C"))
    assert set(s.columns) >= {"tm_K", "dhfus_kJ_mol", "spread_tm", "spread_dh", "source"}
    assert (s.source == F.TIER).all()
