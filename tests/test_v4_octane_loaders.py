"""v4 Task 5: FACE set I/II and Co-Optima bRON loaders.

Fixtures in tests/data/v4 are verbatim copies of the real header + selected real rows
(extras/fuelblend_data/data/raw/{octane,cooptima}). Real-file tests skip when the (gitignored)
raw data is absent.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from fuelblend.data.loaders_octane import load_cooptima_bron, load_face_sets

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "data" / "v4"
RAW = ROOT / "extras" / "fuelblend_data" / "data" / "raw"
FACE_FIX = [FIX / "face_set_I_fixture.csv", FIX / "face_set_II_fixture.csv"]
BRON_FIX = FIX / "bron_fixture.xlsx"
FACE_REAL = [RAW / "octane" / "face_gasoline_set_I.csv", RAW / "octane" / "face_gasoline_set_II.csv"]
BRON_REAL = RAW / "cooptima" / "bron_data_with_cas.xls"
PROV = {"raw_file", "cited_typed", "formula_derived", "target_derived_composition",
        "placeholder_constant", "unknown"}


def _comps(row):
    return json.loads(row["components"])


@pytest.mark.parametrize("paths", [FACE_FIX, pytest.param(
    FACE_REAL, marks=pytest.mark.skipif(not all(p.exists() for p in FACE_REAL), reason="raw FACE absent"))])
def test_face_prf10_parses(paths):
    df = load_face_sets(paths)
    r = df[df.blend_id == "PRF+10"].iloc[0]
    assert r.ron == 55.6 and r.mon == 53.0          # values in the file
    comps = _comps(r)
    # The FACE CSV carries component NAMES only, never fractions: n-heptane / iso-octane are
    # recorded as unresolved members of a PRF lump, ethanol's volume share comes from "+10".
    members = {m.lower() for m in json.loads(r.unresolved_members)}
    assert {"n-heptane", "iso-octane"} <= members
    eth = [c for c in comps if c["name"].lower() == "ethanol"]
    assert len(eth) == 1 and eth[0]["frac"] == pytest.approx(0.10) and eth[0]["basis"] == "volume"
    assert sum(c["frac"] for c in comps) == pytest.approx(1.0, abs=0.01)
    assert any(str(c["inchikey"]).startswith("LUMP:") for c in comps)
    assert not r.target_ok and not r.composition_resolved
    assert r.provenance in PROV and r.composition_provenance in PROV


@pytest.mark.parametrize("paths", [FACE_FIX])
def test_face_rules(paths):
    df = load_face_sets(paths).set_index("blend_id")
    # FACE base gasoline is a lump pseudo-component, never a target
    a = _comps(df.loc["FACE A+10"])
    assert {c["inchikey"] for c in a} == {"LUMP:FACE_A", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"}
    assert next(c for c in a if c["inchikey"] == "LUMP:FACE_A")["frac"] == pytest.approx(0.9)
    assert df.loc["FACE A+10", "base_fuel"] == "FACE_A" and not df.loc["FACE A+10", "target_ok"]
    # "+3" in TPRF 4+3 is a series index, not an ethanol percentage: nothing may be invented
    t = _comps(df.loc["TPRF 4+3"])
    assert len(t) == 1 and t[0]["inchikey"].startswith("LUMP:") and t[0]["frac"] == 1.0
    # no ethanol cell and "+0" -> ethanol 0
    assert all(c["name"].lower() != "ethanol" for c in _comps(df.loc["TPRF 1+0"]))
    # TMB name with embedded commas survives
    assert "1,2,4-Trimethylbenzene" in json.loads(df.loc["TMB+10", "unresolved_members"])
    # every fraction in every row has a basis; no row claims to be a resolved target
    for _, r in df.iterrows():
        assert all(c["basis"] in ("volume", "unresolved") for c in _comps(r))
        assert not r.target_ok


@pytest.mark.parametrize("path", [BRON_FIX, pytest.param(
    BRON_REAL, marks=pytest.mark.skipif(not BRON_REAL.exists(), reason="raw bRON absent"))])
def test_bron_units_and_basis(path):
    raw = pd.read_excel(path)
    # independent expected count: distinct (compound row, vol_frac) over RON and MON data points
    exp = 0
    for _, r in raw.iterrows():
        fr = set()
        for p in ("RON", "MON"):
            for i in range(5):
                v = r.get(f"vol_frac_b{p}_{i}")
                if pd.notna(v) and pd.notna(r.get(f"blend_{p}_b{p}_{i}")):
                    fr.add(round(float(v), 6))
        exp += len(fr)
    df = load_cooptima_bron(path)
    assert len(df) == exp
    assert df.ron.dropna().between(0, 130).all() and df.mon.dropna().between(0, 130).all()
    for _, r in df.iterrows():
        comps = _comps(r)
        assert all(c["basis"] == "volume" for c in comps)
        assert sum(c["frac"] for c in comps) == pytest.approx(1.0, abs=1e-6)
        lump = [c for c in comps if str(c["inchikey"]).startswith("LUMP:")]
        assert any(c["inchikey"].startswith("LUMP:BRON_BASE") for c in lump)
        blendstock_frac = [c["frac"] for c in comps if not c["inchikey"].startswith("LUMP:BRON_BASE")][0]
        assert 0.05 < blendstock_frac < 0.55           # volume fractions, not percent
        assert r.provenance == "raw_file" and r.composition_provenance == "raw_file"
        assert not r.target_ok and r.base_lump_input_only
    # methyl furan: first RON point is 98.1 at 10.29 vol%
    mf = df[df.blend_id.str.contains("534-22-5")].sort_values("blend_id").iloc[0]
    assert mf.ron == pytest.approx(98.1) and mf.base_ron == pytest.approx(90.3)
