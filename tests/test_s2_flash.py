"""v5 S2: mapping coverage, mass->mole conversion, Antoine range audit (scripts/v5/s2_t2_flash.py, audit_antoine_range_v5.py)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


t2 = _load("s2_t2_flash_t", "scripts/v5/s2_t2_flash.py")
aud = _load("audit_antoine_range_v5_t", "scripts/v5/audit_antoine_range_v5.py")


@pytest.fixture(scope="module")
def maps():
    table = t2.parse_rep_table()
    return table, t2.build_maps(table)


def test_rep_table_parsed_and_rule(maps):
    table, _ = maps
    assert table[("alkylbenzene", 9)][0] == "Cc1cc(C)cc(C)c1"
    assert t2.rep_for(table, "n_paraffin", 12)[3] == 0
    assert t2.rep_for(table, "n_paraffin", 30) is None          # beyond the gap limit of 3
    assert t2.rep_for(table, "cyclo_di", 16) is not None         # few-molecule class: any gap
    assert t2.rep_for(table, "aromatic_di", 11) is None          # no representative in the r3 table


def test_every_included_fuel_has_at_least_90_percent_mapped_mass(maps):
    _, (comp, amap, fl, fmap) = maps
    amp = amap[["group", "kind", "carbon_number", "rep_smiles", "mapped", "reason"]]
    fmp = fmap[["gc_bin", "kind", "n_carbon", "rep_smiles", "mapped", "reason"]]
    n_checked = 0
    for fid, g in comp.groupby("fuel_id"):
        _, _, cov = t2.fuel_mixture(g[["group", "kind", "carbon_number", "mass_pct"]], ["group", "kind", "carbon_number"], "mass_pct", 0.01, amp)
        if (1 - cov["mapped_fraction"]) <= t2.MAX_UNMAPPED:       # included fuels
            assert cov["mapped_fraction"] >= 0.90
            n_checked += 1
    for fid, g in fl.groupby("fuel_id"):
        _, _, cov = t2.fuel_mixture(g[["gc_bin", "kind", "n_carbon", "fraction"]], ["gc_bin", "kind", "n_carbon"], "fraction", 1.0, fmp)
        assert cov["mapped_fraction"] >= 0.90
        n_checked += 1
    assert n_checked >= 40            # 38 AV-33-22 + 4 FuelLib in the current data


def test_mapping_excludes_fuel_over_10_percent_unmapped():
    df = pd.DataFrame({"k": ["a", "b"], "m": [0.85, 0.15]})
    mp = pd.DataFrame({"k": ["a", "b"], "rep_smiles": ["CCCC", ""], "mapped": [True, False], "reason": ["", "no rep"]})
    smi, w, cov = t2.fuel_mixture(df, ["k"], "m", 1.0, mp)
    assert smi == ["CCCC"] and np.isclose(w.sum(), 1.0) and np.isclose(cov["mapped_fraction"], 0.85)
    assert (1 - cov["mapped_fraction"]) > t2.MAX_UNMAPPED          # would be excluded


def test_mass_to_mole_sums_to_one_and_is_correct():
    x = t2.mass_to_mole(["CCCC", "CCCCCCCCCC"], [0.5, 0.5])        # butane 58.12, decane 142.28
    assert np.isclose(x.sum(), 1.0)
    assert np.isclose(x[0], (0.5 / 58.123) / (0.5 / 58.123 + 0.5 / 142.286), atol=1e-3)
    assert x[0] > x[1]            # lighter molecule has the larger mole fraction at equal mass
    x3 = t2.mass_to_mole(["C", "CC", "CCC"], [0.2, 0.3, 0.5])
    assert np.isclose(x3.sum(), 1.0)


def test_written_t2_outputs_are_consistent():
    p = ROOT / "results/v5/S2/t2_flash_predictions.csv"
    if not p.exists():
        pytest.skip("run scripts/v5/s2_t2_flash.py first")
    d = pd.read_csv(p)
    assert (d.mapped_fraction >= 0.90).all()
    assert len(d) == 42 and set(d.source) == {"AV-33-22", "FuelLib"}


def test_audit_flags_when_T_outside_range_and_not_when_inside():
    assert aud.flag_outside(40.0, 354.64, 438.87)[0] is True        # 313 K < 354.64 K
    f, d = aud.flag_outside(40.0, 354.64, 438.87)
    assert d == pytest.approx(354.64 - 313.15)
    assert aud.flag_outside(100.0, 354.64, 438.87) == (False, 0.0)  # 373.15 K inside
    assert aud.flag_outside(200.0, 354.64, 438.87)[0] is True        # 473 K > 438.87 K
    assert aud.flag_outside(40.0, np.nan, np.nan) [0] is False       # unknown range is never flagged


def test_written_audit_matches_rule():
    p = ROOT / "results/v5/S2/antoine_range_flags_records.csv"
    if not p.exists():
        pytest.skip("run scripts/v5/audit_antoine_range_v5.py first")
    r = pd.read_csv(p)
    k = r.range_lo_K.notna()
    T = r.T_ideal_C + 273.15
    expect = k & ((T < r.range_lo_K) | (T > r.range_hi_K))
    assert (expect == r.outside_range).all()
