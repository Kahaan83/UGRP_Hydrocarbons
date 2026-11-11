"""S6 diesel tests: data integrity, gate logic, no reconstruction path (PDFs are gitignored; tests use the transcribed CSVs)."""
import importlib.util
import shutil
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("s6_diesel", ROOT / "scripts" / "v5" / "s6_diesel.py")
S6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S6)

V5 = ROOT / "data" / "real_fuels" / "v5"
pytestmark = pytest.mark.skipif(not (V5 / "diesel_face_composition.csv").exists(), reason="diesel CSVs not built")


def test_surrogate_mole_fractions_sum_to_one():
    s = S6.load_surrogates()
    assert set(s.fuel_id) == {"V0a", "V0b", "V1", "V2"}
    sums = s.groupby("fuel_id").mole_fraction.sum()
    assert ((sums - 1.0).abs() < 2e-3).all(), sums


def test_face_composition_sums_in_gate_range_and_is_table_sourced():
    c = S6.load_diesel_composition()
    assert set(c.fuel_id) == set(S6.FUELS)
    assert (c.source_kind == "table").all()
    sums = c.groupby("fuel_id").mass_pct.sum()
    assert ((sums > 99.5) & (sums < 100.5)).all(), sums
    assert c.carbon_number.between(3, 22).all()


def test_class_level_table52_row_totals_consistent():
    t = pd.read_csv(V5 / "diesel_face_class_table52.csv")
    # n-paraffins + isoparaffins = paraffins and known + unknown = total aromatics (as printed, rounding only)
    assert ((t.n_paraffins + t.isoparaffins - t.paraffins).abs() < 0.03).all()
    assert ((t.known_aromatics + t.unknown_aromatics - t.total_aromatics).abs() < 0.03).all()


def test_gate_counts_only_table_sourced_real_fuels():
    g = S6.gate_table()
    assert S6.gate_count(g) == int(g[g.fuel_id.isin(S6.FUELS)].in_gate_count.sum()) == 9
    # surrogates never count even if flagged
    g2 = g.copy()
    g2.loc[g2.fuel_id.isin(["V0a", "V0b", "V1", "V2"]), "in_gate_count"] = True
    assert S6.gate_count(g2) == 9
    # a non-table row never counts
    g3 = g.copy()
    g3.loc[g3.fuel_id == "FD1A", "composition_source_kind"] = "reconstructed"
    assert S6.gate_count(g3) == 8
    assert not g[g.fuel_id.isin(["V0a", "V0b", "V1", "V2"])].in_gate_count.any()


def test_loader_rejects_reconstructed_rows(tmp_path, monkeypatch):
    for f in ("diesel_face_composition.csv", "diesel_face_properties.csv", "diesel_avfl18a_surrogates.csv"):
        shutil.copy(V5 / f, tmp_path / f)
    c = pd.read_csv(tmp_path / "diesel_face_composition.csv")
    c.loc[0, "source_kind"] = "reconstructed"
    c.to_csv(tmp_path / "diesel_face_composition.csv", index=False)
    monkeypatch.setattr(S6, "V5", tmp_path)
    with pytest.raises(ValueError):
        S6.load_diesel_composition()
    s = pd.read_csv(tmp_path / "diesel_avfl18a_surrogates.csv")
    s.loc[0, "source_kind"] = "estimated"
    s.to_csv(tmp_path / "diesel_avfl18a_surrogates.csv", index=False)
    with pytest.raises(ValueError):
        S6.load_surrogates()


def test_no_reconstruction_path_in_script():
    src = (ROOT / "scripts" / "v5" / "s6_diesel.py").read_text(encoding="utf-8")
    # the diesel loader and scorer never call the distillation/bulk-descriptor reconstruction entry points
    for banned in ("reconstruct_fuel", "reconstruct_from", "from_distillation", "bulk_descriptor"):
        assert banned not in src


def test_property_table_units():
    p = S6.load_diesel_properties()
    fl = p[p.property == "flash_point_C"]
    assert fl.value.between(40, 90).all()                      # printed in F, stored in C
    d = p[p.property == "density_kg_m3"]
    assert d.value.between(790, 880).all()
    v = p[p.property == "viscosity_kin_cSt"]
    assert v.value.between(1.0, 4.0).all()


def test_transcription_check_has_no_mismatch_and_enough_values():
    f = ROOT / "results" / "v5" / "S6" / "transcription_check.csv"
    if not f.exists():
        pytest.skip("run scripts/v5/s6_diesel.py check first")
    t = pd.read_csv(f)
    assert t.match.all()
    assert (t.groupby("table").size() >= 8).all()               # NIST Table 11 has 8 values (all compared)
    assert (t.groupby("table").size().drop("NIST TN 1949 Table 11") >= 9).all()
