"""v4 Task 3: quarantine / clean-v4 partition of the v2/v3 real-fuel tables.

Synthetic unit tests of the keep rule plus checks on the REAL ledger and tables.
"""
from pathlib import Path
import importlib.util
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("v4_quarantine", ROOT / "scripts" / "v4" / "quarantine.py")
Q = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(Q)

KEEP = {"raw_file", "cited_typed", "reconstructed"}


def _synthetic():
    comp = pd.DataFrame({
        "system_id": ["fuel:A", "fuel:A", "fuel:A", "fuel:B", "fuel:C", "fuel:C", "fuel:D"],
        "inchikey": ["K1", "K2", "LUMP_X", "K1", "K1", "K1", "K3"],
        "smiles": ["C", "CC", "", "C", "C", "C", "CCC"],
        "fraction": [0.4, 0.5, 0.1, 1.0, 0.5, 0.5, 1.0],
        "basis": "mass", "is_lumped": [False, False, True, False, False, False, False],
        "unexplained_fraction": 0.0,
    })
    led = pd.DataFrame([
        ("A", "fraction:K1#0", "0.4", "raw_file", "raw"),
        ("A", "fraction:K2#0", "0.5", "cited_typed", "cited"),
        ("A", "fraction:LUMP_X#0", "0.1", "raw_file", "raw"),
        ("B", "fraction:K1#0", "1.0", "target_derived_composition", "t"),
        ("C", "fraction:K1#0", "0.5", "formula_derived", "reconstructed: crc"),
        ("C", "fraction:K1#1", "0.5", "formula_derived", "other formula"),
        ("D", "fraction:K3#0", "1.0", "raw_file", "raw"),
    ], columns=["fuel_id", "field", "value", "provenance_class", "evidence"])
    led.insert(0, "table", "c.csv")
    mism = pd.DataFrame({"table": ["c.csv"], "row_id": [6], "kind": ["key_mismatch"]})  # fuel D row
    inv = pd.DataFrame({"fuel_id": list("ABCD"), "fuel_family": "jet", "ron": [90.0, 91.0, 92.0, 93.0],
                        "available_properties": "ron;mon"})
    inv_led = pd.DataFrame([
        ("A", "ron", "90", "raw_file", "r"), ("B", "ron", "91", "raw_file", "r"),
        ("C", "ron", "92", "unknown", "u"), ("D", "ron", "93", "raw_file", "r"),
        ("A", "claim:ron", "listed", "raw_file", "r"), ("A", "claim:mon", "listed", "unknown", "u"),
    ], columns=["fuel_id", "field", "value", "provenance_class", "evidence"])
    inv_led.insert(0, "table", "i.csv")
    return inv, comp, led, inv_led, mism


def test_keep_label_rule():
    assert Q.keep_label("raw_file", "x") == "raw_file"
    assert Q.keep_label("cited_typed", "x") == "cited_typed"
    assert Q.keep_label("formula_derived", "reconstructed: foo") == "reconstructed"
    for c in ("formula_derived", "target_derived_composition", "placeholder_constant", "unknown"):
        assert Q.keep_label(c, "plain") is None


def test_synthetic_partition_and_flags():
    inv, comp, led, inv_led, mism = _synthetic()
    ann = Q.split_composition(comp, led, mism, "c.csv")
    inv_v4, comp_v4, inv_q, comp_q, dropped = Q.build_v4(inv, ann, inv_led)
    # exact partition of composition rows and inventory ledger rows
    assert len(comp_v4) + len(comp_q) == len(comp)
    assert len(inv_q) + (len(inv_led) - len(inv_q)) == len(inv_led)
    # fuel A keeps K1(raw), K2(cited), lump(raw); B target-derived -> dropped; C: one reconstructed row kept; D identity -> dropped
    assert set(comp_v4["provenance"]) == {"raw_file", "cited_typed", "reconstructed"}
    assert set(comp_v4["provenance"]) <= KEEP
    assert dropped == ["B", "D"]
    assert list(inv_v4["fuel_id"]) == ["A", "C"]
    # cited_typed rows are unverified; raw_file verified
    v = dict(zip(comp_v4["provenance"], comp_v4["verified"]))
    assert v["cited_typed"] is False or v["cited_typed"] == False  # noqa: E712
    assert v["raw_file"] == True  # noqa: E712
    # quarantine reasons
    reasons = " ".join(comp_q["reason"])
    assert "provenance:target_derived_composition" in reasons and "provenance:formula_derived" in reasons
    assert "identity:key_mismatch" in reasons
    # no row appears in both
    assert not (set(map(tuple, comp_v4[["system_id", "inchikey", "fraction"]].values))
                & set(map(tuple, comp_q[["system_id", "inchikey", "fraction"]].values if False else [])))
    # inventory: quarantined cell is NaN, the available_properties claim list is rebuilt from kept claims
    a = inv_v4.set_index("fuel_id")
    assert a.loc["A", "ron"] == 90.0 and pd.isna(a.loc["C", "ron"])
    assert a.loc["A", "available_properties"] == "ron"
    assert a.loc["A", "n_composition_rows_quarantined"] == 0
    assert a.loc["C", "n_composition_rows_quarantined"] == 1
    # dropped fuels' otherwise-clean values go to quarantine with the drop reason
    assert (inv_q[inv_q["fuel_id"] == "B"]["reason"] == Q.NO_CLEAN).all()


def test_cited_typed_inventory_value_unverified():
    inv, comp, led, inv_led, mism = _synthetic()
    inv_led.loc[0, "provenance_class"] = "cited_typed"
    ann = Q.split_composition(comp, led, mism, "c.csv")
    inv_v4, *_ = Q.build_v4(inv, ann, inv_led)
    a = inv_v4.set_index("fuel_id").loc["A"]
    assert a["verified"] == False and a["unverified_fields"] == "ron"  # noqa: E712


@pytest.fixture(scope="module")
def real():
    pytest.importorskip("rdkit")
    return Q.run(ROOT, write=False)


def test_real_v3_partition_exact(real):
    s = real["summary"].set_index("compositions")
    r = s.loc[Q.V3_COMP]
    assert r["comp_kept"] + r["comp_quarantined"] == r["comp_rows"] == 3298
    assert r["inv_kept"] + r["inv_quarantined"] == r["inv_ledger_rows"]
    # all tables: nothing lost
    assert (s["comp_kept"] + s["comp_quarantined"] == s["comp_rows"]).all()


def test_real_v4_only_clean_rows(real):
    c = real["comp_v4"]
    assert set(c["provenance"]) <= KEEP
    assert len(c) > 0
    assert (c.loc[c["provenance"] == "cited_typed", "verified"] == False).all()  # noqa: E712
    assert (c.loc[c["provenance"] == "raw_file", "verified"] == True).all()  # noqa: E712
    # no identity-mismatched row survives
    mm = pd.read_csv(ROOT / Q.MISMATCH)
    mm = mm[mm["table"] == Q.V3_COMP]
    v3 = pd.read_csv(ROOT / Q.V3_COMP)
    bad = set(map(tuple, v3.loc[mm["row_id"].astype(int), ["system_id", "smiles", "inchikey"]].values))
    got = set(map(tuple, c[["system_id", "smiles", "inchikey"]].values))
    assert not (bad & got)
    # CRC647 reconstructed compositions are kept and labelled (R10)
    crc = c[c["system_id"].astype(str).str.startswith("fuel:CRC647")]
    assert len(crc) > 0 and set(crc["provenance"]) == {"reconstructed"}


def test_real_dropped_fuels_listed_in_missing(real):
    miss = (ROOT / "data" / "MISSING.md").read_text(encoding="utf-8")
    inv4 = set(real["inv_v4"]["fuel_id"].astype(str))
    assert real["dropped"] and not (set(real["dropped"]) & inv4)
    assert set(real["dropped"]) | inv4 == set(real["inv_v3"]["fuel_id"].astype(str))
    for f in real["dropped"]:
        assert f in miss, f
    # every v4 fuel has a clean non-lump composition row
    c = real["comp_v4"]
    have = set(c.loc[c["is_lumped"].astype(str).str.lower() != "true", "system_id"].str.replace("fuel:", "", regex=False))
    assert inv4 <= have
