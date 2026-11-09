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
        "system_id": ["fuel:A", "fuel:A", "fuel:A", "fuel:B", "fuel:C", "fuel:C", "fuel:D",
                      "fuel:E", "fuel:E", "fuel:F", "fuel:F"],
        "inchikey": ["K1", "K2", "LUMP_X", "K1", "K1", "K1", "K3", "K1", "K2", "K1", "LUMP_Y"],
        "smiles": ["C", "CC", "", "C", "C", "C", "CCC", "C", "CC", "C", ""],
        "fraction": [0.4, 0.5, 0.1, 1.0, 0.97, 0.03, 1.0, 0.5, 0.5, 0.8, 0.2],
        "basis": "mass", "is_lumped": [False, False, True, False, False, False, False,
                                       False, False, False, True],
        "unexplained_fraction": 0.0,
    })
    led = pd.DataFrame([
        ("A", "fraction:K1#0", "0.4", "raw_file", "raw"),
        ("A", "fraction:K2#0", "0.5", "cited_typed", "cited"),
        ("A", "fraction:LUMP_X#0", "0.1", "raw_file", "raw"),
        ("B", "fraction:K1#0", "1.0", "target_derived_composition", "t"),
        ("C", "fraction:K1#0", "0.97", "formula_derived", "reconstructed: crc"),
        ("C", "fraction:K1#1", "0.03", "formula_derived", "other formula"),
        ("D", "fraction:K3#0", "1.0", "raw_file", "raw"),
        ("E", "fraction:K1#0", "0.5", "raw_file", "raw"),
        ("E", "fraction:K2#0", "0.5", "unknown", "u"),
        ("F", "fraction:K1#0", "0.8", "raw_file", "raw"),
        ("F", "fraction:LUMP_Y#0", "0.2", "raw_file", "raw"),
    ], columns=["fuel_id", "field", "value", "provenance_class", "evidence"])
    led.insert(0, "table", "c.csv")
    mism = pd.DataFrame({"table": ["c.csv"], "row_id": [6], "kind": ["key_mismatch"]})  # fuel D row
    inv = pd.DataFrame({"fuel_id": list("ABCDEF"), "fuel_family": "jet", "ron": [90.0, 91.0, 92.0, 93.0, 94.0, 95.0],
                        "unresolved_fraction": [0.0, 0.0, 0.0, 0.0, 0.0, 0.2],
                        "unresolved_handling": ["fully_resolved"] * 5 + ["lumped_heavy_tail"],
                        "available_properties": "ron;mon"})
    inv_led = pd.DataFrame([
        ("A", "ron", "90", "raw_file", "r"), ("B", "ron", "91", "raw_file", "r"),
        ("C", "ron", "92", "unknown", "u"), ("D", "ron", "93", "raw_file", "r"),
        ("A", "claim:ron", "listed", "raw_file", "r"), ("A", "claim:mon", "listed", "unknown", "u"),
        ("E", "ron", "94", "raw_file", "r"), ("F", "ron", "95", "raw_file", "r"),
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
    assert sorted(dropped) == ["B", "D", "E"]
    assert list(inv_v4["fuel_id"]) == ["A", "C", "F"]
    # cited_typed rows are unverified; raw_file verified
    v = dict(zip(comp_v4["provenance"], comp_v4["verified"]))
    assert not v["cited_typed"]
    assert v["raw_file"]
    # quarantine reasons
    reasons = " ".join(comp_q["reason"])
    assert "provenance:target_derived_composition" in reasons and "provenance:formula_derived" in reasons
    assert "identity:key_mismatch" in reasons
    # no row appears in both (compare against the input rows by position-independent key)
    cols = ["system_id", "inchikey", "fraction"]
    kept_keys = set(map(tuple, comp_v4[cols].values))
    quar_keys = set(map(tuple, comp_q[cols].values))
    assert not (kept_keys & quar_keys)
    assert kept_keys | quar_keys == set(map(tuple, comp[cols].values))
    # inventory: quarantined cell is NaN, the available_properties claim list is rebuilt from kept claims
    a = inv_v4.set_index("fuel_id")
    assert a.loc["A", "ron"] == 90.0 and pd.isna(a.loc["C", "ron"])
    assert a.loc["A", "available_properties"] == "ron"
    assert a.loc["A", "n_composition_rows_quarantined"] == 0
    assert a.loc["C", "n_composition_rows_quarantined"] == 1
    # partially quarantined fuel C: unresolved share is 1 - kept sum (not NaN), never verified
    assert abs(a.loc["C", "unresolved_fraction"] - 0.03) < 1e-9
    assert a.loc["C", "unresolved_handling"] == "partial_after_quarantine"
    assert not a.loc["C", "verified"]
    # dropped fuels' otherwise-clean values go to quarantine with the drop reason
    assert (inv_q[inv_q["fuel_id"] == "B"]["reason"] == Q.NO_CLEAN).all()


def test_completeness_gate():
    """R15: kept sum < MIN_KEPT_FRACTION -> dropped, listed with the sum, rows go to quarantine."""
    inv, comp, led, inv_led, mism = _synthetic()
    ann = Q.split_composition(comp, led, mism, "c.csv")
    inv_v4, comp_v4, inv_q, comp_q, dropped = Q.build_v4(inv, ann, inv_led)
    assert Q.MIN_KEPT_FRACTION == 0.95
    assert "incomplete composition after quarantine (kept fraction=0.500)" == dropped["E"]
    assert "E" not in set(inv_v4["fuel_id"]) and "E" not in set(comp_v4["system_id"].str.replace("fuel:", ""))
    # E's raw_file row is clean but is quarantined with the incomplete reason (partition stays exact)
    e = comp_q[comp_q["system_id"] == "fuel:E"]
    assert len(e) == 2 and (e["reason"].str.contains("incomplete_composition|provenance:unknown")).all()
    assert (e["reason"] == Q.INCOMPLETE).sum() == 1
    assert (inv_q[inv_q["fuel_id"] == "E"]["reason"] == Q.INCOMPLETE).all()
    # lumped fuel F: resolved 0.8 + kept lump 0.2 = 1.0 -> kept, no row lost, unresolved = 1 - kept = 0, verified
    f = inv_v4.set_index("fuel_id").loc["F"]
    assert f["unresolved_handling"] == "lumped_heavy_tail" and f["unresolved_fraction"] == 0.0
    assert f["verified"]
    # a lumped fuel whose lump is quarantined and whose resolved sum is < 0.95 is dropped
    led2 = led.copy()
    led2.loc[led2["field"] == "fraction:LUMP_Y#0", "provenance_class"] = "unknown"
    ann2 = Q.split_composition(comp, led2, mism, "c.csv")
    assert "F" in Q.classify_fuels(set(inv["fuel_id"]), ann2)[1]


def test_cited_typed_inventory_value_unverified():
    inv, comp, led, inv_led, mism = _synthetic()
    inv_led.loc[0, "provenance_class"] = "cited_typed"
    ann = Q.split_composition(comp, led, mism, "c.csv")
    inv_v4, *_ = Q.build_v4(inv, ann, inv_led)
    a = inv_v4.set_index("fuel_id").loc["A"]
    assert not a["verified"] and a["unverified_fields"] == "ron"


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
    assert (~c.loc[c["provenance"] == "cited_typed", "verified"].astype(bool)).all()
    assert c.loc[c["provenance"] == "raw_file", "verified"].astype(bool).all()
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


def test_real_gate_and_completeness(real):
    inv = real["inv_v4"]
    assert (inv["kept_composition_fraction"] >= Q.MIN_KEPT_FRACTION).all()
    assert inv["unresolved_fraction"].notna().all()
    # every fuel with a quarantined composition row is flagged unverified
    assert not inv.loc[inv["n_composition_rows_quarantined"] > 0, "verified"].astype(bool).any()
    # the 8 partially-quarantined fuels found in review are dropped, with their kept share in the reason
    eight = ("SURR_Mueller_4comp", "SURR_IDEAL_4comp", "SURR_Peugeot_4comp", "SURR_EN590_European_Diesel",
             "ASTM_D975_ULSD_Reference", "BrunoHuber2010_3comp", "Chen2016_JI", "Honnet2009_Aachen")
    for f in eight:
        assert f in real["dropped"] and "incomplete composition" in real["dropped"][f]
    # NIST_RP1 keeps its legitimate lumped tail
    assert "NIST_RP1" in set(inv["fuel_id"])


def test_real_quarantine_compositions_deduped(real):
    cq = real["comp_quar"]
    assert "also_in" in cq.columns
    # every distinct quarantined value of every source table survives: nothing but exact duplicates collapsed
    key = ["system_id", "smiles", "inchikey", "fraction"]
    out = set(map(tuple, cq[key].fillna("").astype(str).values))
    for comp_p in [p[1] for p in Q.PAIRS]:
        ann_q = real["comp_quar_by_table"][comp_p]
        src = set(map(tuple, ann_q[key].fillna("").astype(str).values))
        assert src <= out, (comp_p, sorted(src - out)[:3])
    # and exact duplicates did collapse
    assert len(cq) < sum(len(v) for v in real["comp_quar_by_table"].values())
    # NASA_LFP3 decane: differing fractions across tables are kept as separate rows
    d = cq[(cq["system_id"] == "fuel:NASA_LFP3") & (cq["smiles"] == "CCCCCCCCCC")]
    assert sorted(d["fraction"]) == [0.023, 0.025]


def test_dedupe_keeps_rows_with_different_fraction():
    a = pd.DataFrame({"table": ["t1", "t2", "t3"], "row_id": [0, 0, 0], "system_id": "fuel:X",
                      "smiles": "C", "inchikey": "K", "fraction": [0.023, 0.03, 0.03], "reason": "r",
                      "evidence": ["e1", "e2", "e3"]})
    out = Q.dedupe_quarantined(a)
    assert sorted(out["fraction"]) == [0.023, 0.03]
    assert out.set_index("fraction").loc[0.03, "also_in"] == "t2;t3"
