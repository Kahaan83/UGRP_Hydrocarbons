"""Tests for scripts/v4/parse_crc647_g1.py (CRC 647 Table G.1 hydrocarbon species parser)."""
import importlib.util
import pathlib

import pytest

SCRIPT = pathlib.Path(__file__).parents[1] / "scripts/v4/parse_crc647_g1.py"
spec = importlib.util.spec_from_file_location("parse_crc647_g1", SCRIPT)
g1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g1)


def test_merge_joins_split_qualifier():
    toks = [("21.61", 255.7, 274.5), ("<", 304.9, 309.3), ("0.1", 309.2, 319.7), ("21.61", 350.2, 369.0)]
    assert g1._merge_fragments(toks) == [("21.61", 255.7, 274.5), ("<0.1", 304.9, 319.7), ("21.61", 350.2, 369.0)]


def test_merge_joins_split_number():
    toks = [("18.12", 255.7, 274.5), ("2", 305.0, 309.2), (".70", 309.1, 319.6), ("20.83", 350.2, 369.0)]
    assert g1._merge_fragments(toks) == [("18.12", 255.7, 274.5), ("2.70", 305.0, 319.6), ("20.83", 350.2, 369.0)]


def test_merge_never_joins_across_columns():
    toks = [("35.0", 446.7, 461.3), ("2050", 492.9, 509.5)]  # adjacent columns, ~31 px gap
    assert g1._merge_fragments(toks) == toks


def test_wide_treats_censored_as_zero_by_default():
    rows = [dict(fuel_code="001", property="a", value_raw="<0.1", value_std="0.1"),
            dict(fuel_code="001", property="b", value_raw="2.70", value_std="2.7")]
    assert g1.wide(rows) == {"001": {"a": 0.0, "b": 2.7}}
    assert g1.wide(rows, censored_as=0.05) == {"001": {"a": 0.05, "b": 2.7}}


def test_parse_real_pdf_end_to_end():
    pdfplumber = pytest.importorskip("pdfplumber")
    pdf_path = next((p for p in g1.PDF_CANDIDATES if p.exists()), None)
    if pdf_path is None:
        pytest.skip("CRC 647 PDF not present")
    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        for pno in sorted(g1.PAGES):
            rows += g1.parse_page(pdf.pages[pno - 1], pno)
    assert len(rows) == 1026
    assert len({r["fuel_code"] for r in rows}) == 57
    raw = {}
    for r in rows:
        raw.setdefault(r["fuel_code"], []).append(r["value_raw"])
    # order per fuel: p125 (6), p126 (6), p127 (6)
    assert raw["002"] == ["19.29", "1.88", "21.17", "78.83", "35.0", "1200", "60.9", "8.3", "8.7", "0.9", "14.5", "4.7",
                          "<0.2", "<0.2", "1.6", "<0.2", "<0.2", "<0.2"]
    assert raw["014"][:6] == ["18.12", "2.70", "20.83", "79.17", "35.0", "2050"]
    assert raw["020"][:12] == ["14.10", "<0.1", "14.10", "85.90", "6.2", "154", "40.2", "10.7", "30.6", "4.2", "8.7", "5.4"]
    assert raw["508"][:12] == ["5.60", "<0.1", "5.60", "94.40", "1.4", "8.5", "92.5", "<0.2", "2.6", "<0.2", "3.0", "2.6"]
    assert raw["721"] == ["21.61", "<0.1", "21.61", "78.39", "1.2", "166", "49.9", "22.0", "6.3", "<0.2", "13.7", "7.9"] + ["<0.2"] * 6
    cens = [r for r in rows if r["fuel_code"] == "721" and r["property"] == "diaromatics_vol_pct"][0]
    assert cens["value_raw"] == "<0.1" and cens["value_std"] == "0.1"
    names = {r["fuel_code"]: r["fuel_name"] for r in rows}
    assert names["018"] == "Seattle/Boeing, WA"
