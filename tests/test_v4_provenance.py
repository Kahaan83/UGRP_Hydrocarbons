"""v4 Task 1: provenance ledger.

Every test runs the classifier on the REAL current files under data/real_fuels/
and the REAL ingest scripts (scripts/v2/ingest_real_fuels_v2.py,
scripts/v3/ingest_phase1_expansion.py). Nothing here is a synthetic fixture:
the point is to prove the audit detects the fabrication that is actually in
the repository, from code + data, without trusting labels such as
`composition_type`.
"""
from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.data.provenance import (  # noqa: E402
    CLASSES,
    LEDGER_COLUMNS,
    classify_fuel_rows,
    load_source_code,
)

RF = ROOT / "data" / "real_fuels"


def _ledger(inv_name: str, comp_name: str) -> pd.DataFrame:
    inv = pd.read_csv(RF / inv_name)
    comp = pd.read_csv(RF / comp_name)
    out = classify_fuel_rows(inv, comp, load_source_code(ROOT))
    assert list(out.columns) == LEDGER_COLUMNS
    assert set(out["provenance_class"]) <= set(CLASSES)
    assert out["evidence"].astype(str).str.len().gt(0).all()
    return out


@pytest.fixture(scope="module")
def v3():
    return _ledger("v3/inventory_v3.csv", "v3/compositions_v3.csv")


@pytest.fixture(scope="module")
def v2():
    return _ledger("inventory.csv", "compositions.csv")


@pytest.fixture(scope="module")
def recon():
    return _ledger("inventory_reconstructed.csv", "compositions_reconstructed.csv")


def _fractions(led, pattern):
    m = led["fuel_id"].str.fullmatch(pattern) & led["field"].str.startswith("fraction:")
    return led[m]


def test_ledger_flags_face_composition_as_target_derived(v3, v2):
    for led in (v3, v2):
        rows = _fractions(led, r"CRC_FACE_[A-J]")
        assert rows["fuel_id"].nunique() == 10
        bad = rows[rows["provenance_class"] != "target_derived_composition"]
        assert bad.empty, bad.to_string()
    # the evidence must come from the data relation, not from a label
    toluene = _fractions(v3, r"CRC_FACE_D")
    toluene = toluene[toluene["field"].str.contains("UAEPNZWRGJTJPN")]
    assert toluene["evidence"].str.contains("dens").all()


def test_ledger_flags_ethanol_dvpe_formula(v3):
    r = v3[(v3["fuel_id"] == "CRC_FACE_A_E10") & (v3["field"] == "dvpe_kpa")]
    assert len(r) == 1
    assert r["provenance_class"].iloc[0] == "formula_derived"
    assert "7.0" in r["evidence"].iloc[0]


def test_ledger_flags_constant_hoc(v3):
    m = (v3["fuel_id"].str.startswith("CRC_FACE_")
         & (v3["field"] == "heat_of_combustion_MJ_kg")
         & (pd.to_numeric(v3["value"], errors="coerce") == 43.5))
    rows = v3[m]
    assert len(rows) == 10
    assert (rows["provenance_class"] == "placeholder_constant").all(), rows.to_string()


def test_ledger_flags_diesel_flash_formula(v3):
    r = v3[(v3["fuel_id"] == "CRC_FACE_Diesel_1") & (v3["field"] == "flash_point_C")]
    assert len(r) == 1
    assert r["provenance_class"].iloc[0] == "formula_derived"
    assert "0.25" in r["evidence"].iloc[0]


def test_ledger_flags_posf_48_32_20(v3):
    rows = _fractions(v3, r"POSF_(10264|10289|11498|12344|12345)")
    assert rows["fuel_id"].nunique() == 5
    assert (rows["provenance_class"] == "placeholder_constant").all(), rows.to_string()
    assert rows["evidence"].str.contains("identical").all()


def test_ledger_keeps_crc647_measured(recon, v3):
    # The v2 CRC 647 inventory only *claims* freezing_point_C (no stored value). The measured values
    # exist in the PDF-ingested raw table, and the evidence says so -- but a claim without a stored,
    # matching value is never raw_file.
    r = recon[recon["fuel_id"].str.startswith("CRC647_") & (recon["field"] == "claim:freezing_point_C")]
    assert len(r) == 57
    assert (r["provenance_class"] == "unknown").all(), r.to_string()
    assert r["evidence"].str.contains("raw exists:").all()
    assert r["evidence"].str.contains("crc647_fuel_properties.csv").all()
    vc = v3[v3["fuel_id"].str.startswith("CRC647_") & (v3["field"] == "claim:freezing_point_C")]
    assert len(vc) == 57 and (vc["provenance_class"] != "raw_file").all()
    # The v3 copy of those fuels carries one constant (-48.0) instead of the
    # measured values: that is a placeholder, not a measurement.
    v = v3[v3["fuel_id"].str.startswith("CRC647_") & (v3["field"] == "freezing_point_C")]
    assert len(v) == 57
    assert (v["provenance_class"] == "placeholder_constant").all()
    assert v["evidence"].str.contains("crc647_fuel_properties.csv").all()


def test_ledger_is_deterministic_and_covers_every_row(v3):
    again = _ledger("v3/inventory_v3.csv", "v3/compositions_v3.csv")
    pd.testing.assert_frame_equal(v3.reset_index(drop=True), again.reset_index(drop=True))
    inv = pd.read_csv(RF / "v3" / "inventory_v3.csv")
    comp = pd.read_csv(RF / "v3" / "compositions_v3.csv")
    assert set(inv["fuel_id"]) <= set(v3["fuel_id"])
    assert (v3["field"].str.startswith("fraction:")).sum() == len(comp)
    assert not v3.duplicated(["fuel_id", "field"]).any()


# --------------------------------------------------------------------------- fix round 1
def test_replayed_comprehension_sees_env_variables():
    """Interpreter-independent replay: a comprehension must see the replay env (pre-3.12 scoping)."""
    import ast
    from fuelblend.data.provenance import _safe_eval
    node = ast.parse("[f * norm_factor for f in raw_fracs]", mode="eval").body
    assert _safe_eval(node, {"raw_fracs": [1.0, 2.0], "norm_factor": 0.5}) == [0.5, 1.0]


def test_data_rules_stand_without_code():
    """Detection (b) alone: no source code at all, real v3 files."""
    inv = pd.read_csv(RF / "v3" / "inventory_v3.csv")
    comp = pd.read_csv(RF / "v3" / "compositions_v3.csv")
    led = classify_fuel_rows(inv, comp, {})
    face = _fractions(led, r"CRC_FACE_[A-J]")
    assert face["fuel_id"].nunique() == 10
    assert (face["provenance_class"] == "target_derived_composition").all(), face.to_string()
    diesel = _fractions(led, r"CRC_FACE_Diesel_\d")
    assert (diesel["provenance_class"] == "target_derived_composition").all()
    hoc = led[led["fuel_id"].str.fullmatch(r"CRC_FACE_[A-J]") & (led["field"] == "heat_of_combustion_MJ_kg")]
    assert len(hoc) == 10 and (hoc["provenance_class"] == "placeholder_constant").all()
    fz = led[(led["fuel_id"] == "CRC647_335") & (led["field"] == "freezing_point_C")]
    assert fz["provenance_class"].tolist() == ["placeholder_constant"]


_SNIPPET = '''
def main():
    table = [("F_A", "fuel a", 700.0), ("F_B", "fuel b", 740.0)]
    props = {"G_A": {"dens": 705.0}, "G_B": {"dens": 745.0}}
    inv, comp = [], []
    for fuel, label, rho in table:                      # renamed property variable
        aro = min(0.45, max(0.05, (rho - 680.0) / 200.0))
        inv.append({"fuel_id": fuel, "fuel_name": label, "density_kg_m3": rho})
        comp.append({"system_id": "fuel:" + fuel, "inchikey": "TOLUENE", "smiles": "Cc1ccccc1",
                     "fraction": round(aro, 5), "basis": "mass", "is_lumped": False})
    for fuel in ["G_A", "G_B"]:
        q = props[fuel]["dens"]                         # property unpacked by subscript
        comp.append({"system_id": "fuel:" + fuel, "inchikey": "TOLUENE", "smiles": "Cc1ccccc1",
                     "fraction": round((q - 680.0) / 200.0, 5), "basis": "mass", "is_lumped": False})
'''


def test_renamed_or_unpacked_property_still_target_derived():
    inv = pd.DataFrame({"fuel_id": ["F_A", "F_B", "G_A", "G_B"], "fuel_name": ["fuel a", "fuel b", "g", "g"],
                        "fuel_family": "gasoline", "density_kg_m3": [700.0, 740.0, None, None]})
    comp = pd.DataFrame({"system_id": ["fuel:F_A", "fuel:F_B", "fuel:G_A", "fuel:G_B"],
                         "inchikey": "TOLUENE", "smiles": "Cc1ccccc1",
                         "fraction": [0.1, 0.3, 0.125, 0.325], "basis": "mass", "is_lumped": False,
                         "unexplained_fraction": 0.0})
    led = classify_fuel_rows(inv, comp, {"synthetic/ingest.py": _SNIPPET})
    fr = led[led["field"].str.startswith("fraction:")].set_index("fuel_id")
    assert (fr["provenance_class"] == "target_derived_composition").all(), fr.to_string()
    assert fr.loc["F_A", "evidence"].count("density_kg_m3") >= 1
