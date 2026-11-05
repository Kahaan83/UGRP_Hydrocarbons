from pathlib import Path
import json
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_phase1_deliverables_exist():
    assert (ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv").exists()
    assert (ROOT / "data" / "real_fuels" / "v3" / "compositions_v3.csv").exists()
    assert (ROOT / "results" / "v3" / "S4" / "coverage_matrix.csv").exists()
    assert (ROOT / "splits" / "v3" / "real_fuels_v3_splits.json").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "PHASE_1_SUMMARY.md").exists()

def test_s1_gasoline_count():
    inv = pd.read_csv(ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv")
    gasolines = inv[inv["fuel_family"] == "gasoline"]
    assert len(gasolines) >= 30, f"Expected >= 30 gasolines, got {len(gasolines)}"
    assert "ethanol_vol_pct" in inv.columns or "oxygenate_pct" in inv.columns

def test_s2_diesel_count():
    inv = pd.read_csv(ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv")
    diesels = inv[inv["fuel_family"] == "diesel"]
    assert len(diesels) >= 15, f"Expected >= 15 diesels, got {len(diesels)}"
    # Verify cetane labeling exists
    assert "cetane_test_type" in inv.columns or any("cetane" in c.lower() for c in inv.columns)

def test_s4_coverage_matrix():
    cov = pd.read_csv(ROOT / "results" / "v3" / "S4" / "coverage_matrix.csv")
    assert "fuels_with_composition_and_measurement" in cov.columns
    assert "status" in cov.columns
    # Check that flags are present if count < 10
    below_10 = cov[cov["fuels_with_composition_and_measurement"] < 10]
    if len(below_10) > 0:
        assert (below_10["status"] == "FLAG_BELOW_10").all()

def test_s5_split_no_leakage():
    splits_path = ROOT / "splits" / "v3" / "real_fuels_v3_splits.json"
    with open(splits_path, "r") as f:
        data = json.load(f)
    assert "outer_folds" in data or "folds" in data or "held_out_fuels" in data
    # Verify 0% overlap between train and test/val in folds if folds exist
    folds = data.get("outer_folds") or data.get("folds", [])
    for fold in folds:
        train_set = set(fold.get("train", []))
        test_set = set(fold.get("test", []) + fold.get("val", []))
        assert len(train_set.intersection(test_set)) == 0, "Train and held-out overlap detected!"

def test_s6_oxygenates_registered():
    comp = pd.read_csv(ROOT / "data" / "unified" / "compounds.csv")
    ethanol = comp[comp["smiles"].str.upper() == "CCO"]
    assert len(ethanol) >= 1, "Ethanol missing from compounds.csv"
    row = ethanol.iloc[0]
    assert pd.notna(row["mw"])
    assert pd.notna(row["antoine_A"])
    assert pd.notna(row["lfl_vol_pct"])
