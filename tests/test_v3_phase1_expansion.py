from pathlib import Path
import json
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

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
