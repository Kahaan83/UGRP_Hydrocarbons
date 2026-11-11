"""Tests for Phase B: Real-Fuel Validation Data (Tasks B1 - B6)."""
from pathlib import Path
import json
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REAL_FUELS_DIR = ROOT / "data" / "real_fuels"
SPLITS_PATH = ROOT / "splits" / "real_fuels_nested_v2.json"
MATRIX_PATH = ROOT / "results" / "v2" / "B6" / "real_fuels_property_matrix.csv"

from fuelblend.data.validate import validate_composition_df


def test_b1_candidate_inventory():
    cand_path = REAL_FUELS_DIR / "candidates_v2.csv"
    assert cand_path.exists(), f"Missing {cand_path}"
    df = pd.read_csv(cand_path)
    counts = df["fuel_family"].value_counts()
    assert counts.get("jet", 0) >= 20
    assert counts.get("gasoline", 0) >= 15
    assert counts.get("diesel", 0) >= 15


def test_b2_inventory_and_compositions_validation():
    inv_path = REAL_FUELS_DIR / "inventory.csv"
    comp_path = REAL_FUELS_DIR / "compositions.csv"
    assert inv_path.exists(), f"Missing {inv_path}"
    assert comp_path.exists(), f"Missing {comp_path}"

    df_inv = pd.read_csv(inv_path)
    df_comp = pd.read_csv(comp_path)

    assert len(df_inv) >= 90
    required_inv_cols = [
        "fuel_id", "fuel_name", "source", "fuel_family", "composition_type",
        "number_of_resolved_components", "unresolved_fraction", "unresolved_handling"
    ]
    for col in required_inv_cols:
        assert col in df_inv.columns

    # Verify every fuel has resolved component count >= 1 and unresolved fraction in [0, 1]
    assert (df_inv["number_of_resolved_components"] >= 1).all()
    assert (df_inv["unresolved_fraction"] >= 0.0).all()
    assert (df_inv["unresolved_fraction"] <= 1.0).all()

    # Verify composition validator passes with 0 errors
    val_report = validate_composition_df(df_comp, strict=False)
    assert val_report["valid"] is True
    assert len(val_report["errors"]) == 0


def test_b3_strict_segregation_measured_vs_reconstructed():
    inv_meas_path = REAL_FUELS_DIR / "inventory_measured.csv"
    inv_recon_path = REAL_FUELS_DIR / "inventory_reconstructed.csv"
    comp_meas_path = REAL_FUELS_DIR / "compositions_measured.csv"
    comp_recon_path = REAL_FUELS_DIR / "compositions_reconstructed.csv"

    assert inv_meas_path.exists()
    assert inv_recon_path.exists()
    assert comp_meas_path.exists()
    assert comp_recon_path.exists()

    inv_meas = pd.read_csv(inv_meas_path)
    inv_recon = pd.read_csv(inv_recon_path)

    # 1. No overlap between measured and reconstructed inventories
    meas_ids = set(inv_meas["fuel_id"])
    recon_ids = set(inv_recon["fuel_id"])
    assert len(meas_ids.intersection(recon_ids)) == 0

    # 2. Reconstructed inventory must be 100% reconstructed_composition (CRC 647)
    assert (inv_recon["composition_type"] == "reconstructed_composition").all()
    assert len(inv_recon) == 57

    # 3. Measured inventory must be 100% measured_dha_gcxgc
    assert (inv_meas["composition_type"] == "measured_dha_gcxgc").all()
    assert len(inv_meas) >= 30

    # 4. Compositions must align strictly with inventories
    comp_meas = pd.read_csv(comp_meas_path)
    comp_recon = pd.read_csv(comp_recon_path)
    meas_sys = set("fuel:" + inv_meas["fuel_id"])
    recon_sys = set("fuel:" + inv_recon["fuel_id"])
    assert set(comp_meas["system_id"]).issubset(meas_sys)
    assert set(comp_recon["system_id"]).issubset(recon_sys)
    assert len(set(comp_meas["system_id"]).intersection(set(comp_recon["system_id"]))) == 0


def test_b4_blend_links():
    links_path = REAL_FUELS_DIR / "blend_links.csv"
    assert links_path.exists()
    df_links = pd.read_csv(links_path)
    assert len(df_links) >= 10
    required_cols = ["blend_id", "parent_fuel_1", "target_fuel", "measured_properties"]
    for col in required_cols:
        assert col in df_links.columns


def test_b5_nested_splits_no_overlap():
    assert SPLITS_PATH.exists()
    with open(SPLITS_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    assert "outer_folds" in data
    assert len(data["outer_folds"]) == 5

    all_test_fuels = []
    for fold in data["outer_folds"]:
        test_fuels = fold["outer_test_fuels"]
        train_fuels = fold["outer_train_fuels"]

        # No overlap within outer fold
        overlap = set(test_fuels).intersection(set(train_fuels))
        assert len(overlap) == 0, f"Overlap in outer fold {fold['outer_fold']}: {overlap}"
        all_test_fuels.extend(test_fuels)

        # Check inner folds
        assert len(fold["inner_nested_folds"]) == 4
        for in_fold in fold["inner_nested_folds"]:
            in_tr = in_fold["inner_train_fuels"]
            in_val = in_fold["inner_val_fuels"]
            assert len(set(in_tr).intersection(set(in_val))) == 0

    # Outer test folds must partition all real fuels exactly once
    assert len(all_test_fuels) == len(set(all_test_fuels))


