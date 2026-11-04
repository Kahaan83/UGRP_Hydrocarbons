"""Unit tests for Phase 3 Real-Fuel Dataset (T3.1 - T3.5)."""
import json
from pathlib import Path
import pandas as pd

from fuelblend.data.validate import validate_composition_df

ROOT = Path(__file__).resolve().parents[1]
INV_PATH = ROOT / "data" / "real_fuels" / "inventory.csv"
COMP_PATH = ROOT / "data" / "real_fuels" / "compositions.csv"
LINKS_PATH = ROOT / "data" / "real_fuels" / "blend_links.csv"
SPLIT_PATH = ROOT / "splits" / "real_fuels_holdout.json"


def test_real_fuels_inventory():
    assert INV_PATH.exists(), f"Missing {INV_PATH}"
    df = pd.read_csv(INV_PATH)
    assert len(df) >= 60

    # T3.4 strict segregation check:
    types = set(df["composition_type"].unique())
    assert types == {"measured_dha_gcxgc", "reconstructed_composition"}

    crc_fuels = df[df["source"] == "crc647"]
    assert (crc_fuels["composition_type"] == "reconstructed_composition").all()

    rp1 = df[df["fuel_id"] == "NIST_RP1"].iloc[0]
    assert rp1["composition_type"] == "measured_dha_gcxgc"
    assert rp1["number_of_resolved_components"] == 33


def test_real_fuels_compositions_validation():
    assert COMP_PATH.exists(), f"Missing {COMP_PATH}"
    df = pd.read_csv(COMP_PATH)
    res = validate_composition_df(df, strict=True)
    assert res["valid"] is True


def test_blend_links_exists():
    assert LINKS_PATH.exists(), f"Missing {LINKS_PATH}"
    df = pd.read_csv(LINKS_PATH)
    assert len(df) >= 3
    assert "target_fuel" in df.columns
    assert "parent_fuel_1" in df.columns


def test_real_fuels_holdout_split():
    assert SPLIT_PATH.exists(), f"Missing {SPLIT_PATH}"
    data = json.loads(SPLIT_PATH.read_text(encoding="utf-8"))
    assert "holdout_real_fuel_ids" in data
    assert len(data["holdout_real_fuel_ids"]) >= 60
    assert "reconstructed_composition_fuels" in data
