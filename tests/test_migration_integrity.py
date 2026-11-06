"""Integrity verification test for migrated core measurements (T1.5)."""
import hashlib
from pathlib import Path
import pandas as pd

from fuelblend.data.validate import validate_measurement_df

ROOT = Path(__file__).resolve().parents[1]
MIGRATED_PATH = ROOT / "data" / "processed" / "measurements_migrated_core.csv"
MODEL_TABLE_PATH = ROOT / "data" / "processed" / "model_table.csv"


def test_migrated_data_integrity():
    assert MIGRATED_PATH.exists(), f"Missing {MIGRATED_PATH}"
    df_migrated = pd.read_csv(MIGRATED_PATH, low_memory=False)

    # 1. Validation test
    res = validate_measurement_df(df_migrated, strict=True)
    assert res["valid"] is True

    # 2. Match with model_table core rows
    df_src = pd.read_csv(MODEL_TABLE_PATH, low_memory=False)
    target_props = ["viscosity_kin_cSt", "viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"]
    core_src = df_src[df_src["property"].isin(target_props)].reset_index(drop=True)

    assert len(df_migrated) == len(core_src)

    hash_migrated = hashlib.sha256(df_migrated["value"].round(6).astype(str).str.cat().encode("utf-8")).hexdigest()
    hash_src = hashlib.sha256(core_src["value"].round(6).astype(str).str.cat().encode("utf-8")).hexdigest()

    assert hash_migrated == hash_src
