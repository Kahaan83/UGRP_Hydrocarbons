from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
ANCHORED_REG_PATH = ROOT / "data" / "processed" / "v2" / "compounds_anchored_vp.csv"
FLASH_EVAL_PATH = ROOT / "results" / "v2" / "improvement" / "flash_anchored_loso.csv"


def test_anchored_files_exist():
    assert ANCHORED_REG_PATH.exists()
    assert FLASH_EVAL_PATH.exists()


def test_isododecane_ambient_vapor_pressure():
    from fuelblend.thermo.vapor_pressure_anchoring import get_anchored_vapor_pressure_kpa
    # Isododecane (2,2,4,6,6-pentamethylheptane, CAS 13475-82-6)
    # Old WebBook was 0.0029 kPa (erroneous 91x underestimation)
    # Anchored vapor pressure must be ~0.2-0.3 kPa at 25 C (298.15 K)
    p_kpa = get_anchored_vapor_pressure_kpa("VKPSKYDESGTTFR-UHFFFAOYSA-N", 298.15)
    assert 0.15 <= p_kpa <= 0.40


def test_flash_anchored_loso_error_reduction():
    df = pd.read_csv(FLASH_EVAL_PATH)
    assert len(df) == 126

    old_max = df["old_base_error"].max()
    new_max = df["new_anchored_error"].max()

    # Verify maximum over-prediction dropped from > 40 C down to <= 11 C
    assert old_max > 40.0
    assert new_max <= 11.0

    # Verify dangerous over-prediction rate > 5 C dropped below 2%
    over_5c = (df["new_anchored_error"] > 5.0).mean()
    assert over_5c < 0.02
