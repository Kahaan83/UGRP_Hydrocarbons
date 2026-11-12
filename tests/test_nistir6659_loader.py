"""NISTIR 6659 density loader (fixture = 10 real rows of the extracted table)."""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.loaders_nistir6659 import COLUMNS, load_nistir6659  # noqa: E402

FIX = ROOT / "tests/data/nistir6659_10rows.csv"


def test_loader_schema_and_units():
    d = load_nistir6659(FIX)
    assert list(d.columns) == COLUMNS and len(d) == 10
    assert (d.posf == "4658").all() and (d.property == "density_kg_m3").all()
    assert d.value.between(600, 900).all() and d.T_K.between(250, 500).all() and (d.P_kPa > 100).all()   # kg/m3, K, kPa


def test_loader_rejects_mpa_as_kpa(tmp_path):
    raw = pd.read_csv(FIX)
    raw["P_kPa"] = raw.P_MPa
    p = tmp_path / "bad.csv"
    raw.to_csv(p, index=False)
    with pytest.raises(ValueError):
        load_nistir6659(p)
