import hashlib
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("bu", ROOT / "scripts/v9/build_unified_v9.py")
bu = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bu)
D9 = ROOT / "data/unified_v9/measurements_dedup.csv"
KEY = ["property", "system_id", "x_key", "T_key", "P_key"]


def test_canonical_id_matches_ajftd_label():
    assert bu.canonical_id("boehm2022", "A-2 (POSF 10325)") == "fuel:2"
    assert bu.canonical_id("boehm2022", "A-3 (POSF10289)") == "fuel:39"
    assert bu.canonical_id("boehm2022", "C-3 (POSF 12959)").startswith("fuel:boehm2022:")


def test_v5_rows_unchanged_unless_merged():
    a = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", low_memory=False)
    b = pd.read_csv(D9, low_memory=False)
    for df in (a, b):
        for c in ("x_key", "T_key", "P_key"):
            df[c] = df[c].astype(str)
    j = a.merge(b[~b.v9_new], on=KEY, suffixes=("", "_9"))
    same = j[~j.v9_merged.astype(bool)]
    assert len(same) >= len(a) - 400
    assert (same.value.sub(same.value_9).abs() < 1e-9).all()


def test_v9_new_groups_have_no_v5_row():
    a = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=KEY, low_memory=False).astype(str)
    b = pd.read_csv(D9, low_memory=False)
    new = b[b.v9_new][KEY].astype(str)
    assert len(new.merge(a, on=KEY)) == 0


def test_new_freezing_rows_land_on_known_fuel_ids():
    b = pd.read_csv(D9, low_memory=False)
    rf = b[b.v9_new & (b.kind == "real_fuel")]
    assert (rf.property.isin(["viscosity_kin_cSt", "density_kg_m3"])).all()
    assert rf.system_id.isin(["fuel:2", "fuel:39", "fuel:32", "fuel:24570", "fuel:143", "fuel:304"]).sum() > 0


def test_unified_v5_untouched():
    h = hashlib.sha256((ROOT / "data/unified/measurements_dedup.csv").read_bytes()).hexdigest()
    assert len(h) == 64
    assert not (ROOT / "data/unified/fuel_identity_map.csv").exists()
