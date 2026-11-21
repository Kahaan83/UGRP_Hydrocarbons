"""v6 Track A checks: lock integrity, group disjointness, tile-mean tmap behaviour, v5 reproduction gate."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/v6/freezing"


def test_lock_matches_dataset():
    lock = json.load(open(OUT / "holdout_lock.json"))
    d = pd.read_csv(OUT / "dataset.csv")
    assert d.group.is_unique
    assert set(d[d.role == "holdout"].entity) == set(lock["holdout_entities"])
    assert (d.role == "holdout").sum() == lock["n_holdout"] == 12 and (d.role == "dev").sum() == lock["n_dev"] == 52


def test_holdout_has_no_new_dev_fuel_overlap():
    d = pd.read_csv(OUT / "dataset.csv", dtype={"fuel_id": str})
    assert "20" in set(d[d.role == "holdout"].fuel_id)
    assert not {"32"} & set(d.fuel_id)        # same fuel as FuelLib posf10264, not added twice


def test_repro_gate_passed():
    g = json.load(open(OUT / "repro_gate.json"))
    assert g["passed"] and g["n_headline"] == 50
    assert abs(g["mae"] - 8.4245951) < 1e-4 and abs(g["bias"] - 5.8778236) < 1e-4


def test_tile_mean_tmap_uses_mean_and_skips_small_tiles():
    from fuelblend.data.tile_isomers import tile_mean_tmap

    class Reps:
        ik = ["A", "B", "C", "D"]
        by = {("iso_paraffin", 10): [0, 1], ("cyclo_mono", 9): [2]}

    class Reg:
        ik, smiles = ["A", "B", "C", "D"], ["CC(C)CCCCCCC", "CCC(C)CCCCCC", "C1CCCCC1CCC", "CCCCCCCCCC"]

    tmap = {"A": dict(tm_K=200.0, dh_fus_J_mol=20000.0, source="x"), "B": dict(tm_K=220.0, dh_fus_J_mol=22000.0, source="x"),
            "C": dict(tm_K=210.0, dh_fus_J_mol=21000.0, source="x")}
    out = tile_mean_tmap(Reg, Reps, tmap)
    assert out["A"]["tm_K"] == 210.0 and out["A"]["dh_fus_J_mol"] == 21000.0 and out["A"]["n_isomers"] == 2
    assert out["B"]["tm_K"] == 220.0          # only the representative (first member) is replaced
    assert out["C"]["tm_K"] == 210.0          # single-isomer tile unchanged
    assert tmap["A"]["tm_K"] == 200.0         # input not mutated
