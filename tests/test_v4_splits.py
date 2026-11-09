"""v4 Task 7: real-fuel group splits (5-fold GroupKFold + locked test set)."""
import importlib.util
import json
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("v4_splits", ROOT / "scripts" / "v4" / "build_splits_v4.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

V4 = ROOT / "data" / "real_fuels" / "v4"
pytestmark = pytest.mark.skipif(not (V4 / "inventory_v4.csv").exists(), reason="v4 corpus absent")


@pytest.fixture(scope="module")
def built():
    return S.build()


def _roles(sp):
    return sp["fuel_role"]  # fuel_id -> "locked" | "fold0".."fold4"


def test_no_fuel_in_two_folds(built):
    sp = built["splits"]
    seen = {}
    for k, fuels in sp["folds"].items():
        for f in fuels:
            assert f not in seen, f"{f} in {seen[f]} and {k}"
            seen[f] = k
    assert len(sp["folds"]) == 5
    # a GROUP never straddles folds either
    for g, members in sp["groups"].items():
        roles = {_roles(sp)[m] for m in members}
        assert len(roles) == 1, (g, roles)


def test_locked_set_disjoint(built):
    sp = built["splits"]
    locked = set(sp["locked"])
    assert locked
    dev = {f for fuels in sp["folds"].values() for f in fuels}
    assert not locked & dev
    assert locked | dev == set(_roles(sp))
    assert sp["seed"] == 20261030
    frac = len(locked) / len(_roles(sp))
    assert 0.12 <= frac <= 0.28, frac
    # locked group list is chosen from groups, folds built over the remainder only
    assert not set(sp["locked_groups"]) & {g for g, m in sp["groups"].items() if _roles(sp)[m[0]] != "locked"}


def test_union_covers_r1_sources_and_uses_included_jet(built):
    from fuelblend.data.loaders_jet import included_jet_fuel_ids
    inv = set(pd.read_csv(V4 / "inventory_v4.csv").fuel_id)
    jet = pd.read_csv(V4 / "jet_class_compositions.csv")
    gas = set(pd.read_csv(V4 / "gasoline_octane_blends.csv").blend_id)
    included = set(included_jet_fuel_ids(jet))
    ids = set(_roles(built["splits"]))
    assert ids == inv | included | gas
    assert "JP5_PotterSimmons" not in ids and "JP8_PotterSimmons" not in ids
    assert "JetA_API14avg" in ids
    assert not inv & included


def test_no_v3_v4_duplicate_fuel(built):
    ids = list(_roles(built["splits"]))
    assert not [i for i in ids if i.startswith("POSF_")]  # v3 spelling absent
    norm = {}
    for i in ids:
        k = re.sub(r"[-_\s]", "", i).upper()
        assert k not in norm, (i, norm[k])
        norm[k] = i
    v3 = pd.read_csv(ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv").fuel_id.unique()
    for v in v3:
        if str(v).startswith("POSF_"):
            assert v not in ids


def test_surrogates_stay_with_real_rp1(built):
    g = built["splits"]["fuel_group"]
    assert len({g[f] for f in ["NIST_RP1", "RP1_Binary_2comp", "RP1_Edwards_3comp", "RP1_NIST_14comp"]}) == 1


def test_bron_blendstock_family_not_split(built):
    sp = built["splits"]
    gas = pd.read_csv(V4 / "gasoline_octane_blends.csv")
    bron = gas[gas.blend_id.str.startswith("BRON:")]
    for cas, sub in bron.groupby(bron.blend_id.str.split(":").str[1]):
        assert {sp["fuel_group"][b] for b in sub.blend_id} == {f"BRON:{cas}"}, cas
        assert len({_roles(sp)[b] for b in sub.blend_id}) == 1, cas
    for base, sub in gas[~gas.blend_id.str.startswith("BRON:")].groupby("base_fuel"):
        assert len({_roles(sp)[b] for b in sub.blend_id}) == 1, base


def test_locked_fraction_recorded_and_bounded(built):
    c = built["splits"]["counts"]
    assert 0.15 <= c["locked_groups"] / c["groups"] <= 0.25
    assert c["locked_fraction_of_fuels"] == pytest.approx(c["locked_fuels"] / c["fuels"], abs=1e-3)
    assert c["fold_size_max_over_median"] <= 2.0


def test_locked_has_ron_and_mon_targets(built):
    cov, sp = built["coverage"], built["splits"]
    lk = cov[cov.fuel_id.isin(sp["locked"]) & (cov.role == "target") & cov.present]
    assert lk[lk.property == "ron"].fuel_id.nunique() >= 1
    assert lk[lk.property == "mon"].fuel_id.nunique() >= 1
    assert sp["counts"]["target_fuels_by_split_and_property"]["locked"].get("ron", 0) >= 1


def test_different_seed_changes_locked_and_stored_json_matches(built):
    other = S.build(seed=1)["splits"]
    assert other["locked"] != built["splits"]["locked"]
    stored = S.OUT_SPLITS
    if stored.exists():
        js = json.loads(stored.read_text(encoding="utf-8"))
        assert js["locked"] == built["splits"]["locked"]
        assert js["folds"] == built["splits"]["folds"]


def test_headline_note_states_dev_fold(built):
    sp = built["splits"]
    assert "NOT held out" in sp["headline_note"] or "locked set" in sp["headline_note"]
    assert "POSF-10325" in sp["headline_note"]


def test_headline_flags_carried(built):
    sp, cov = built["splits"], built["coverage"]
    assert sp["headline_set"] == ["POSF-10325"]
    assert sp["headline_statistics_allowed"] is False
    row = cov[cov.fuel_id == "POSF-10325"]
    assert row.headline_ok.all() and row.cn_split_class_independent.all()
    assert not cov[cov.fuel_id != "POSF-10325"].headline_ok.any()
    assert {"fuel_id", "property", "present", "source", "provenance", "headline_ok"} <= set(cov.columns)
    assert set(cov.fuel_id) == set(_roles(sp))


def test_deterministic(built):
    again = S.build()
    assert json.dumps(again["splits"], sort_keys=True) == json.dumps(built["splits"], sort_keys=True)
    pd.testing.assert_frame_equal(again["coverage"], built["coverage"])


def test_diesel_has_no_qualifying_fuels_and_sd02_recorded():
    assert S.diesel_qualifying_fuels() == []
    txt = (ROOT / "docs" / "scope_decisions.md").read_text(encoding="utf-8")
    assert "SD-02" in txt and "Task 8" in txt
