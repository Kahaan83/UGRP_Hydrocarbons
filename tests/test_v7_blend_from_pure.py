import math
import numpy as np
from fuelblend.api.predict_route import predict_fuel
from fuelblend.api import blend_from_pure as B

BLEND = [("CCCCCCCCCCCC", .5), ("Cc1ccccc1", .5)]


def test_existing_keys_unchanged_when_exploratory_toggled():
    a, b = predict_fuel(BLEND, T_K=298.15), predict_fuel(BLEND, T_K=298.15, exploratory=False)
    assert set(a) - set(b) == {"dcn_blend_exploratory", "ron_blend_exploratory", "mon_blend_exploratory"}
    for k in b:
        if k not in ("components", "warnings"):
            assert str(a[k]) == str(b[k]), k


def test_exploratory_entries_are_uncertified_and_mon_is_no_skill():
    r = predict_fuel(BLEND)
    for k in ("dcn_blend_exploratory", "ron_blend_exploratory"):
        d = r[k]
        assert d["uncertified"] and not d["in_domain"] and math.isnan(d["lo"]) and d["margin"] is None
        assert d["status"] in ("exploratory", "model_missing")
    assert r["mon_blend_exploratory"]["status"] == "no_skill_shown" and r["mon_blend_exploratory"]["value"] is None


def test_value_is_linear_mole_fraction_of_components():
    d = predict_fuel(BLEND)["dcn_blend_exploratory"]
    if d["status"] != "exploratory":
        return                                           # D-MPNN weights absent on this machine
    c = d["components"]
    assert abs(sum(x["x_mole"] for x in c) - 1) < 1e-9
    assert abs(d["value"] - sum(x["x_mole"] * x["pure_value"] for x in c)) < 1e-9


def test_evidence_matches_committed_summary():
    ev = B._evidence()
    assert ev["dcn_leenhouts"]["n_systems"] == 31 and bool(ev["dcn_leenhouts"]["usable_exploratory"])
    assert not bool(ev["mon_leenhouts"]["usable_exploratory"])
