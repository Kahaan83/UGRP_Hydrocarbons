"""v5 delta_freeze_margin: verified Tm / dHfus map, liquidus scoring readings, yaml entry, API and design margin mechanics."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from rdkit import Chem
from rdkit.Chem import inchi

from fuelblend.api import predict_route as pr
from fuelblend.design import engine_v4 as E
from fuelblend.design import library_route as LR
from fuelblend.eval.composition_route import (COUTINHO_SOURCE, VERIFIED_SOURCE, build_tmap_designated, build_tmap_verified, ideal_sle_C)

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "v5" / "BASE" / "delta_freeze_margin"
KEY = "freezing_point_C.lab_liquidus.svl_verified"
UNDECANE, DODECANE, BCH, TOLUENE = "CCCCCCCCCCC", "CCCCCCCCCCCC", "CCCCC1CCCCC1", "CCCCCCc1ccccc1"   # TOLUENE now hexylbenzene: toluene gained verified constants 2026-11-14; hexylbenzene has no dHfus source


def _ik(s):
    return inchi.MolToInchiKey(Chem.MolFromSmiles(s))


def test_verified_map_content():
    t = build_tmap_verified()
    assert len(t) == 44                                           # 46 csv rows (39 + 7 upgraded 2026-11-14); stereo-unspecified decalin and 1-hexyne have no constants
    u = t[_ik(UNDECANE)]
    assert u["tm_K"] == pytest.approx(247.59) and u["dh_fus_J_mol"] == pytest.approx(29037.0)    # total incl. solid-solid for the odd alkane
    assert t[_ik("c1ccc2c(c1)Cc1ccccc1-2")]["tm_K"] == pytest.approx(387.94)       # fluorene: not benzene's 278.7 K (Coutinho defect)
    assert t[_ik("Cc1ccc2ccccc2c1")]["tm_K"] == pytest.approx(307.73)              # 2-methylnaphthalene: not naphthalene's 353.4 K


def test_designated_map_prefers_verified_and_tags_source():
    base = {"A": {"tm_K": 200.0, "dh_fus_J_mol": 1.0e4}, _ik(UNDECANE): {"tm_K": 247.6, "dh_fus_J_mol": 33400.0}}
    d = build_tmap_designated(base, build_tmap_verified())
    assert d["A"]["source"] == COUTINHO_SOURCE
    assert d[_ik(UNDECANE)]["source"] == VERIFIED_SOURCE and d[_ik(UNDECANE)]["dh_fus_J_mol"] == pytest.approx(29037.0)
    assert d[_ik(BCH)]["source"] == VERIFIED_SOURCE                                  # butylcyclohexane: absent from the Coutinho table


def test_pure_compound_svl_is_its_melting_point():
    t = build_tmap_verified()
    ik = _ik(UNDECANE)
    assert ideal_sle_C([1.0], [ik], t) == pytest.approx(247.59 - 273.15)


def test_readings_are_consistent_with_summary_and_rule():
    r = json.loads((D / "readings.json").read_text(encoding="utf-8"))
    s = pd.read_csv(D / "summary.csv").query("scope == 'S_F' and model == 'svl_verified'").iloc[0]
    assert r["R-FM-1"]["n_systems"] == int(s.n_systems) and r["R-FM-1"]["passed"] == (int(s.n_systems) >= 9)
    assert r["R-FM-2"]["passed"] == bool(s.mae < s.mae_const)
    assert r["R-FM-3"]["margin_C"] == pytest.approx(float(s.margin_alpha10)) and r["R-FM-3"]["passed"] == math.isfinite(float(s.margin_alpha10))
    assert r["rule_write_margin"] == (r["R-FM-1"]["passed"] and r["R-FM-2"]["passed"] and r["R-FM-3"]["passed"])


def test_scope_rows_obey_declared_scope():
    rows = pd.read_csv(D / "rows.csv")
    sf = rows[rows.in_SF]
    assert sf.all_hc.all() and sf.all_c8_c16.all() and sf.all_verified.all() and (sf["T_K"] <= 273.15).all()
    assert sf.pred_verified.notna().all()


def test_margin_recomputes_from_rows():
    from fuelblend.uq.conformal import one_sided_margin
    rows = pd.read_csv(D / "rows.csv")
    sf = rows[rows.in_SF]
    m = one_sided_margin((sf.pred_verified - sf.measured).to_numpy(), sf.group.to_numpy(), alpha=0.10, direction="under")
    assert m == pytest.approx(json.loads((D / "readings.json").read_text(encoding="utf-8"))["R-FM-3"]["margin_C"])


def test_yaml_entry_written_iff_rule_passed():
    doc = yaml.safe_load((ROOT / "configs" / "uncertainty_v5.yaml").read_text(encoding="utf-8"))["margins"]
    rd = json.loads((D / "readings.json").read_text(encoding="utf-8"))
    assert (KEY in doc) == rd["rule_write_margin"]
    assert doc["flash_point_C.lab_mixtures.ideal"]["applied_margin"] == 2.206932808388615
    assert doc["smoke_point_mm.lab_mixtures.bulk_regression"]["applied_margin"] == 7.8009965628100595
    if KEY in doc:
        e = doc[KEY]
        assert e["direction"] == "under" and e["unit"] == "C" and e["applied_margin"] == pytest.approx(rd["R-FM-3"]["margin_C"])
        assert e["n_groups"] == rd["R-FM-1"]["n_systems"]


def test_api_margin_known_all_verified():
    m = yaml.safe_load((ROOT / "configs" / "uncertainty_v5.yaml").read_text(encoding="utf-8"))["margins"][KEY]["applied_margin"]
    f = pr.predict_fuel({UNDECANE: 0.4, DODECANE: 0.2, BCH: 0.4}, T_K=288.15, fuel_type="jet")["freezing_point_C"]
    assert f["margin"] == pytest.approx(m) and f["margin_direction"] == "under" and f["in_domain"] is True and f["uncertified"] is False
    assert f["hi"] == pytest.approx(f["value"] + m) and f["safe_bound"] == f["hi"] and f["lo"] == pytest.approx(f["value"] - m)
    assert f["variants"]["coutinho_constants"]["value"] is None or math.isfinite(f["variants"]["coutinho_constants"]["value"])


def test_api_margin_withheld_for_reconstructed_and_for_unverified_component():
    rec = pr.predict_fuel({UNDECANE: 0.5, DODECANE: 0.5}, T_K=288.15, fuel_type="jet", composition_kind="reconstructed")["freezing_point_C"]
    assert rec["uncertified"] is True and rec["in_domain"] is False and math.isnan(rec["hi"]) and rec["value"] is not None
    tol = pr.predict_fuel({UNDECANE: 0.5, DODECANE: 0.3, TOLUENE: 0.2}, T_K=288.15, fuel_type="jet")["freezing_point_C"]
    assert tol["uncertified"] is True and math.isnan(tol["hi"]) and tol["value"] is not None       # hexylbenzene has no verified constant


def test_api_designated_value_uses_verified_constants():
    f = pr.predict_fuel({UNDECANE: 0.5, DODECANE: 0.5}, T_K=288.15, fuel_type="jet")["freezing_point_C"]
    t = build_tmap_verified()
    from rdkit.Chem import Descriptors
    x = np.array([0.5 / Descriptors.MolWt(Chem.MolFromSmiles(UNDECANE)), 0.5 / Descriptors.MolWt(Chem.MolFromSmiles(DODECANE))])
    x = x / x.sum()
    assert f["value"] == pytest.approx(ideal_sle_C(x, [_ik(UNDECANE), _ik(DODECANE)], t), abs=1e-6)


def test_design_model_carries_margin_key_and_flags_coutinho_members():
    names = ["n-undecane", "n-dodecane", "butylcyclohexane"]
    models = LR.build_route_models(names)
    fz = models["props"]["freezing"]
    assert fz.margin_keys == (KEY,) and fz.margin_missing_reason == "" and fz.certified_kinds == {"jet"}
    margins = E.load_margins(ROOT / "configs" / "uncertainty_v5.yaml")
    c = E.evaluate_constraint("freezing_point_C", {"max": 0.0}, models["pool"], np.array([0.4, 0.2, 0.4]), models, margins, kind="jet")
    assert c.applied_margin == pytest.approx(margins[KEY]["applied_margin"]) and c.status in ("pass", "fail")
    assert c.slack == pytest.approx((0.0 - (c.prediction + c.applied_margin)) / 1.0) or c.status == "fail"
    # a pool member on a Coutinho entry (n-nonane has no verified constant) keeps freezing uncertified
    m2 = LR.build_route_models(names + ["C" * 9])
    c2 = E.evaluate_constraint("freezing_point_C", {"max": 0.0}, m2["pool"], np.array([0.3, 0.2, 0.3, 0.2]), m2, margins, kind="jet")
    assert c2.status == "uncertified" and "non-verified constant" in c2.reason
