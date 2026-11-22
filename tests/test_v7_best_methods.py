import numpy as np
import pandas as pd
import pytest
from fuelblend.api import best_methods as BM, predict_pure as PP


def test_registry_covers_ten_pure_properties_and_reads_back_committed_mae():
    assert len(BM.PURE) == 10
    for p in BM.PURE:
        assert BM.evidence_mae(p) > 0


def test_choice_is_lowest_committed_mae_or_declared_tied():
    rf = pd.read_csv(BM.RF_SRC[0]).set_index("target").mae_rf
    dm = pd.read_csv(BM.DM_SRC[0]).set_index("target").mae_long
    for p, e in BM.PURE.items():
        best = min(rf[p], dm[p])
        chosen = rf[p] if e.method == "rf_rdkit2d" else dm[p]
        assert chosen <= best + 1e-9 or "tied" in e.significance, p


def test_blend_registry_keeps_v5_for_designated_routes():
    for k in ("viscosity", "density", "flash_point", "freezing_point"):
        assert BM.best_method(k, "blend")["status"] == "adopted_v5"
    with pytest.raises(KeyError):
        BM.best_method("nope")


def test_predict_pure_rf_and_invalid_smiles(tmp_path, monkeypatch):
    monkeypatch.setattr(PP, "CACHE", tmp_path)
    r = PP.predict_pure(["CCCCCCCC", "not_a_smiles"], props=["dhfus_kJ_mol"])
    ok, bad = r["CCCCCCCC"]["dhfus_kJ_mol"], r["not_a_smiles"]["dhfus_kJ_mol"]
    assert ok["status"] == "ok" and np.isfinite(ok["value"]) and 0 < ok["value"] < 60
    assert bad["status"] == "invalid_smiles" and np.isnan(bad["value"])


def test_dmpnn_missing_weights_reports_model_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(PP, "CACHE", tmp_path)
    r = PP.predict_pure("CCCCCCCC", props=["dcn"])["CCCCCCCC"]["dcn"]
    assert r["status"] == "model_missing" and np.isnan(r["value"])
