"""Tests for the v5 delta_ml_ysi gap-filler (ridge group-contribution YSI)."""
import json
from pathlib import Path

import numpy as np

from fuelblend.ml.ysi_model import canon, features, predict_ysi

ROOT = Path(__file__).resolve().parents[1]


def test_measured_is_reported_measured_and_unchanged():
    r = predict_ysi(["CCCCCCCCCCCC", "c1ccccc1"])
    assert list(r.source) == ["measured", "measured"]
    assert abs(r.ysi[0] - 71.7) < 1e-9 and abs(r.ysi[1] - 100.3) < 1e-9


def test_fixture_pinned():
    r = predict_ysi(["CCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCC"])
    assert abs(r.ysi[0] - 93.0496860815) < 1e-6 and r.source[0] == "ml_predicted"
    assert abs(r.ysi[1] - 109.7492087895) < 1e-6 and r.source[1] == "ml_predicted_extrapolated"


def test_n_alkanes_monotone_and_extrapolated_label():
    r = predict_ysi(["C" * n for n in range(13, 19)])
    assert np.all(np.diff(r.ysi.to_numpy()) > 0)
    assert r.source.iloc[-1] == "ml_predicted_extrapolated"
    assert np.all(r.lo <= r.ysi) and np.all(r.ysi <= r.hi)


def test_stereo_free_lookup_and_non_hc():
    assert canon("C1CC[C@H]2CCCC[C@@H]2C1") == canon("C1CCC2CCCCC2C1")
    r = predict_ysi(["CCO"])
    assert r.source[0] == "out_of_domain" and np.isnan(r.ysi[0])
    assert features("CCCC").shape == (16,)


def test_readings_keys():
    d = json.loads((ROOT / "results/v5/BASE/delta_ml_ysi/readings.json").read_text(encoding="utf-8"))
    for k in ("R-MLY-1", "R-MLY-1x", "R-MLY-2", "R-MLY-3", "R-MLY-4", "R-MLY-5", "adoption"):
        assert k in d
