import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("lightgbm")
from fuelblend.ml import density_pair as dp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
S = ["CCCCCCCCCC", "c1ccccc1C", "C1CCCCC1C"]
RHO = [726.0, 862.0, 765.0]


def test_pure_component_zero_correction():
    assert dp.excess_volume(S[:1], [1.0], 298.15, rho_i=RHO[:1]) == 0.0
    assert dp.excess_volume(S, [1.0, 0.0, 0.0], 298.15, rho_i=RHO) == 0.0


def test_order_invariant():
    a = dp.excess_volume(S, [0.3, 0.4, 0.3], 298.15, rho_i=RHO)
    b = dp.excess_volume(S[::-1], [0.3, 0.4, 0.3][::-1], 298.15, rho_i=RHO[::-1])
    assert abs(a - b) < 1e-12


def test_density_between_ideal_and_physical():
    r = dp.density_ml(S, [0.3, 0.4, 0.3], 298.15, rho_i=RHO)
    ideal = 1.0 / sum(w / rho for w, rho in zip(
        np.array([0.3, 0.4, 0.3]) * np.array([142.28, 92.14, 98.19]) / np.dot([0.3, 0.4, 0.3], [142.28, 92.14, 98.19]), RHO))
    assert abs(r["rho_ideal"] - ideal) < 0.1
    assert abs(r["rho_ml"] - r["rho_ideal"]) < 15.0


def test_readings_keys_and_meta():
    r = json.loads((ROOT / "results/v5/BASE/delta_ml_density/readings.json").read_text())
    for k in ("R-MLD-1", "R-MLD-2", "R-MLD-3", "R-MLD-4", "R-MLD-5", "adoption"):
        assert k in r
    assert (ROOT / "models/v5/ml_density_pair/meta.json").exists()
    assert dp.domain_status(S)["in_domain"]
