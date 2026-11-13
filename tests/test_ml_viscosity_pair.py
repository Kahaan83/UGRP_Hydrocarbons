import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("lightgbm")
from fuelblend.ml import viscosity_pair as vp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
S = ["CCCCCCCCCC", "c1ccccc1C", "C1CCCCC1C"]


def test_pure_component_zero():
    assert vp.ln_correction(S[:1], [1.0], 298.15) == 0.0
    assert vp.ln_correction(S, [1.0, 0.0, 0.0], 298.15) == 0.0


def test_order_invariant_and_symmetric():
    a = vp.ln_correction(S, [0.3, 0.4, 0.3], 298.15)
    b = vp.ln_correction(S[::-1], [0.3, 0.4, 0.3][::-1], 298.15)
    assert abs(a - b) < 1e-12
    G = vp.interaction_matrix(S, 298.15)
    assert np.allclose(G, G.T) and np.allclose(np.diag(G), 0)


def test_fixture_prediction_pinned():
    assert abs(vp.ln_correction(S, [0.3, 0.4, 0.3], 298.15) - (-0.11807490658002913)) < 1e-9
    assert abs(vp.ln_correction(S[:2], [0.5, 0.5], 298.15) - (-0.09307318764893947)) < 1e-9


def test_readings_keys_and_meta():
    r = json.loads((ROOT / "results/v5/BASE/delta_ml_viscosity/readings.json").read_text())
    for k in ("R-MLV-1", "R-MLV-2", "R-MLV-3", "R-MLV-4", "R-MLV-6", "adoption"):
        assert k in r
    assert (ROOT / "models/v5/ml_viscosity_pair/meta.json").exists()
    assert vp.domain_status(S, 298.15)["in_domain"]
