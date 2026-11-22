import numpy as np
import pytest

pytest.importorskip("transformers")
from pathlib import Path
from fuelblend.thermo import gamma_hanna as G

pytestmark = pytest.mark.skipif(not (G._paths() / "Model/HANNA_Val.pt").exists(), reason="HANNA repo not cloned")


def test_hanna_pure_limit_gamma_one():
    g = G.gamma_hanna(["CCCCCCCC", "Cc1ccccc1"], [1 - 1e-6, 1e-6], 300.0)
    assert abs(g[0] - 1.0) < 1e-3
    g = G.gamma_hanna(["CCCCCCCC", "Cc1ccccc1"], [1e-6, 1 - 1e-6], 300.0)
    assert abs(g[1] - 1.0) < 1e-3


def test_hanna_rejects_ternary():
    with pytest.raises(NotImplementedError):
        G.gamma_hanna(["C", "CC", "CCC"], [0.3, 0.3, 0.4], 300.0)
    with pytest.raises(NotImplementedError):
        G.make_hanna_gamma_fn(["C", "CC", "CCC"])
