import importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("v7_baseline_d", ROOT / "scripts/v7/freezing/baseline_d.py")
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)


def test_limiting_species_single_component():
    # a pure compound "fuel" (x = 1) returns its own melting point (degC)
    assert abs(B.limiting_species_C(np.array([1.0]), [(250.0, 20000.0)]) - (250.0 - 273.15)) < 1e-9


def test_limiting_species_is_max():
    x = np.array([0.6, 0.3, 0.1])
    consts = [(240.0, 15000.0), (260.0, 25000.0), (300.0, 30000.0)]
    t = [1.0 / (1.0 / tm - B.R * np.log(xi) / dh) for xi, (tm, dh) in zip(x, consts)]
    assert abs(B.limiting_species_C(x, consts) - (max(t) - 273.15)) < 1e-9
    # species without constants do not crystallise and are skipped
    assert abs(B.limiting_species_C(x, [None, consts[1], None]) - (t[1] - 273.15)) < 1e-9
    assert np.isnan(B.limiting_species_C(x, [None, None, None]))
