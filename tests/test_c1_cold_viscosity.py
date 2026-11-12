"""v5 C1: offset fitting, bootstrap sign convention, non-finite exclusion, success rule, Walther reference (scripts/v5/c1_cold_viscosity_correction.py)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("c1_t", ROOT / "scripts/v5/c1_cold_viscosity_correction.py")
c1 = importlib.util.module_from_spec(spec)
sys.modules["c1_t"] = c1
spec.loader.exec_module(c1)


def test_offset_fitting_recovers_known_offsets():
    T = np.array([233.15, 233.15, 253.15, 253.15, 313.15, 313.15])
    true = 0.35 + (-0.10) * 1000.0 / T                      # ln(meas/pred) linear in 1000/T
    ab = c1.fit_k2(T, true)
    assert np.allclose(ab, (0.35, -0.10), atol=1e-9)
    assert np.allclose(c1.apply_k2(ab, [233.15, 400.0]), 0.35 - 0.10 * 1000.0 / np.array([233.15, 400.0]))
    k1 = c1.fit_k1(T, np.array([0.5, 0.7, 0.2, 0.3, -0.1, 0.1]))
    assert np.isclose(k1[233.15], 0.6) and np.isclose(k1[253.15], 0.25) and np.isclose(k1[313.15], 0.0)
    assert np.isnan(c1.apply_k1(k1, [273.15])[0])             # K1 is undefined at an unseen temperature
    assert c1.fit_k2([253.15, 253.15], [0.1, 0.2]) is None    # one temperature cannot give a slope


def test_bootstrap_positive_means_correction_better():
    g = np.arange(20)
    e0 = np.full(20, -0.6) + np.linspace(-0.05, 0.05, 20)     # raw: large error
    eK = np.linspace(-0.05, 0.05, 20)                         # corrected: small error
    gain, lo, hi, p = c1.gain_bootstrap(g, e0, eK, n=500, seed=0)
    assert gain > 0.4 and lo > 0 and hi > lo and p == 1.0
    gain2, lo2, hi2, _ = c1.gain_bootstrap(g, eK, e0, n=500, seed=0)   # swapped: correction worse -> negative
    assert gain2 < 0 and hi2 < 0
    assert c1.gain_bootstrap(g, e0, eK, n=500, seed=0) == c1.gain_bootstrap(g, e0, eK, n=500, seed=0)   # seeded


def test_non_finite_predictions_are_excluded_and_counted():
    ok, n = c1.split_finite([2.0, np.nan, np.inf, -1.0, 0.0, 5.0])
    assert ok.tolist() == [True, False, False, False, False, True] and n == 4


def test_success_rule_requires_n_interval_and_bias():
    assert c1.verdict(7, 0.2, 0.0)[0] == "n too small"
    assert c1.verdict(10, -0.01, 0.0)[1] is False            # interval includes 0
    assert c1.verdict(10, 0.2, -0.15)[1] is False            # residual bias too large
    assert c1.verdict(10, 0.2, 0.10)[1] is True              # boundary: |bias| <= 0.10 passes
    assert c1.verdict(8, 0.05, -0.02) == ("works", True)


def test_walther_line_recovers_and_extrapolates():
    A, B = 8.0, -1.4
    T = np.array([288.15, 313.15, 373.15])
    nu = np.exp(np.exp(A + B * np.log(T))) - 0.7
    assert np.isclose(c1.walther_extrapolate(T, nu, 253.15), np.exp(np.exp(A + B * np.log(253.15))) - 0.7, rtol=1e-9)
    assert np.isnan(c1.walther_extrapolate([313.15], [2.0], 253.15))
    assert np.isnan(c1.margin_for(300.0)) and c1.margin_for(253.15) == 0.4154 and c1.margin_for(273.15) == 0.6084
