"""Tests for scripts.v9.refit_antoine_v9.

Verifies:
- Fewer than 5 points -> failed:too_few_points
- Range missing the 40-100 degC flash window -> failed:missing_flash_window (trap 3)
- Falling or flat pressure data does not crash -> failed:falling_pressure
- WebBook constant converts with +2 (A_kPa = A_bar + 2)
- Valid data covering the flash window refits to refit_measured
"""
import numpy as np
import pytest

from scripts.v9.refit_antoine_v9 import (
    refit_compound,
    convert_webbook_antoine,
    antoine_log10_p,
)


def test_webbook_constant_converts_with_plus_two():
    """WebBook Antoine parameter A in bar converts to kPa with +2."""
    A_bar, B, C = 4.16082, 1826.948, -78.148
    A_kpa, B_kpa, C_kpa = convert_webbook_antoine(A_bar, B, C)
    assert pytest.approx(A_kpa, rel=1e-6) == 6.16082
    assert B_kpa == B
    assert C_kpa == C


def test_refit_too_few_points_fails():
    """Trap: fewer than 5 points gets failed:too_few_points."""
    T = np.array([320.0, 340.0, 360.0])
    P = np.array([0.1, 0.5, 1.5])
    res = refit_compound("test_comp", "CC", T, P)
    assert res["status"] == "failed:too_few_points"
    assert res["n_points"] == 3


def test_refit_missing_flash_window_fails():
    """Trap 3: range missing the 40-100 degC (313.15-373.15 K) flash window gets failed:missing_flash_window."""
    # Data high above 100 degC only (450 to 520 K)
    T = np.array([450.0, 465.0, 480.0, 495.0, 510.0, 520.0])
    P = np.array([10.0, 18.0, 30.0, 50.0, 80.0, 101.3])
    res = refit_compound("test_high_T", "CC", T, P)
    assert res["status"] == "failed:missing_flash_window"


def test_refit_falling_or_flat_pressure_does_not_crash():
    """Falling or flat pressure data does not crash scipy and gets failed:falling_pressure."""
    # Flat pressure
    T_flat = np.array([315.0, 325.0, 335.0, 345.0, 355.0, 365.0])
    P_flat = np.array([5.0, 5.0, 5.0, 5.0, 5.0, 5.0])
    res_flat = refit_compound("test_flat", "CC", T_flat, P_flat)
    assert res_flat["status"] == "failed:falling_pressure"

    # Inverted / falling pressure
    T_fall = np.array([315.0, 325.0, 335.0, 345.0, 355.0, 365.0])
    P_fall = np.array([10.0, 8.0, 6.0, 4.0, 2.0, 1.0])
    res_fall = refit_compound("test_fall", "CC", T_fall, P_fall)
    assert res_fall["status"] == "failed:falling_pressure"


def test_refit_valid_compound_succeeds():
    """Valid data covering the 40-100 degC flash window produces refit_measured."""
    # Synthetic clean data generated from known Antoine
    T = np.linspace(300.0, 380.0, 9)
    # Known Antoine: A=6.2, B=1800, C=-80
    P = 10.0 ** antoine_log10_p(T, 6.2, 1800.0, -80.0)
    res = refit_compound("test_valid", "Cc1ccccc1", T, P)
    assert res["status"] == "refit_measured"
    assert pytest.approx(res["A_kPa"], abs=0.1) == 6.2
    assert pytest.approx(res["B"], abs=20.0) == 1800.0
    assert pytest.approx(res["C"], abs=5.0) == -80.0
    assert res["rms_log10"] < 0.01
