import numpy as np
import pytest

from fuelblend.thermo.antoine_interp import T_GRID, grid_interp, log10p, param_interp, refit_abc, rms_log10


def _series(n):
    # smooth synthetic homologous series (kPa form), roughly n-alkane-like
    return (5.9 + 0.06 * n, 900.0 + 95.0 * n, -20.0 - 4.0 * n)


def test_refit_recovers_curve():
    abc = _series(12)
    fit = refit_abc(T_GRID, log10p(abc, T_GRID))
    assert np.max(np.abs(log10p(fit, T_GRID) - log10p(abc, T_GRID))) < 1e-4


def test_grid_interp_recovers_smooth_series_across_a_4_carbon_gap():
    knots = {n: _series(n) for n in (8, 9, 10, 12, 16, 20, 21)}
    abc, lp = grid_interp(knots, 14)
    truth = log10p(_series(14), T_GRID)
    assert np.sqrt(np.mean((log10p(abc, T_GRID) - truth) ** 2)) < 0.03
    assert lp.shape == T_GRID.shape


def test_knots_never_overwritten_and_no_extrapolation():
    knots = {n: _series(n) for n in (8, 10, 12)}
    with pytest.raises(ValueError):
        grid_interp(knots, 10)
    with pytest.raises(ValueError):
        grid_interp(knots, 14)
    with pytest.raises(ValueError):
        param_interp(knots, 8)


def test_rms_log10_ignores_nonpositive():
    abc = _series(10)
    T = np.array([300.0, 350.0, 400.0])
    P = 10 ** log10p(abc, T)
    P[1] = 0.0
    assert rms_log10(abc, T, P) < 1e-12
