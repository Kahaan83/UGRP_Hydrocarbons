"""Unit tests for Heat of Combustion Physics Baseline (T2.2)."""
import pytest
import numpy as np

from fuelblend.baselines.heat_of_combustion import (
    mass_fractions_from_mole,
    mole_fractions_from_mass,
    gcv_to_ncv,
    ncv_to_gcv,
    mixture_heat_of_combustion,
    LATENT_HEAT_WATER_FACTOR,
)


def test_basis_conversion_roundtrip():
    # Equimolar binary of heptane (MW ~ 100.2) and toluene (MW ~ 92.14)
    x_mole = [0.5, 0.5]
    mw = [100.205, 92.141]

    w_mass = mass_fractions_from_mole(x_mole, mw)
    assert np.isclose(w_mass.sum(), 1.0)
    assert w_mass[0] > w_mass[1]  # Heptane is heavier

    # Roundtrip back to mole
    x_recov = mole_fractions_from_mass(w_mass, mw)
    assert np.allclose(x_mole, x_recov)


def test_gross_to_net_conversion():
    gcv = 46.20  # MJ/kg
    w_H = 14.0   # wt% H
    expected_ncv = gcv - LATENT_HEAT_WATER_FACTOR * w_H
    calc_ncv = gcv_to_ncv(gcv, w_H)
    assert np.isclose(calc_ncv, expected_ncv)

    # Reverse
    recov_gcv = ncv_to_gcv(calc_ncv, w_H)
    assert np.isclose(recov_gcv, gcv)


def test_mixture_heat_of_combustion_mole_vs_mass():
    # Pure A: 48 MJ/kg, MW=100
    # Pure B: 42 MJ/kg, MW=200
    gcv_pure = [48.0, 42.0]
    mw = [100.0, 200.0]
    w_H_pure = [15.0, 10.0]

    # Equal mass
    w_mass = [0.5, 0.5]
    res_mass = mixture_heat_of_combustion(w_mass, gcv_pure, mw, pure_w_H=w_H_pure, basis="mass")
    assert np.isclose(res_mass["gcv_MJ_kg"], 45.0)  # (48 + 42)/2
    assert np.isclose(res_mass["w_H_wt_pct"], 12.5) # (15 + 10)/2

    # In mole basis, equal mass corresponds to 2:1 molar ratio (since MW_B = 2 * MW_A)
    x_mole = [2.0 / 3.0, 1.0 / 3.0]
    res_mole = mixture_heat_of_combustion(x_mole, gcv_pure, mw, pure_w_H=w_H_pure, basis="mole")
    assert np.isclose(res_mole["gcv_MJ_kg"], res_mass["gcv_MJ_kg"])
    assert np.isclose(res_mole["ncv_MJ_kg"], res_mass["ncv_MJ_kg"])
