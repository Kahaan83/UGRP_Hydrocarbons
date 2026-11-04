"""Thermodynamic Solid-Liquid Equilibrium (SLE) Baseline for Wax Appearance Temperature (T2.11).

Governing thermodynamics:
  Schroeder-van Laar ideal-solution SLE equation for precipitating solid solute i:
      ln(x_i) = -(Delta H_fus,i / R) * (1/T_precip,i - 1/T_m,i)
      1/T_precip,i = 1/T_m,i - (R / Delta H_fus,i) * ln(x_i)

The Wax Appearance Temperature (WAT / cloud point / incipient freezing point)
is the maximum precipitation temperature across all precipitating solutes:
      WAT_K = max_i (T_precip,i)
      WAT_C = WAT_K - 273.15
"""
from __future__ import annotations

from typing import Sequence, Dict, Any
import numpy as np

R_GAS = 8.314462618  # J / (mol * K)


def schroeder_van_laar_T(x_mole: float, tm_K: float, hfus_J_mol: float) -> float:
    """Calculate precipitation temperature (K) for a single precipitating solute."""
    if x_mole <= 1e-12 or not np.isfinite(tm_K) or not np.isfinite(hfus_J_mol) or hfus_J_mol <= 0:
        return 0.0
    inv_T = (1.0 / tm_K) - (R_GAS / hfus_J_mol) * np.log(x_mole)
    if inv_T <= 0:
        return 0.0
    return float(1.0 / inv_T)


def wax_appearance_temperature_C(
    x_mole: Sequence[float],
    tm_K: Sequence[float],
    hfus_J_mol: Sequence[float],
) -> float:
    """Calculate mixture Wax Appearance Temperature (WAT, degC) using ideal SLE.

    Parameters:
      x_mole: Mole fractions of components
      tm_K: Pure component melting temperatures (Kelvin)
      hfus_J_mol: Pure component fusion enthalpies (Joules/mole)

    Returns:
      WAT in degC (or -273.15 if no component precipitates)
    """
    x = np.asarray(x_mole, dtype=float)
    tm = np.asarray(tm_K, dtype=float)
    hf = np.asarray(hfus_J_mol, dtype=float)

    if len(x) != len(tm) or len(x) != len(hf):
        raise ValueError("Lengths of x_mole, tm_K, and hfus_J_mol must match.")

    # Normalize mole fractions
    x_norm = x / x.sum()

    t_precip = []
    for xi, tmi, hfi in zip(x_norm, tm, hf):
        if np.isfinite(tmi) and np.isfinite(hfi) and hfi > 0:
            tp = schroeder_van_laar_T(xi, tmi, hfi)
            if tp > 0:
                t_precip.append(tp)

    if not t_precip:
        return -273.15

    wat_K = float(np.max(t_precip))
    return float(wat_K - 273.15)
