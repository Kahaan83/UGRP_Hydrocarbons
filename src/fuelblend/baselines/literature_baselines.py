"""Classical literature physical mixing rules and baselines (T5.3).

Implements:
  1. Grunberg-Nissan mixing rule for liquid dynamic viscosity:
      ln(eta_mix) = sum_i (x_i * ln(eta_i)) + sum_{i < j} x_i * x_j * G_ij
  2. Kendall-Monroe mixing rule for liquid dynamic viscosity:
      eta_mix^(1/3) = sum_i (x_i * eta_i^(1/3))
  3. Redlich-Kister excess molar volume model for liquid density:
      V_E = x_1 * x_2 * [ A_0 + A_1 * (x_1 - x_2) ]
      rho_mix = (sum_i x_i * M_i) / (V_ideal + V_E)
"""
from __future__ import annotations

from typing import Sequence, Dict, Any
import numpy as np


def grunberg_nissan_viscosity(
    x_mole: Sequence[float],
    pure_eta_mPas: Sequence[float],
    G_interaction: float = 0.0,
) -> float:
    """Grunberg-Nissan mixing rule for liquid dynamic viscosity (mPa*s)."""
    x = np.asarray(x_mole, dtype=float)
    x = x / x.sum()
    eta = np.asarray(pure_eta_mPas, dtype=float)

    # Base log-linear term
    ln_eta_ideal = float(np.dot(x, np.log(eta)))

    # Interaction term for binary/multicomponent
    interact = 0.0
    n = len(x)
    if n >= 2 and G_interaction != 0.0:
        for i in range(n):
            for j in range(i + 1, n):
                interact += x[i] * x[j] * G_interaction

    return float(np.exp(ln_eta_ideal + interact))


def kendall_monroe_viscosity(
    x_mole: Sequence[float],
    pure_eta_mPas: Sequence[float],
) -> float:
    """Kendall-Monroe cube-root mixing rule for liquid dynamic viscosity (mPa*s).

    eta_mix^(1/3) = sum_i (x_i * eta_i^(1/3))
    """
    x = np.asarray(x_mole, dtype=float)
    x = x / x.sum()
    eta = np.asarray(pure_eta_mPas, dtype=float)
    cube_root_sum = np.dot(x, eta ** (1.0 / 3.0))
    return float(cube_root_sum ** 3.0)


def redlich_kister_density(
    x_mole: Sequence[float],
    mw: Sequence[float],
    pure_rho_kg_m3: Sequence[float],
    A0: float = 0.0,
    A1: float = 0.0,
) -> float:
    """Redlich-Kister excess molar volume model for liquid mixture density (kg/m3).

    V_ideal = sum_i (x_i * M_i / rho_i)
    V_E = x_1 * x_2 * [A0 + A1 * (x_1 - x_2)]  (in cm3/mol = 1e-6 m3/mol)
    rho_mix = (sum_i x_i * M_i) / (V_ideal + V_E)
    """
    x = np.asarray(x_mole, dtype=float)
    x = x / x.sum()
    m = np.asarray(mw, dtype=float)       # g/mol
    rho = np.asarray(pure_rho_kg_m3, dtype=float) # kg/m3 = g/cm3 * 1000 -> rho in g/cm3 is rho / 1000

    # Molar volumes V_i in cm3/mol
    rho_g_cm3 = rho / 1000.0
    V_i = m / rho_g_cm3
    V_ideal = float(np.dot(x, V_i)) # cm3/mol

    # Excess volume V_E in cm3/mol
    V_E = 0.0
    if len(x) == 2:
        V_E = x[0] * x[1] * (A0 + A1 * (x[0] - x[1]))

    V_total = V_ideal + V_E
    if V_total <= 0:
        V_total = V_ideal

    m_avg = float(np.dot(x, m)) # g/mol
    rho_mix_g_cm3 = m_avg / V_total
    return float(rho_mix_g_cm3 * 1000.0)
