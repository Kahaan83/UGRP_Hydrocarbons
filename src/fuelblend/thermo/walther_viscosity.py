"""ASTM D341 / Walther Temperature-Viscosity Curvature Model.

Implements the standard ASTM D341 double-logarithmic viscosity equation:
    log10(log10(nu + 0.7)) = A - B * log10(T)
where nu is kinematic viscosity in cSt (mm2/s) and T is absolute temperature in K.

Replaces the crude linear 1000/T assumption with standard industry temperature curvature,
substantially reducing cryogenic prediction error at -40 C.
"""
from __future__ import annotations

import math
from typing import Sequence
import numpy as np


def walther_viscosity(T_K: float, A: float, B: float) -> float:
    """Compute kinematic viscosity in cSt using Walther ASTM D341 equation."""
    logT = math.log10(T_K)
    Y = A - B * logT
    Z = 10.0 ** (10.0 ** Y)
    return max(Z - 0.7, 0.1)


def fit_walther_parameters(T_K_list: Sequence[float], nu_cSt_list: Sequence[float]) -> tuple[float, float]:
    """Fit Walther parameters A and B from temperature-viscosity pairs."""
    T = np.asarray(T_K_list, float)
    nu = np.asarray(nu_cSt_list, float)
    
    ok = (T > 100.0) & (nu > 0.0)
    T, nu = T[ok], nu[ok]
    if len(T) < 2:
        raise ValueError("At least 2 points required to fit Walther parameters.")
        
    Z = nu + 0.7
    Y = np.log10(np.log10(Z))
    X = np.log10(T)
    
    # Linear fit: Y = A - B * X  =>  Y = A + (-B) * X
    p = np.polyfit(X, Y, 1)
    B = -float(p[0])
    A = float(p[1])
    return A, B


def extrapolate_cryogenic_walther(
    T_target_C: float,
    T1_C: float,
    nu1_cSt: float,
    T2_C: float,
    nu2_cSt: float,
) -> float:
    """Extrapolate kinematic viscosity to cryogenic temperatures along the ASTM D341 line."""
    T1_K = T1_C + 273.15
    T2_K = T2_C + 273.15
    T_target_K = T_target_C + 273.15
    
    A, B = fit_walther_parameters([T1_K, T2_K], [nu1_cSt, nu2_cSt])
    return walther_viscosity(T_target_K, A, B)
