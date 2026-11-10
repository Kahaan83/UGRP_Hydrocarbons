"""
Flash point baseline models using Le Chatelier's rule, Raoult's law, and Antoine vapor pressure.

Model:
    H(T) = sum_i [ x_i * P_i^sat(T) / (P_total * LFL_i(T)) ] = 1.0
where:
    P_total = 101.325 kPa
    P_i^sat(T) = 10 ** (A_i - B_i / (T_K + C_i))  [kPa]
    LFL_i(T) = LFL_i(25 C) * [1 - 0.000721 * (T_C - 25)]  (optional Zabetakis temperature correction)
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple
import numpy as np
from scipy.optimize import brentq

P_TOTAL_KPA = 101.325
ZABETAKIS_COEFF = 0.000721


def psat_antoine(T_K: float, A: float, B: float, C: float) -> float:
    """Compute saturation vapor pressure in kPa from Antoine parameters.
    
    Convention: log10(P/kPa) = A - B / (T_K + C)
    Guards against T_K + C <= 0.
    """
    d = T_K + C
    if d <= 1.0:
        return 0.0
    return float(10.0 ** (A - B / d))


def lfl_temperature_correction(lfl_25_vol_pct: float, T_K: float, min_factor: float = 0.05) -> float:
    """Compute temperature-dependent LFL (as volume fraction) using Zabetakis equation.
    
    LFL(T) = LFL(25 C) * max(1 - 0.000721 * (T_C - 25), min_factor)
    """
    T_C = T_K - 273.15
    factor = max(1.0 - ZABETAKIS_COEFF * (T_C - 25.0), min_factor)
    return (lfl_25_vol_pct / 100.0) * factor


def hazard_function(
    T_K: float,
    x_mole: Sequence[float],
    antoine_params: Sequence[Tuple[float, float, float]],
    lfl_vol_pct: Sequence[float],
    lfl_temp: bool = False,
    P_total: float = P_TOTAL_KPA,
    gamma_fn=None,
) -> float:
    """Compute Le Chatelier flammability hazard ratio H(T) - 1.0.
    
    H(T) = sum_i [ y_i / LFL_i ] = sum_i [ x_i * P_i^sat(T) / (P_total * LFL_i) ]
    Flash point occurs where H(T) = 1.0 (hazard_function == 0.0).

    ``gamma_fn(T_K, x_mole) -> sequence of activity coefficients`` (optional, v4 Task 9) replaces the
    ideal-solution assumption gamma_i = 1 (modified Raoult: y_i P = x_i gamma_i Psat_i).
    """
    gammas = [1.0] * len(x_mole) if gamma_fn is None else list(gamma_fn(T_K, x_mole))
    total_ratio = 0.0
    for xi, gi, (A, B, C), lfl_val in zip(x_mole, gammas, antoine_params, lfl_vol_pct):
        if xi <= 0.0:
            continue
        p_sat = psat_antoine(T_K, A, B, C) * float(gi)
        if lfl_temp:
            lfl_frac = lfl_temperature_correction(lfl_val, T_K)
        else:
            lfl_frac = lfl_val / 100.0
        
        if lfl_frac <= 0.0:
            continue
        total_ratio += (xi * p_sat) / (P_total * lfl_frac)
    return total_ratio - 1.0


def predict_flash_le_chatelier(
    x_mole: Sequence[float],
    antoine_params: Sequence[Tuple[float, float, float]],
    lfl_vol_pct: Sequence[float],
    lfl_temp: bool = False,
    lo_K: float = 150.0,
    hi_K: float = 650.0,
    gamma_fn=None,
) -> Tuple[float, str]:
    """Solve for flash point in deg C using Le Chatelier's rule.
    
    Returns:
        (predicted_flash_C, status_string)
    """
    if len(x_mole) != len(antoine_params) or len(x_mole) != len(lfl_vol_pct):
        return np.nan, "mismatched_component_lengths"
    if sum(x_mole) <= 0.0:
        return np.nan, "zero_mole_fraction_sum"
    
    # Normalize mole fractions
    tot = sum(x_mole)
    x_norm = [xi / tot for xi in x_mole]
    
    f = lambda T: hazard_function(T, x_norm, antoine_params, lfl_vol_pct, lfl_temp=lfl_temp, gamma_fn=gamma_fn)
    
    try:
        f_lo = f(lo_K)
        f_hi = f(hi_K)
        if not np.isfinite(f_lo) or f_lo > 0.0:
            return np.nan, f"below_bracket_{lo_K}K"
        if not np.isfinite(f_hi) or f_hi < 0.0:
            return np.nan, f"above_bracket_{hi_K}K"
        
        T_root = brentq(f, lo_K, hi_K, xtol=1e-4)
        return float(T_root - 273.15), "ok"
    except Exception as e:
        return np.nan, f"solver_error_{type(e).__name__}"


def check_antoine_validity(
    T_flash_C: float,
    antoine_Tlo: Sequence[float],
    antoine_Thi: Sequence[float],
    tolerance_K: float = 0.0,
) -> Tuple[bool, float]:
    """Check if predicted flash temperature is within Antoine validity limits.
    
    Returns:
        (in_range_bool, max_extrapolation_distance_K)
    """
    T_K = T_flash_C + 273.15
    max_dist = 0.0
    in_range = True
    
    for t_lo, t_hi in zip(antoine_Tlo, antoine_Thi):
        if np.isnan(t_lo) or np.isnan(t_hi):
            in_range = False
            continue
        if T_K < (t_lo - tolerance_K):
            in_range = False
            dist = (t_lo - tolerance_K) - T_K
            if dist > max_dist:
                max_dist = dist
        elif T_K > (t_hi + tolerance_K):
            in_range = False
            dist = T_K - (t_hi + tolerance_K)
            if dist > max_dist:
                max_dist = dist
    return in_range, max_dist
