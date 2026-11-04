"""Physics baselines for Net and Gross Heat of Combustion (T2.2).

Governing physics:
  - First Law of Thermodynamics: Enthalpy is an extensive state function.
  - Heat of combustion of ideal mixtures blends linearly on a mass basis:
      Delta H_comb,mix = sum_i (w_i * Delta H_comb,i)
  - Gross Calorific Value (GCV / HHV, water liquid) to Net Calorific Value (NCV / LHV, water vapor):
      NCV = GCV - (h_fg,H2O * (m_H2O / m_fuel))
      NCV = GCV - 0.2122 * w_H  (with w_H in wt% hydrogen, constant = 0.2122 MJ/kg per wt% H)
"""
from __future__ import annotations

from typing import Sequence, Dict, Any, Union
import numpy as np


LATENT_HEAT_WATER_FACTOR = 0.2122  # MJ/kg per wt% hydrogen in fuel (ASTM D3338 / D4809)
DEFAULT_JET_HYDROGEN_WT_PCT = 13.8  # Standard default for aviation turbine fuel


def mass_fractions_from_mole(x_mole: Sequence[float], mw: Sequence[float]) -> np.ndarray:
    """Convert mole fractions to mass fractions."""
    x = np.asarray(x_mole, dtype=float)
    m = np.asarray(mw, dtype=float)
    if len(x) != len(m):
        raise ValueError("Lengths of x_mole and mw must match.")
    mass = x * m
    total_mass = mass.sum()
    if total_mass <= 0:
        raise ValueError("Total molecular mass must be positive.")
    return mass / total_mass


def mole_fractions_from_mass(w_mass: Sequence[float], mw: Sequence[float]) -> np.ndarray:
    """Convert mass fractions to mole fractions."""
    w = np.asarray(w_mass, dtype=float)
    m = np.asarray(mw, dtype=float)
    if len(w) != len(m):
        raise ValueError("Lengths of w_mass and mw must match.")
    moles = w / m
    total_moles = moles.sum()
    if total_moles <= 0:
        raise ValueError("Total moles must be positive.")
    return moles / total_moles


def gcv_to_ncv(gcv_MJ_kg: float, w_H_wt_pct: float = DEFAULT_JET_HYDROGEN_WT_PCT) -> float:
    """Convert Gross Calorific Value (MJ/kg) to Net Calorific Value (MJ/kg)."""
    return float(gcv_MJ_kg - LATENT_HEAT_WATER_FACTOR * w_H_wt_pct)


def ncv_to_gcv(ncv_MJ_kg: float, w_H_wt_pct: float = DEFAULT_JET_HYDROGEN_WT_PCT) -> float:
    """Convert Net Calorific Value (MJ/kg) to Gross Calorific Value (MJ/kg)."""
    return float(ncv_MJ_kg + LATENT_HEAT_WATER_FACTOR * w_H_wt_pct)


def mixture_heat_of_combustion_mass(
    w_mass: Sequence[float],
    pure_gcv: Sequence[float],
) -> float:
    """Calculate mixture Gross Calorific Value (MJ/kg) from mass fractions."""
    w = np.asarray(w_mass, dtype=float)
    h = np.asarray(pure_gcv, dtype=float)
    if len(w) != len(h):
        raise ValueError("Lengths of w_mass and pure_gcv must match.")
    # Renormalize mass fractions if slight drift
    w_norm = w / w.sum()
    return float(np.dot(w_norm, h))


def mixture_heat_of_combustion(
    fractions: Sequence[float],
    pure_gcv: Sequence[float],
    mw: Sequence[float],
    pure_w_H: Sequence[float] | None = None,
    basis: str = "mass",
    target: str = "ncv",
) -> Dict[str, float]:
    """Complete physics baseline for mixture heat of combustion.

    Parameters:
      fractions: Component fractions in given basis
      pure_gcv: Pure component GCV (MJ/kg)
      mw: Component molecular weights (g/mol)
      pure_w_H: Pure component hydrogen content (wt%), optional
      basis: 'mass' or 'mole'
      target: 'gcv' (Gross/HHV) or 'ncv' (Net/LHV)

    Returns:
      {'gcv_MJ_kg': float, 'ncv_MJ_kg': float, 'w_H_wt_pct': float}
    """
    if basis.lower() == "mole":
        w = mass_fractions_from_mole(fractions, mw)
    elif basis.lower() == "mass":
        w = np.asarray(fractions, dtype=float)
        w = w / w.sum()
    else:
        raise ValueError(f"Unsupported composition basis: {basis}")

    gcv_mix = mixture_heat_of_combustion_mass(w, pure_gcv)

    if pure_w_H is not None:
        w_H_mix = float(np.dot(w, np.asarray(pure_w_H, dtype=float)))
    else:
        w_H_mix = DEFAULT_JET_HYDROGEN_WT_PCT

    ncv_mix = gcv_to_ncv(gcv_mix, w_H_mix)

    return {
        "gcv_MJ_kg": gcv_mix,
        "ncv_MJ_kg": ncv_mix,
        "w_H_wt_pct": w_H_mix,
    }
