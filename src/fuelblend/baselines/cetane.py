"""Physics and empirical blending baselines for Derived Cetane Number (DCN) (T2.15).

Models supported:
  1. Linear-by-mole blending: DCN_mix = sum(x_i * DCN_i)
  2. Linear-by-mass blending: DCN_mix = sum(w_i * DCN_i)
  3. Linear-by-volume blending: DCN_mix = sum(v_i * DCN_i)
"""
from __future__ import annotations

from typing import Sequence, Dict, Any
import numpy as np


def blend_dcn_linear(fractions: Sequence[float], pure_dcn: Sequence[float]) -> float:
    """Compute linear combination of pure component DCN values."""
    f = np.asarray(fractions, dtype=float)
    c = np.asarray(pure_dcn, dtype=float)
    if len(f) != len(c):
        raise ValueError("Lengths of fractions and pure_dcn must match.")
    f_norm = f / f.sum()
    return float(np.dot(f_norm, c))


def predict_mixture_dcn(
    pure_dcn: Sequence[float],
    x_mole: Sequence[float] | None = None,
    w_mass: Sequence[float] | None = None,
    v_vol: Sequence[float] | None = None,
    basis: str = "mole",
) -> float:
    """Predict mixture DCN given pure component DCNs and fractions in requested basis.

    Parameters:
      pure_dcn: Pure component DCN values
      x_mole: Mole fractions
      w_mass: Mass fractions
      v_vol: Volume fractions
      basis: 'mole', 'mass', or 'volume'
    """
    b = basis.lower()
    if b == "mole":
        if x_mole is None:
            raise ValueError("x_mole required for mole basis.")
        return blend_dcn_linear(x_mole, pure_dcn)
    elif b == "mass":
        if w_mass is None:
            raise ValueError("w_mass required for mass basis.")
        return blend_dcn_linear(w_mass, pure_dcn)
    elif b == "volume":
        if v_vol is None:
            # If volume not explicitly provided, fall back to mass
            if w_mass is not None:
                return blend_dcn_linear(w_mass, pure_dcn)
            raise ValueError("v_vol or w_mass required for volume basis.")
        return blend_dcn_linear(v_vol, pure_dcn)
    else:
        raise ValueError(f"Unknown basis: {basis}")
