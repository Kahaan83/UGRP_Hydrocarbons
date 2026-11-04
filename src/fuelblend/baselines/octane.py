"""Octane (RON, MON, Octane Sensitivity) Blending Baselines (T2.17, T2.18).

Baselines supported:
  1. Linear-by-volume blending:
      ON_mix = sum(v_i * ON_i)
  2. Non-linear interaction-based blending (BON - Blending Octane Number):
      ON_mix = sum(v_i * ON_i) + sum_{i < j} v_i * v_j * I_{ij}
      where I_{ij} captures paraffin-aromatic or paraffin-olefin synergistic/antagonistic octane response.
  3. Derived Octane Sensitivity:
      Sensitivity = RON - MON
"""
from __future__ import annotations

from typing import Sequence, Dict, Any
import numpy as np


def blend_octane_linear_volume(v_vol: Sequence[float], pure_on: Sequence[float]) -> float:
    """Linear-by-volume octane blending rule (ASTM D2699/D2700 standard baseline)."""
    v = np.asarray(v_vol, dtype=float)
    on = np.asarray(pure_on, dtype=float)
    if len(v) != len(on):
        raise ValueError("Lengths of v_vol and pure_on must match.")
    v_norm = v / v.sum()
    return float(np.dot(v_norm, on))


def blend_octane_interaction(
    v_vol: Sequence[float],
    pure_on: Sequence[float],
    class_labels: Sequence[str],
    interaction_matrix: Dict[tuple, float] | None = None,
) -> float:
    """Interaction-based non-linear blending octane number (BON) model.

    Captures synergistic boosting (e.g. aromatic + alkane) and antagonistic blending.
    """
    linear_val = blend_octane_linear_volume(v_vol, pure_on)
    v = np.asarray(v_vol, dtype=float) / np.sum(v_vol)
    n = len(v)

    if interaction_matrix is None:
        # Default interaction parameters based on literature (e.g. Ghosh et al., Anderson et al.)
        interaction_matrix = {
            ("aromatic", "n_paraffin"): 3.5,     # Synergistic boost
            ("aromatic", "iso_paraffin"): 1.8,
            ("alkene_alkyne", "n_paraffin"): 2.2,
            ("cycloalkane", "aromatic"): -1.0,
        }

    interaction_term = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            pair_key = tuple(sorted([class_labels[i], class_labels[j]]))
            coeff = interaction_matrix.get(pair_key, 0.0)
            interaction_term += v[i] * v[j] * coeff

    return float(linear_val + interaction_term)


def calculate_octane_sensitivity(ron: float, mon: float) -> float:
    """Calculate derived octane sensitivity (RON - MON)."""
    return float(ron - mon)
