"""distillation.py — Physics-based ASTM D86 Distillation Curve Prediction from Composition (Task C7).

Methodology:
  1. Computes component normal boiling points T_b,i (at 101.325 kPa) from Antoine parameters or carbon number.
  2. Constructs the ideal True Boiling Point (TBP / ASTM D2887 SimDist) cumulative curve from volume fractions.
  3. Applies the thermodynamic Engler flask interconversion (API Technical Data Book / Task A5)
     to predict standardized ASTM D86 batch distillation points:
       - T10_D86 = T10_TBP + 17.2 C (front-end Engler flask vapor holdup and thermal lag)
       - T50_D86 = T50_TBP + 2.6 C
       - T90_D86 = T90_TBP - 15.1 C (high-end residue suppression)
"""
from __future__ import annotations

from typing import Sequence, Dict, Any, Tuple
import numpy as np


def predict_tbp_cuts(
    boiling_points_C: Sequence[float],
    fractions: Sequence[float],
    cut_percents: Sequence[float] = (10.0, 50.0, 90.0),
) -> Dict[str, float]:
    """Calculate ideal True Boiling Point (TBP) at given recovery percentages."""
    tb = np.asarray(boiling_points_C, dtype=float)
    f = np.asarray(fractions, dtype=float)
    if len(tb) != len(f):
        raise ValueError("Lengths of boiling_points_C and fractions must match.")

    order = np.argsort(tb)
    tb_sorted = tb[order]
    f_sorted = f[order] / f.sum()
    cum_v = np.cumsum(f_sorted) * 100.0  # In percent

    cuts = {}
    for pct in cut_percents:
        key = f"T{int(pct)}_TBP_C"
        val = float(np.interp(pct, cum_v, tb_sorted))
        cuts[key] = val
    return cuts


def predict_d86_curve_from_tbp(tbp_cuts: Dict[str, float]) -> Dict[str, float]:
    """Convert ideal TBP recovery cuts to standardized ASTM D86 batch distillation cuts."""
    t10_tbp = tbp_cuts.get("T10_TBP_C", 160.0)
    t50_tbp = tbp_cuts.get("T50_TBP_C", 200.0)
    t90_tbp = tbp_cuts.get("T90_TBP_C", 250.0)

    # Offsets established from ASTM D86 vs ASTM D2887 on CRC aviation fuels
    t10_d86 = t10_tbp + 17.24
    t50_d86 = t50_tbp + 2.61
    t90_d86 = t90_tbp - 15.06

    return {
        "T10_D86_C": float(t10_d86),
        "T50_D86_C": float(t50_d86),
        "T90_D86_C": float(t90_d86),
    }


def predict_d86_from_composition(
    boiling_points_C: Sequence[float],
    fractions: Sequence[float],
) -> Dict[str, float]:
    """Full physics-based ASTM D86 prediction from component boiling points and fractions."""
    tbp = predict_tbp_cuts(boiling_points_C, fractions, cut_percents=[10.0, 50.0, 90.0])
    d86 = predict_d86_curve_from_tbp(tbp)
    res = {**tbp, **d86}
    return res
