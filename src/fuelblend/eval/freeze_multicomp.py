"""Liquidus models compared in v5 delta_ws2_freeze_multicomp (descriptive; the designated freezing model stays ideal SvL).

All take mole fractions ``x``, SMILES and per-component constants ``consts`` = [(tm_K, dh_J_mol) or None, ...] aligned to the components.
A component with None does not crystallise (same treatment as ``composition_route.ideal_sle_C``).  Results in degC.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem
from scipy.optimize import brentq

R = 8.314462618
X_MIN = 1e-4


def class_mix(smiles) -> str:
    """alkane_only | alkane_naphthene | with_aromatic (declared in the delta README)."""
    has_arom = has_ring_sat = False
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if any(a.GetIsAromatic() for a in m.GetAtoms()):
            has_arom = True
        elif m.GetRingInfo().NumRings() > 0:
            has_ring_sat = True
    if has_arom:
        return "with_aromatic"
    return "alkane_naphthene" if has_ring_sat else "alkane_only"


def is_n_alkane(smiles: str) -> bool:
    m = Chem.MolFromSmiles(smiles)
    return (m.GetRingInfo().NumRings() == 0 and all(a.GetSymbol() == "C" for a in m.GetAtoms())
            and all(b.GetBondType() == Chem.BondType.SINGLE for b in m.GetBonds())
            and all(a.GetDegree() <= 2 for a in m.GetAtoms()))


def _sat_x(tm, dh, T):
    """Ideal solubility (mole fraction) of a pure solid at T."""
    return np.exp(-dh / R * (1.0 / T - 1.0 / tm))


def svl_ideal_C(x, consts) -> float:
    best = np.nan
    for xi, c in zip(x, consts):
        if c is None or xi <= X_MIN:
            continue
        tm, dh = c
        inv = 1.0 / tm - R / dh * np.log(xi)
        if inv > 0:
            T = 1.0 / inv - 273.15
            best = T if not np.isfinite(best) else max(best, T)
    return best


def svl_activity_C(x, consts, gamma_fn) -> float:
    """Pure solids, liquid activity from ``gamma_fn(T, x)``: per component solve ln(x_i g_i(T)) = -dH/R (1/T - 1/Tm)."""
    x = np.asarray(x, float)
    best = np.nan
    for i, c in enumerate(consts):
        if c is None or x[i] <= X_MIN:
            continue
        tm, dh = c

        def f(T):
            return np.log(x[i] * gamma_fn(T, x)[i]) + dh / R * (1.0 / T - 1.0 / tm)

        lo, hi = 60.0, tm + 60.0
        try:
            if f(lo) * f(hi) > 0:
                continue
            T = brentq(f, lo, hi, xtol=1e-4) - 273.15
        except (ValueError, ZeroDivisionError, FloatingPointError):
            continue
        best = T if not np.isfinite(best) else max(best, T)
    return best


def solid_solution_nalkanes_C(x, smiles, consts) -> float:
    """Labelled stand-in (not predictive UNIQUAC): n-alkanes crystallise as ONE ideal solid solution (liquidus where
    sum_i x_i / x_sat_i(T) = 1 over n-alkanes); every other component as a pure solid (ideal SvL).  Maximum of the two."""
    x = np.asarray(x, float)
    alk = [i for i, s in enumerate(smiles) if is_n_alkane(s) and consts[i] is not None and x[i] > 0]
    other = [consts[i] if i not in alk else None for i in range(len(x))]
    t_other = svl_ideal_C(x, other)
    t_ss = np.nan
    if alk:
        def g(T):
            return sum(x[i] / _sat_x(consts[i][0], consts[i][1], T) for i in alk) - 1.0
        hi = max(consts[i][0] for i in alk)
        try:
            if g(60.0) > 0 > g(hi):
                t_ss = brentq(g, 60.0, hi, xtol=1e-4) - 273.15
        except ValueError:
            pass
    vals = [v for v in (t_other, t_ss) if np.isfinite(v)]
    return max(vals) if vals else np.nan
