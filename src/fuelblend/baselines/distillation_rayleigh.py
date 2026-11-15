"""distillation_rayleigh.py - ASTM D86 as an ideal single-stage (Rayleigh) batch distillation.

Zero fitted parameters. Model (declared in results/v5/BASE/delta_d86_rayleigh/README.md):

* Liquid of n_i moles boils at its bubble point T_b(x):  sum_i x_i P_i^sat(T) = P  (Raoult, gamma = 1, P = 101.325 kPa),
  with log10(P_i^sat / kPa) = A_i - B_i / (T/K + C_i)  (the same Antoine sets the Engler route uses for Tb).
* The vapour leaving is in equilibrium with the liquid, y_i = x_i P_i^sat(T_b) / P, and is removed at once (one equilibrium stage, no reflux):
      d n_i / d N = y_i          (N = total liquid moles, decreasing).
* Recovered volume  V_d(N) = sum_i (n_i0 - n_i) * v_i,  v_i = M_i / rho_i(288.15 K)  (liquid molar volume at 15 C, ideal volumes).
* T_q = bubble point of the residual liquid when V_d / V_0 = q/100, for q in 10, 50, 90 and the FBP cut (99, as the Engler route).

Integrated with ``solve_ivp`` (RK45, tight tolerances) and exact event location at each cut, so the cut temperatures are smooth functions
of the recipe (the optimiser uses finite-difference gradients).
"""
from __future__ import annotations

import math
from typing import Dict, Sequence

import numpy as np
from scipy.integrate import solve_ivp

P_ATM_KPA = 101.325
LN10 = math.log(10.0)


def bubble_T_K(x: np.ndarray, A: np.ndarray, B: np.ndarray, C: np.ndarray, P_kPa: float = P_ATM_KPA, T0: float | None = None) -> float:
    """Bubble temperature (K) of liquid mole fractions x (Newton on ln sum x_i P_i)."""
    lnP = math.log(P_kPa)
    if T0 is None:
        tb = B / (A - math.log10(P_kPa)) - C
        T0 = float(np.dot(x, tb))
    T = T0
    for _ in range(60):
        lp = LN10 * (A - B / (T + C))
        m = lp.max()
        e = x * np.exp(lp - m)
        s = e.sum()
        g = math.log(s) + m - lnP
        dg = float(np.dot(e, LN10 * B / (T + C) ** 2) / s)
        step = g / dg
        step = max(min(step, 50.0), -50.0)
        T -= step
        if abs(step) < 1e-10:
            break
    return float(T)


def rayleigh_d86(A: Sequence[float], B: Sequence[float], C: Sequence[float], moles: Sequence[float], molar_volume: Sequence[float],
                 cuts_pct: Sequence[float] = (10.0, 50.0, 90.0, 99.0), P_kPa: float = P_ATM_KPA,
                 rtol: float = 1e-9, atol_rel: float = 1e-12) -> Dict[str, float]:
    """Cut temperatures (C) of an ideal Rayleigh distillation.  Returns {'T10': .., 'T50': .., ..., 'IBP': T at 0 %}."""
    A, B, C = (np.asarray(v, float) for v in (A, B, C))
    n0 = np.asarray(moles, float)
    v = np.asarray(molar_volume, float)
    keep = n0 > 0
    A, B, C, n0, v = A[keep], B[keep], C[keep], n0[keep], v[keep]
    N0 = n0.sum()
    n0 = n0 / N0                                 # normalise: N goes 1 -> 0
    V0 = float(np.dot(n0, v))
    out: Dict[str, float] = {"IBP": bubble_T_K(n0, A, B, C, P_kPa) - 273.15}
    if len(n0) == 1:
        for q in cuts_pct:
            out[f"T{int(q)}"] = out["IBP"]
        return out

    last = {"T": out["IBP"] + 273.15}

    def rhs(N, n):
        n = np.clip(n, 0.0, None)
        s = n.sum()
        x = n / s
        T = bubble_T_K(x, A, B, C, P_kPa, T0=last["T"])
        last["T"] = T
        y = x * np.exp(LN10 * (A - B / (T + C))) / P_kPa
        return y / y.sum()

    def make_event(q):
        def ev(N, n):
            return float(np.dot(n0 - n, v)) / V0 - q / 100.0
        ev.terminal = False
        ev.direction = 0
        return ev

    events = [make_event(q) for q in cuts_pct]
    N_end = 1e-6
    sol = solve_ivp(rhs, (1.0, N_end), n0, method="RK45", rtol=rtol, atol=atol_rel, events=events, dense_output=False)
    for q, ys in zip(cuts_pct, sol.y_events):
        if len(ys) == 0:
            out[f"T{int(q)}"] = float("nan")
            continue
        n = np.clip(ys[0], 0.0, None)
        out[f"T{int(q)}"] = bubble_T_K(n / n.sum(), A, B, C, P_kPa) - 273.15
    return out


def rayleigh_quantities(A, B, C, w_mass, mw, rho15, fbp_cut: float = 99.0) -> Dict[str, float]:
    """Mass fractions -> dict T10, T50, T90, FBP, T50_T10, T90_T10 (C).  rho15 in any consistent unit (only ratios matter)."""
    w = np.asarray(w_mass, float)
    mw = np.asarray(mw, float)
    rho = np.asarray(rho15, float)
    moles = w / mw
    vm = mw / rho
    r = rayleigh_d86(A, B, C, moles, vm, cuts_pct=(10.0, 50.0, 90.0, fbp_cut))
    q = dict(T10=r["T10"], T50=r["T50"], T90=r["T90"], FBP=r[f"T{int(fbp_cut)}"], IBP=r["IBP"])
    q["T50_T10"] = q["T50"] - q["T10"]
    q["T90_T10"] = q["T90"] - q["T10"]
    return q
