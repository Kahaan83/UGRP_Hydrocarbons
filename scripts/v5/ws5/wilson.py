"""WS5: Wilson / Margules activity coefficients from infinite-dilution activity coefficients (pure numpy, no GH-GNN import).

Wilson (Lambda_ii = 1):  ln g_i = 1 - ln(sum_j x_j L_ij) - sum_k x_k L_ki / sum_j x_j L_kj
Infinite dilution of i in j:  ln g_i_inf = 1 - ln L_ij - L_ji.
Given a = ln g_i_inf and b = ln g_j_inf, solve for u = ln L_ij, v = ln L_ji:
    f1 = 1 - u - exp(v) - a = 0,  f2 = 1 - v - exp(u) - b = 0.
Pairs without a solution (residual > TOL) get L = 1 (u = v = 0) and are counted (declared failure rule).
Two-parameter Margules (binary): ln g1 = x2^2 [A12 + 2 (A21 - A12) x1],  ln g2 = x1^2 [A21 + 2 (A12 - A21) x2], A12 = ln g1_inf, A21 = ln g2_inf.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import fsolve

TOL = 1e-8


def _resid(u, v, a, b):
    return 1.0 - u - np.exp(np.clip(v, -50, 50)) - a, 1.0 - v - np.exp(np.clip(u, -50, 50)) - b


def wilson_invert(a, b):
    """Vectorised inversion.  a, b arrays (ln g_i_inf, ln g_j_inf).  Returns (u, v, ok) with ok a bool array."""
    a = np.atleast_1d(np.asarray(a, float))
    b = np.atleast_1d(np.asarray(b, float))
    u = -0.5 * a + 1e-3
    v = -0.5 * b - 1e-3
    for _ in range(80):
        f1, f2 = _resid(u, v, a, b)
        eu, ev = np.exp(np.clip(u, -50, 50)), np.exp(np.clip(v, -50, 50))
        det = (-1.0) * (-1.0) - (-ev) * (-eu)             # J = [[-1, -ev], [-eu, -1]]
        det = np.where(np.abs(det) < 1e-12, 1e-12, det)
        du = (-1.0 * (-f1) - (-ev) * (-f2)) / det           # J^-1 (-f)
        dv = ((-1.0) * (-f2) - (-eu) * (-f1)) / det
        n0 = np.maximum(np.abs(f1), np.abs(f2))
        step = np.ones_like(u)
        for _k in range(8):                                 # backtracking
            un, vn = u + step * du, v + step * dv
            g1, g2 = _resid(un, vn, a, b)
            worse = np.maximum(np.abs(g1), np.abs(g2)) > n0
            if not worse.any():
                break
            step = np.where(worse, step * 0.5, step)
        u, v = u + step * du, v + step * dv
        if np.max(np.maximum(np.abs(f1), np.abs(f2))) < TOL * 0.1:
            break
    f1, f2 = _resid(u, v, a, b)
    ok = (np.maximum(np.abs(f1), np.abs(f2)) <= TOL) & np.isfinite(u) & np.isfinite(v)
    for i in np.where(~ok)[0]:                              # fsolve restarts
        for s in ((0.5, 0.5), (-0.5, -0.5), (1.0, -1.0), (-1.0, 1.0), (2.0, 2.0)):
            sol, info, ier, _ = fsolve(lambda z: _resid(z[0], z[1], a[i], b[i]), s, full_output=True, xtol=1e-12)
            r = np.max(np.abs(_resid(sol[0], sol[1], a[i], b[i])))
            if ier == 1 and r <= TOL:
                u[i], v[i], ok[i] = sol[0], sol[1], True
                break
    u = np.where(ok, u, 0.0)
    v = np.where(ok, v, 0.0)
    return u, v, ok


def wilson_gamma(x, lnL):
    """x: (n,), lnL: (n, n) with lnL[i, j] = ln Lambda_ij (diagonal 0).  Returns gamma (n,)."""
    x = np.asarray(x, float)
    L = np.exp(lnL)
    S = L @ x                                               # S_i = sum_j L_ij x_j
    ln_g = 1.0 - np.log(S) - L.T @ (x / S)
    return np.exp(ln_g)


def margules_gamma_binary(x, a12, a21):
    x1, x2 = float(x[0]), float(x[1])
    lg1 = x2 ** 2 * (a12 + 2.0 * (a21 - a12) * x1)
    lg2 = x1 ** 2 * (a21 + 2.0 * (a12 - a21) * x2)
    return np.exp([lg1, lg2])


class PairTable:
    """ln g_ij_inf(T) = A_ij + B_ij / T_K for ordered pairs (solute i in solvent j) keyed by SMILES strings."""

    def __init__(self, A: dict, B: dict):
        self.A, self.B = A, B

    def lng_inf(self, si: str, sj: str, T_K: float) -> float:
        return self.A[(si, sj)] + self.B[(si, sj)] / T_K

    def lnL_matrix(self, smiles, T_K: float, ad_ok: dict | None = None):
        """Returns (lnL, n_fail, n_pairs).  ad_ok: optional {(frozenset or tuple pair): bool}; pairs flagged False are set to L = 1."""
        n = len(smiles)
        iu, ju = np.triu_indices(n, 1)
        a = np.array([self.lng_inf(smiles[i], smiles[j], T_K) for i, j in zip(iu, ju)])
        b = np.array([self.lng_inf(smiles[j], smiles[i], T_K) for i, j in zip(iu, ju)])
        u, v, ok = wilson_invert(a, b)
        if ad_ok is not None:
            keep = np.array([bool(ad_ok.get(frozenset((smiles[i], smiles[j])), False)) for i, j in zip(iu, ju)])
            u, v = np.where(keep, u, 0.0), np.where(keep, v, 0.0)
        lnL = np.zeros((n, n))
        lnL[iu, ju] = u
        lnL[ju, iu] = v
        return lnL, int((~ok).sum()), int(len(iu))
