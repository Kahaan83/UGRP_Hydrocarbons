#!/usr/bin/env python
"""Refit Antoine constants for candidate design compounds (v9 Target 3).

Inputs:
  data/raw/pure_v9/vapour_pressure_draft.csv
  results/v5/BASE/delta_representatives/candidate_gaps.csv
Output:
  results/v9/antoine_refit/refits.csv

Columns:
  name, smiles, A_kPa, B, C, T_min_K, T_max_K, n_points, rms_log10, status

Rules and Traps:
- Units: log10(P / kPa) = A - B / (T / K + C).
- NIST WebBook constants are in bar: A_kPa = A_bar + 2.
- Flash window: 313.15 to 373.15 K (40 to 100 degC).
- Trap 3: A refit with fewer than 5 points or a range missing the 40-100 degC flash window gets 'failed:'.
- Falling or flat pressure data does not crash (caught and marked 'failed:falling_pressure').
- Starting guesses clipped to bounds LO = [-5, 50, -250], HI = [20, 12000, 150].
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
RAW_VP = ROOT / "data" / "raw" / "pure_v9" / "vapour_pressure_draft.csv"
CANDIDATES_CSV = ROOT / "results" / "v5" / "BASE" / "delta_representatives" / "candidate_gaps.csv"
OUT_DIR = ROOT / "results" / "v9" / "antoine_refit"
OUT_CSV = OUT_DIR / "refits.csv"

FLASH_WINDOW = (313.15, 373.15)
LO = np.array([-5.0, 50.0, -250.0])
HI = np.array([20.0, 12000.0, 150.0])


def antoine_log10_p(T: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
    """log10(P/kPa) = A - B / (T/K + C)."""
    return A - B / (np.asarray(T, dtype=float) + C)


def convert_webbook_antoine(A_bar: float, B: float, C: float) -> tuple[float, float, float]:
    """Convert NIST WebBook Antoine parameters from bar to kPa: A_kPa = A_bar + 2."""
    return float(A_bar + 2.0), float(B), float(C)


def fit_antoine_parameters(
    T: np.ndarray,
    P: np.ndarray,
    weights: np.ndarray | None = None
) -> tuple[float, float, float, float]:
    """Fit 3-parameter Antoine equation to (T, P_kPa) data.
    
    Returns:
        (A, B, C, rms_log10)
    Raises:
        ValueError if data is non-physical or fitting fails.
    """
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)
    ok = np.isfinite(T) & np.isfinite(P) & (P > 0)
    T, P = T[ok], P[ok]
    
    if len(T) < 5:
        raise ValueError("too_few_points: fewer than 5 valid points")
    
    logP = np.log10(P)
    w = np.ones_like(T) if weights is None else np.asarray(weights, dtype=float)[ok]

    # Check for flat or falling pressure with increasing temperature
    if np.ptp(T) <= 0 or np.ptp(logP) <= 0:
        raise ValueError("falling_pressure: zero temperature or pressure range")
    
    corr = np.corrcoef(T, logP)[0, 1]
    if np.isnan(corr) or corr <= 0.5:
        raise ValueError("falling_pressure: vapour pressure does not increase with temperature")

    best_res = None
    for C0 in (-120.0, -80.0, -50.0, -20.0, 0.0):
        try:
            slope, icpt = np.polyfit(1.0 / (T + C0), logP, 1)
            p0 = np.clip([icpt, -slope, C0], LO + 1e-5, HI - 1e-5)
            res = least_squares(
                lambda p: w * (antoine_log10_p(T, *p) - logP),
                p0,
                bounds=(LO, HI),
                loss="soft_l1",
                f_scale=0.02
            )
            if best_res is None or res.cost < best_res.cost:
                best_res = res
        except Exception:
            continue

    if best_res is None or not best_res.success:
        raise ValueError("optimization_failed: no starting point converged")

    A, B, C = best_res.x
    pred_logP = antoine_log10_p(T, A, B, C)
    rms = float(np.sqrt(np.mean((pred_logP - logP) ** 2)))
    return float(A), float(B), float(C), rms


def refit_compound(
    name: str,
    smiles: str,
    T: np.ndarray,
    P: np.ndarray
) -> dict:
    """Evaluate and refit a single compound against project criteria."""
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)
    ok = np.isfinite(T) & np.isfinite(P) & (P > 0)
    T_clean, P_clean = T[ok], P[ok]
    n_pts = len(T_clean)

    if n_pts < 5:
        return {
            "name": name,
            "smiles": smiles,
            "A_kPa": np.nan,
            "B": np.nan,
            "C": np.nan,
            "T_min_K": float(T_clean.min()) if n_pts else np.nan,
            "T_max_K": float(T_clean.max()) if n_pts else np.nan,
            "n_points": n_pts,
            "rms_log10": np.nan,
            "status": "failed:too_few_points"
        }

    t_min = float(T_clean.min())
    t_max = float(T_clean.max())

    # Check 40-100 degC (313.15 to 373.15 K) flash window coverage (Trap 3)
    # Range must either encompass the window or have points within it
    covers_window = (t_min <= FLASH_WINDOW[0] + 5.0) and (t_max >= FLASH_WINDOW[1] - 5.0)
    has_window_points = np.any((T_clean >= FLASH_WINDOW[0]) & (T_clean <= FLASH_WINDOW[1]))
    if not (covers_window or has_window_points):
        return {
            "name": name,
            "smiles": smiles,
            "A_kPa": np.nan,
            "B": np.nan,
            "C": np.nan,
            "T_min_K": t_min,
            "T_max_K": t_max,
            "n_points": n_pts,
            "rms_log10": np.nan,
            "status": "failed:missing_flash_window"
        }

    try:
        A, B, C, rms = fit_antoine_parameters(T_clean, P_clean)
        return {
            "name": name,
            "smiles": smiles,
            "A_kPa": A,
            "B": B,
            "C": C,
            "T_min_K": t_min,
            "T_max_K": t_max,
            "n_points": n_pts,
            "rms_log10": rms,
            "status": "refit_measured"
        }
    except ValueError as e:
        reason = str(e).split(":")[0]
        return {
            "name": name,
            "smiles": smiles,
            "A_kPa": np.nan,
            "B": np.nan,
            "C": np.nan,
            "T_min_K": t_min,
            "T_max_K": t_max,
            "n_points": n_pts,
            "rms_log10": np.nan,
            "status": f"failed:{reason}"
        }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Load candidate gaps
    if not CANDIDATES_CSV.exists():
        raise FileNotFoundError(f"Missing {CANDIDATES_CSV}")
    candidates_df = pd.read_csv(CANDIDATES_CSV)

    # Load draft vapour pressure
    if not RAW_VP.exists():
        raise FileNotFoundError(f"Missing {RAW_VP}")
    vp_df = pd.read_csv(RAW_VP)

    results = []
    for row in candidates_df.itertuples():
        c_name = row.name
        c_smi = row.smiles
        c_vp = vp_df[vp_df["name"] == c_name]
        
        T = c_vp["T_K"].to_numpy(dtype=float) if not c_vp.empty else np.array([])
        P = c_vp["P_kPa"].to_numpy(dtype=float) if not c_vp.empty else np.array([])
        
        res = refit_compound(c_name, c_smi, T, P)
        results.append(res)

    out_df = pd.DataFrame(results)
    out_df.to_csv(OUT_CSV, index=False)
    print(f"Antoine refit results saved to {OUT_CSV}:")
    print(out_df[["name", "n_points", "rms_log10", "status"]])


if __name__ == "__main__":
    main()
