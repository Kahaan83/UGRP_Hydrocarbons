"""Coutinho-Style Multi-Component Non-Ideal Wax Crystallization Model.

Replaces the crude flat offset with a composition-dependent eutectic depression
model calibrated on multi-component n-alkane systems (Dauphin et al., Coutinho SP-Wax).
Accounts for solid-solution formation and n-alkane carbon-number polydispersity:
    T_freeze = T_ideal_WAT - Delta_T_eutectic(sigma_NC, non_paraffin_fraction)
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
COMPOUNDS_PATH = ROOT / "data" / "unified" / "compounds.csv"
R_GAS = 8.314462618  # J / (mol * K)

# Empirical Coutinho eutectic depression parameters calibrated on Dauphin & CRC multi-component data
COUTINHO_ALPHA = 3.20   # K per unit of carbon-number standard deviation
COUTINHO_BETA = 14.50   # K for solvent dilution (non-paraffin matrix)
COUTINHO_GAMMA = 0.60   # baseline solid-solution depression


def get_pure_thermal_properties(comp_df: pd.DataFrame | None = None) -> dict[str, dict[str, float]]:
    """Return dictionary of inchikey -> {Tm_K, dHfus_J_mol, n_carbon, cls}."""
    if comp_df is None:
        comp_df = pd.read_csv(COMPOUNDS_PATH).set_index("inchikey")

    res = {}
    for ik, row in comp_df.iterrows():
        nc = float(row.get("n_carbon", 0.0))
        cls = str(row.get("cls", "other"))

        # Pure alkane melting enthalpy & temperature
        if cls == "n_paraffin" and nc >= 6:
            # Broadhurst / NIST correlations for n-alkanes
            # Refined for common alkanes up to C36
            known_tm = {
                6: 177.8, 7: 182.6, 8: 216.4, 9: 219.7, 10: 243.5, 11: 247.6, 12: 263.6,
                13: 267.8, 14: 279.0, 15: 283.1, 16: 291.3, 17: 295.1, 18: 301.3, 19: 305.0, 20: 309.6,
                21: 313.6, 22: 317.5, 23: 320.7, 24: 324.0, 25: 327.1, 26: 329.5, 27: 332.6, 28: 334.5,
                29: 336.8, 30: 338.9, 31: 341.1, 32: 342.8, 34: 345.9, 36: 349.1,
            }
            tm_default = 414.3 * (nc - 1.5) / (nc + 5.0)
            tm_k = known_tm.get(int(nc), tm_default)
            dh_fus = (4.0 * nc - 10.6) * 1000.0  # J/mol
        elif cls == "aromatic":
            smi = str(row.get("smiles", ""))
            # Known high-melting aromatics
            if "Cc1cc(C)c(C)cc1C" in smi: # durene
                tm_k, dh_fus = 352.4, 21000.0
            elif "c1ccc2ccccc2c1" in smi: # naphthalene
                tm_k, dh_fus = 353.4, 19000.0
            elif "Cc1ccc(C)cc1" in smi: # p-xylene
                tm_k, dh_fus = 286.4, 17100.0
            elif "c1ccccc1" in smi: # benzene
                tm_k, dh_fus = 278.7, 9870.0
            else:
                tm_k, dh_fus = np.nan, np.nan
        else:
            tm_k, dh_fus = np.nan, np.nan

        res[ik] = {
            "tm_K": tm_k,
            "dh_fus_J_mol": dh_fus,
            "n_carbon": nc,
            "cls": cls,
        }
    return res


def predict_freezing_point_coutinho(
    fractions: Sequence[float],
    inchikeys: Sequence[str],
    thermal_map: dict | None = None,
    basis: str = "mass",
) -> tuple[float, float, dict]:
    """Predict freezing point / wax appearance temperature using Coutinho non-ideal wax model.

    Returns:
        (T_freeze_C, Delta_T_eutectic, metadata_dict)
    """
    if thermal_map is None:
        thermal_map = get_pure_thermal_properties()

    w = np.asarray(fractions, float)
    w_sum = w.sum()
    if w_sum <= 0:
        return np.nan, 0.0, {}
    w = w / w_sum

    t_solutes = []
    n_carbons = []
    w_alkanes = []

    for wi, ik in zip(w, inchikeys):
        if ik not in thermal_map:
            continue
        info = thermal_map[ik]
        tm = info["tm_K"]
        dh = info["dh_fus_J_mol"]
        nc = info["n_carbon"]
        cls = info["cls"]

        if np.isfinite(tm) and np.isfinite(dh) and dh > 0 and wi > 1e-4:
            # Schroeder-van Laar ideal solubility equation:
            # 1/T = 1/Tm - (R/dH) * ln(x_i)
            inv_t = (1.0 / tm) - (R_GAS / dh) * math.log(wi)
            if inv_t > 0:
                t_solutes.append(1.0 / inv_t)
            if cls == "n_paraffin":
                n_carbons.append(nc)
                w_alkanes.append(wi)

    if not t_solutes:
        return -999.0, 0.0, {"reason": "no_crystallizing_components"}

    ideal_wat_K = max(t_solutes)
    ideal_wat_C = ideal_wat_K - 273.15

    # Carbon-number polydispersity of crystallizing n-alkanes
    if len(w_alkanes) >= 2:
        wa = np.array(w_alkanes) / sum(w_alkanes)
        na = np.array(n_carbons)
        mean_nc = np.sum(wa * na)
        std_nc = math.sqrt(float(np.sum(wa * (na - mean_nc) ** 2)))
    else:
        std_nc = 0.0

    paraffin_fraction = sum(w_alkanes)
    non_paraffin_fraction = 1.0 - paraffin_fraction

    # Multi-component Coutinho eutectic depression
    delta_eutectic = (
        COUTINHO_GAMMA
        + COUTINHO_ALPHA * min(std_nc, 4.0)
        + COUTINHO_BETA * non_paraffin_fraction
    )

    t_freeze_C = ideal_wat_C - delta_eutectic

    meta = {
        "ideal_wat_C": round(ideal_wat_C, 2),
        "delta_eutectic": round(delta_eutectic, 2),
        "paraffin_fraction": round(paraffin_fraction, 4),
        "paraffin_std_nc": round(std_nc, 2),
    }

    return round(t_freeze_C, 2), round(delta_eutectic, 2), meta


def ideal_wat_schroeder(
    fractions: Sequence[float],
    inchikeys: Sequence[str],
    thermal_map: dict | None = None,
) -> float:
    """Compute single-solute ideal Wax Appearance Temperature (WAT) via Schroeder-van Laar equation in °C."""
    _, _, meta = predict_freezing_point_coutinho(fractions, inchikeys, thermal_map)
    val = meta.get("ideal_wat_C", np.nan)
    return float(val) if val is not None else np.nan

