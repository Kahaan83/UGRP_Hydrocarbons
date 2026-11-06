"""Activity-Coefficient DVPE Model for Gasoline and Oxygenated Blends (G5-G6).

Governing Physics:
    1. ASTM D5191 Standard Vapor Pressure Correlation:
       DVPE (kPa) = 0.965 * P_tot(37.8 °C) - 3.78 kPa
    2. Non-ideal Ethanol Azeotropy (Margules Activity Model):
       ln(gamma_EtOH) = A_EtOH * (1 - x_EtOH)^2
       where A_EtOH in [1.5, 2.2] (calibrated default = 1.80).
       At low dilution (x_EtOH -> 0), gamma_EtOH ~ 4.5 to 6.0, generating
       the characteristic +7 to +10 kPa vapor pressure jump in E10-E20 blends.
    3. Light-End Speciation:
       Incorporates C4-C5 butane/pentane volatility which carries >50% of finished
       gasoline vapor pressure.
"""
from __future__ import annotations

from typing import Dict, Any, Tuple, Sequence
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneOut

T_TEST_K = 310.95  # 37.8 °C (100 °F) ASTM D5191 test temperature

# Saturation vapor pressures at 37.8 °C [kPa]
PSAT_37_8C_KPA = {
    "ethanol": 17.5,
    "methanol": 35.0,
    "isooctane": 7.0,
    "n_heptane": 6.2,
    "toluene": 6.8,
    "1_hexene": 38.0,
    "light_ends_c4_c5": 140.0,  # Effective butane/isopentane mixture
    "heavy_ends": 5.0,
}


def calculate_ethanol_activity_coefficient(x_ethanol: float, A_param: float = 1.80) -> float:
    """Margules 1-parameter activity coefficient for ethanol in hydrocarbon matrix."""
    if x_ethanol <= 0.0 or x_ethanol >= 1.0:
        return 1.0
    return float(np.exp(A_param * ((1.0 - x_ethanol) ** 2)))


def predict_dvpe_ideal_raoult(fractions: Sequence[float], smiles_list: Sequence[str]) -> float:
    """Calculate baseline DVPE using zero-parameter ideal Raoult's law (no azeotrope/light-end correction)."""
    p_partials = []
    w_sum = sum(fractions) or 1.0
    for w, smi in zip(fractions, smiles_list):
        if smi in ("CCO", "CO"):
            psat = PSAT_37_8C_KPA["ethanol"]
        elif smi == "CC(C)CC(C)(C)C":
            psat = PSAT_37_8C_KPA["isooctane"]
        elif smi == "CCCCCCC":
            psat = PSAT_37_8C_KPA["n_heptane"]
        elif smi == "Cc1ccccc1":
            psat = PSAT_37_8C_KPA["toluene"]
        elif smi == "CCCCC=C":
            psat = PSAT_37_8C_KPA["1_hexene"]
        else:
            psat = PSAT_37_8C_KPA["heavy_ends"]
        p_partials.append((w / w_sum) * psat)

    p_tot = sum(p_partials)
    return float(max(0.0, 0.965 * p_tot - 3.78))


def extract_dvpe_features(fuel_comp: pd.DataFrame, A_param: float = 1.80) -> Dict[str, float]:
    """Extract composition and non-ideal activity features for DVPE modeling."""
    w_eth = 0.0
    w_iso = 0.0
    w_hep = 0.0
    w_tol = 0.0
    w_ole = 0.0
    w_tail = 0.0

    for _, row in fuel_comp.iterrows():
        smi = str(row.get("smiles", ""))
        w = float(row.get("fraction", 0.0))
        if smi in ("CCO", "CO") or "ETHANOL" in str(row.get("inchikey", "")):
            w_eth += w
        elif smi == "CC(C)CC(C)(C)C":
            w_iso += w
        elif smi == "CCCCCCC":
            w_hep += w
        elif smi == "Cc1ccccc1" or "c1ccccc1" in smi.lower():
            w_tol += w
        elif "=" in smi:
            w_ole += w
        else:
            w_tail += w

    # Mole fraction of ethanol (MW_eth = 46.07 vs average MW_HC ~ 100)
    n_eth = w_eth / 46.07
    n_hc = (1.0 - w_eth) / 100.0 if (1.0 - w_eth) > 0 else 0.0
    x_eth = n_eth / (n_eth + n_hc) if (n_eth + n_hc) > 0 else 0.0

    gamma_eth = calculate_ethanol_activity_coefficient(x_eth, A_param)
    p_eth_partial = x_eth * gamma_eth * PSAT_37_8C_KPA["ethanol"]

    return {
        "w_eth": w_eth,
        "x_eth": x_eth,
        "gamma_eth": gamma_eth,
        "p_eth_partial": p_eth_partial,
        "w_iso": w_iso,
        "w_hep": w_hep,
        "w_tol": w_tol,
        "w_ole": w_ole,
        "w_tail": w_tail,
    }
