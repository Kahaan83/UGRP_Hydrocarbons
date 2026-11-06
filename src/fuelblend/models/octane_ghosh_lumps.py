"""Ghosh-Style Lump-Based Octane Model (G1-G4).

References:
    - Ghosh et al. (2006) Ind. Eng. Chem. Res. 45:337-345.
    - Anderson et al. (2012) SAE Int. J. Fuels Lubr. 5:1001-1014.

Groups gasoline components into chemical lumps:
    - n-paraffins (NP): n-heptane, light n-alkanes (RON_0 = 0.0, MON_0 = 0.0)
    - iso-paraffins (IP): isooctane, alkylates (RON_0 = 100.0, MON_0 = 100.0)
    - aromatics (A): toluene, reformate aromatics (RON_0 = 120.1, MON_0 = 103.5)
    - olefins (O): 1-hexene, FCC olefins (RON_0 = 76.4, MON_0 = 63.4)
    - naphthenes (N): cycloalkanes (RON_0 = 80.0, MON_0 = 75.0)
    - ethanol (EtOH): fuel ethanol (RON_0 = 108.6, MON_0 = 89.7)
    - unresolved tail (tail): heavy ends / trace DHA (RON_0 = 82.0, MON_0 = 78.0)

Incorporates binary class-interactions with L2 shrinkage regularization:
    - Ethanol-paraffin synergistic boost (k_EtOH_P)
    - Ethanol saturation curve (e * (1 - e))
    - Aromatic-olefin interaction (k_A_O)
    - Aromatic-paraffin interaction (k_A_P)
    - Olefin-paraffin interaction (k_O_P)
"""
from __future__ import annotations

from typing import Sequence, Dict, Any, Tuple
import numpy as np
import pandas as pd
from sklearn.linear_model import HuberRegressor, Ridge
from sklearn.model_selection import LeaveOneOut

PURE_LUMP_PROPERTIES = {
    "isooctane": {"ron": 100.0, "mon": 100.0, "mw": 114.23, "class": "iso_paraffin"},
    "n_heptane": {"ron": 0.0, "mon": 0.0, "mw": 100.20, "class": "n_paraffin"},
    "toluene": {"ron": 120.1, "mon": 103.5, "mw": 92.14, "class": "aromatic"},
    "olefin": {"ron": 76.4, "mon": 63.4, "mw": 84.16, "class": "alkene"},
    "naphthene": {"ron": 80.0, "mon": 75.0, "mw": 98.0, "class": "cycloalkane"},
    "ethanol": {"ron": 108.6, "mon": 89.7, "mw": 46.07, "class": "oxygenate"},
    "tail": {"ron": 82.0, "mon": 78.0, "mw": 110.0, "class": "other"},
}


def classify_component_to_lump(smiles: str | None, inchikey: str | None) -> str:
    """Map a component SMILES or InChIKey to one of the 7 Ghosh lumps."""
    if not isinstance(smiles, str) or pd.isna(smiles) or not smiles:
        if isinstance(inchikey, str) and "ETHANOL" in inchikey:
            return "ethanol"
        return "tail"

    smi = smiles.strip()
    if smi == "CCO" or smi == "CO":
        return "ethanol"
    elif smi == "CC(C)CC(C)(C)C":
        return "isooctane"
    elif smi == "CCCCCCC":
        return "n_heptane"
    elif smi == "Cc1ccccc1" or "c1ccccc1" in smi.lower():
        return "toluene"
    elif "=" in smi:
        return "olefin"
    elif "1" in smi and "c" not in smi:
        return "naphthene"
    elif "(" in smi:
        return "isooctane"
    elif len(smi) > 3:
        return "n_heptane"
    return "tail"


def extract_lump_profile(fuel_comp: pd.DataFrame) -> Dict[str, float]:
    """Aggregate composition fractions into normalized mass and mole lump fractions."""
    lump_mass = {k: 0.0 for k in PURE_LUMP_PROPERTIES}
    lump_moles = {k: 0.0 for k in PURE_LUMP_PROPERTIES}

    for _, row in fuel_comp.iterrows():
        lump = classify_component_to_lump(row.get("smiles"), row.get("inchikey"))
        w = float(row.get("fraction", 0.0))
        lump_mass[lump] += w
        mw = PURE_LUMP_PROPERTIES[lump]["mw"]
        lump_moles[lump] += w / mw

    w_tot = sum(lump_mass.values()) or 1.0
    n_tot = sum(lump_moles.values()) or 1.0

    profile = {}
    for k in PURE_LUMP_PROPERTIES:
        profile[f"w_{k}"] = lump_mass[k] / w_tot
        profile[f"x_{k}"] = lump_moles[k] / n_tot
    return profile


def compute_linear_baseline(profile: Dict[str, float], basis: str = "mole") -> Tuple[float, float]:
    """Calculate standard linear unadjusted octane blending baseline."""
    prefix = "x_" if basis == "mole" else "w_"
    ron = sum(profile[f"{prefix}{k}"] * PURE_LUMP_PROPERTIES[k]["ron"] for k in PURE_LUMP_PROPERTIES)
    mon = sum(profile[f"{prefix}{k}"] * PURE_LUMP_PROPERTIES[k]["mon"] for k in PURE_LUMP_PROPERTIES)
    return float(ron), float(mon)


def build_interaction_features(profile: Dict[str, float]) -> np.ndarray:
    """Build Ghosh interaction feature vector from lump fractions."""
    e = profile["w_ethanol"]
    tol = profile["w_toluene"]
    ole = profile["w_olefin"]
    iso = profile["w_isooctane"]
    hep = profile["w_n_heptane"]
    par = iso + hep

    # Non-linear synergistic terms:
    # 1. Ethanol-paraffin synergy
    # 2. Ethanol saturation curve: e * (1 - e)
    # 3. Aromatic-olefin antagonism
    # 4. Aromatic-paraffin boost
    # 5. Olefin-paraffin interaction
    # 6. Low-dilution oxygenate curvature
    return np.array([
        e * (1.0 - e),
        e * par,
        e * tol,
        tol * ole,
        tol * par,
        ole * par,
        e * (1.0 - np.exp(-5.0 * e)),
    ], dtype=float)


def compute_blending_octane_value(
    base_on: float,
    blend_on: float,
    blendstock_vol_fraction: float,
) -> float:
    """Compute Blending Octane Value (BOV) of a component in a base fuel.
    
    BOV = (ON_blend - (1 - v_B) * ON_base) / v_B
    """
    if blendstock_vol_fraction <= 0.0 or blendstock_vol_fraction >= 1.0:
        return blend_on
    return float((blend_on - (1.0 - blendstock_vol_fraction) * base_on) / blendstock_vol_fraction)
