"""Fuel family classification for mixtures and real fuels (T1.4).

Classes:
  - 'jet': kerosene-range mixtures and aviation turbine fuels (typical C8-C16)
  - 'diesel': middle-distillate fuels, cetane reference blends, and diesel surrogates (typical C10-C22)
  - 'gasoline': spark-ignition fuels, octane reference blends (PRF/TRF), FACE gasolines (typical C4-C11)
  - 'model-mixture': laboratory binary/ternary/model mixtures without specific fuel-surrogate intent
"""
from __future__ import annotations

import json
from typing import Dict, List, Any
import pandas as pd


def classify_mixture(
    system_id: str,
    n_comp: int,
    comp_inchikeys: List[str],
    properties: List[str],
    comp_c_nums: List[int],
    source: str = "",
) -> Dict[str, Any]:
    """Classify a mixture into fuel_family and real_fuel boolean.

    Returns:
      {'fuel_family': str, 'real_fuel': bool, 'reason': str}
    """
    is_real = str(system_id).startswith("fuel:") or any(s in str(source).lower() for s in ["ajftd", "crc647", "crc_"])

    if is_real:
        # Check source or system annotations
        src_lower = str(source).lower()
        if "diesel" in src_lower:
            return {"fuel_family": "diesel", "real_fuel": True, "reason": "real diesel source"}
        if "gasoline" in src_lower or "face" in src_lower:
            return {"fuel_family": "gasoline", "real_fuel": True, "reason": "real gasoline source"}
        return {"fuel_family": "jet", "real_fuel": True, "reason": "aviation fuel source (CRC/AJFTD)"}

    # Speciated model or surrogate mixtures
    props_set = set(p.lower() for p in properties)
    avg_c = sum(comp_c_nums) / len(comp_c_nums) if comp_c_nums else 0.0

    # 1. Octane measurements or detailed gasoline surrogates
    if any(p in props_set for p in ["ron", "mon", "octane_sensitivity"]):
        return {"fuel_family": "gasoline", "real_fuel": False, "reason": "octane property (RON/MON)"}

    if n_comp >= 8 and avg_c <= 11.0:
        return {"fuel_family": "gasoline", "real_fuel": False, "reason": "multi-component gasoline surrogate"}

    # 2. Cetane / DCN measurements
    if "dcn" in props_set:
        if avg_c >= 10.0:
            return {"fuel_family": "diesel", "real_fuel": False, "reason": "cetane property with C>=10 range"}
        else:
            return {"fuel_family": "gasoline", "real_fuel": False, "reason": "cetane property with light C range"}

    # 3. Dedicated jet surrogates (e.g. flash-point multi-alkanes or ternary cold-flow in C10-C16)
    if "flash_point_C" in props_set and n_comp >= 3 and 10.0 <= avg_c <= 15.0:
        return {"fuel_family": "jet", "real_fuel": False, "reason": "multi-component kerosene-range flash blend"}

    # 4. Standard academic binary / ternary mixtures
    return {"fuel_family": "model-mixture", "real_fuel": False, "reason": f"academic {n_comp}-component mixture"}
