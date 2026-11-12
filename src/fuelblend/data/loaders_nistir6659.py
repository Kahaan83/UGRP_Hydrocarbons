"""Loader: NISTIR 6659 / Outcalt, Laesecke & Freund (Energy & Fuels 2009) compressed-liquid jet-fuel densities.

Input: the extracted table (columns fuel_id, fuel_name, category, T_K, P_MPa, P_kPa, density_kg_m3, uncertainty_kg_m3, source_table, reference).
Output columns: fuel_id, posf, property, value, T_K, P_kPa, source_file, provenance.  Units: kg/m3, K, kPa.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

COLUMNS = ["fuel_id", "posf", "property", "value", "T_K", "P_kPa", "source_file", "provenance"]


def _posf(fuel_id: str, name: str) -> str:
    m = re.search(r"POSF\s*(\d{3,5})", str(name)) or re.search(r"-(\d{4,5})$", str(fuel_id))
    return m.group(1) if m else ""


def load_nistir6659(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path)
    need = {"fuel_id", "fuel_name", "T_K", "P_MPa", "P_kPa", "density_kg_m3", "source_table", "reference"}
    if not need <= set(d.columns):
        raise ValueError(f"missing columns: {sorted(need - set(d.columns))}")
    if not np.allclose(d.P_kPa, d.P_MPa * 1000.0, rtol=1e-6, atol=0.6):
        raise ValueError("P_kPa is not 1000 * P_MPa")
    if not d.density_kg_m3.between(500, 1100).all() or not d.T_K.between(150, 700).all():
        raise ValueError("density (kg/m3) or T (K) outside the plausible liquid range: unit error?")
    out = pd.DataFrame(dict(fuel_id=d.fuel_id, posf=[_posf(f, n) for f, n in zip(d.fuel_id, d.fuel_name)], property="density_kg_m3",
                            value=d.density_kg_m3.astype(float), T_K=d.T_K.astype(float), P_kPa=d.P_kPa.astype(float),
                            source_file=Path(path).name, provenance=d.reference.astype(str) + " " + d.source_table.astype(str)))
    return out[COLUMNS]
