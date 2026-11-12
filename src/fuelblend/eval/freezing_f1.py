"""v5 API: F1x freezing recalibration (a + b * T_first, T_first from the complete Coutinho n-alkane table).

``fit_f1x`` is the same least squares the S5 ``F1x`` column uses (``np.polyfit`` of ``freezing_C`` on ``F0x`` = the first-species
temperature built with the complete Coutinho n-alkane table), fitted once over all S5 headline groups.  ``F1X_LOGO_MAE`` is the
pooled leave-one-group-out MAE read from BASE ``delta_freezing/delta.csv``: it is the only freezing score that may be quoted;
the in-sample fit error must not be.
"""
from __future__ import annotations

import importlib.util
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
S5_PREDICTIONS = ROOT / "results" / "v5" / "S5" / "s5_predictions.csv"
DELTA_FREEZING = ROOT / "results" / "v5" / "BASE" / "delta_freezing" / "delta.csv"
F1X_MODEL_NAME = "F1x_coutinho_complete"


def _read_logo_mae() -> float:
    d = pd.read_csv(DELTA_FREEZING)
    r = d[(d.property == "freezing_point_C") & (d.tier == "T2_measured") & (d.set == "pooled") & (d.model == F1X_MODEL_NAME)]
    if len(r) != 1 or int(r.iloc[0].n_groups) != 49:
        raise ValueError(f"{DELTA_FREEZING}: expected one pooled {F1X_MODEL_NAME} row with 49 groups")
    return float(r.iloc[0].mae)


F1X_LOGO_MAE: float = _read_logo_mae()


@lru_cache(maxsize=1)
def _s5():
    spec = importlib.util.spec_from_file_location("s5_freezing_nalkane_api", ROOT / "scripts" / "v5" / "s5_freezing_nalkane.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def t_first(moles_by_cn: dict, total_moles: float, consts: dict):
    """Re-export of S5 ``t_first``: (max_i T_i in K, carbon number of that species)."""
    return _s5().t_first(moles_by_cn, total_moles, consts)


def nalkane_mw(n: int) -> float:
    return _s5().nalkane_mw(n)


@lru_cache(maxsize=1)
def coutinho_consts() -> dict:
    """carbon number -> {tm_K, dh_J_mol} for C6-C36 n-alkanes (complete Coutinho table; same call as S5 ``coutinho_nalkane_constants``)."""
    from fuelblend.thermo.coutinho_wax import get_pure_thermal_properties
    s5 = _s5()
    ik = {n: s5.nalkane_ik(n) for n in range(6, 37)}
    info = pd.DataFrame({"smiles": ["C" * n for n in ik], "cls": "n_paraffin", "n_carbon": [float(n) for n in ik]}, index=list(ik.values()))
    t = get_pure_thermal_properties(info)
    return {n: dict(tm_K=t[k]["tm_K"], dh_J_mol=t[k]["dh_fus_J_mol"], tier="coutinho_table") for n, k in ik.items()}


@lru_cache(maxsize=1)
def fit_f1x() -> tuple[float, float]:
    """(a, b) of freezing_C = a + b * T_first_C over all S5 headline groups (one row per group) with a finite T_first."""
    d = pd.read_csv(S5_PREDICTIONS)
    h = d[(d.role == "headline") & d.F0x.notna() & d.freezing_C.notna()]
    if h.group.nunique() != len(h):
        raise ValueError("S5 headline rows are not one per group")
    b, a = np.polyfit(h.F0x.to_numpy(float), h.freezing_C.to_numpy(float), 1)
    return float(a), float(b)
