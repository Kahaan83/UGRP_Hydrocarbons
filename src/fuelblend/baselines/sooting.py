"""Sooting / smoke-point baselines (v4 Task 14).

* ysi_mix: linear-by-mole blending of the unified Yield Sooting Index (Das et al.).
* smoke_point_from_ysi(ysi, mw): TSI-form relation MW/SP = c0 + c1*YSI fitted on pure
  compounds having BOTH a measured (uncensored) ASTM D1322 smoke point and a YSI.
* smoke_point_from_bulk(aromatics_vol, diaromatics_vol): OLS on AJFTD jet fuels.

Censoring: ASTM D1322 lamp limit is 50 mm; rows with censored=True mean SP > 50 mm and are
NEVER used as exact targets. Predictions are capped at 50 mm ("at least 50").

Provenance (see scripts/v4/build_sooting_inputs.py): YSI = raw_file; pure-compound smoke points
= cited_typed (hand-typed from Mensch 2010 / Li & Sunderland 2012 / Boehm 2022, not checkable
from files in this repo); AJFTD = raw_file. The 175 extras tsi_table rows whose smoke point was
computed from YSI (derived=True) are circular and are not used.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import warnings

import numpy as np
import pandas as pd
from rdkit import Chem

DATA = Path(__file__).resolve().parents[3] / "data" / "processed" / "v4" / "sooting"
SP_CAP_MM = 50.0


def _canon(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    return Chem.MolToSmiles(m)


@lru_cache(maxsize=1)
def ysi_table() -> dict[str, float]:
    df = pd.read_csv(DATA / "ysi_hc.csv")
    return {_canon(s): float(y) for s, y in zip(df.smiles, df.ysi)}


def ysi_mix(smiles, x_mole, table: dict[str, float] | None = None) -> float:
    """Mole-fraction-weighted YSI. Raises KeyError for a component with no stored YSI."""
    tab = ysi_table() if table is None else table
    x = np.asarray(x_mole, float)
    if len(smiles) != len(x) or np.any(x < 0) or x.sum() <= 0:
        raise ValueError("bad composition")
    x = x / x.sum()
    vals = []
    for s in smiles:
        c = _canon(s)
        if c not in tab:
            raise KeyError(f"no YSI for {s}")
        vals.append(tab[c])
    return float(np.dot(x, vals))


def mixture_mw(smiles, x_mole) -> float:
    """Mole-averaged molecular weight (g/mol), for smoke_point_from_ysi on a mixture."""
    from rdkit.Chem import Descriptors
    x = np.asarray(x_mole, float)
    x = x / x.sum()
    return float(sum(xi * Descriptors.MolWt(Chem.MolFromSmiles(s)) for xi, s in zip(x, smiles)))


class OutOfDomainWarning(UserWarning):
    """Input outside the training domain of a smoke-point model."""


def exact_rows(df: pd.DataFrame, value_col: str = "smoke_point_mm") -> pd.DataFrame:
    """Rows usable as exact targets: not censored and strictly below the lamp limit."""
    ok = ~df["censored"].astype(bool) & (df[value_col] < SP_CAP_MM)
    return df[ok]


def fit_ysi_sp(df: pd.DataFrame) -> tuple[float, float]:
    """OLS of MW/SP on YSI using uncensored rows only. Returns (c0, c1)."""
    d = exact_rows(df).dropna(subset=["ysi", "mw", "smoke_point_mm"])
    if len(d) < 3:
        raise ValueError("need >=3 uncensored pure compounds with YSI and smoke point")
    y = (d.mw / d.smoke_point_mm).to_numpy()
    A = np.c_[np.ones(len(d)), d.ysi.to_numpy()]
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(c[0]), float(c[1])


def predict_ysi_sp(ysi, mw, coef) -> np.ndarray:
    ysi, mw = np.asarray(ysi, float), np.asarray(mw, float)
    den = coef[0] + coef[1] * ysi
    with np.errstate(divide="ignore", invalid="ignore"):
        sp = np.where(den > 0, mw / den, SP_CAP_MM)
    return np.clip(sp, 0.0, SP_CAP_MM)


@lru_cache(maxsize=1)
def _pure_coef() -> tuple[float, float]:
    return fit_ysi_sp(pd.read_csv(DATA / "smoke_point_pure_citedtyped.csv"))


def smoke_point_from_ysi(ysi: float, mw: float) -> float:
    return float(predict_ysi_sp(ysi, mw, _pure_coef()))


BULK_FEATURES = ("aromatics_pct_vol", "diaromatics_pct_vol")


def fit_bulk(df: pd.DataFrame, features=BULK_FEATURES, value_col="value"):
    """OLS of smoke point on bulk composition. Returns (coef incl. intercept, feature means)."""
    d = df[~df["censored"].astype(bool)] if "censored" in df else df
    d = d.dropna(subset=[value_col, features[0]])
    X = d[list(features)].copy()
    means = X.mean()
    X = X.fillna(means)
    A = np.c_[np.ones(len(d)), X.to_numpy()]
    c, *_ = np.linalg.lstsq(A, d[value_col].to_numpy(), rcond=None)
    return c, means


def predict_bulk(df: pd.DataFrame, coef, means, features=BULK_FEATURES) -> np.ndarray:
    X = df[list(features)].fillna(means)
    return np.clip(np.c_[np.ones(len(df)), X.to_numpy()] @ coef, 0.0, SP_CAP_MM)


@lru_cache(maxsize=1)
def _bulk_model():
    return fit_bulk(pd.read_csv(DATA / "smoke_point_ajftd.csv"))


def smoke_point_from_bulk(aromatics_vol: float, diaromatics_vol: float | None = None) -> float:
    """Smoke point (mm) of a jet fuel from D1319 aromatics and diaromatics (vol %).

    Trained on 100 AJFTD jet fuels (aromatics ~8-23 vol %, SP 19.3-25 mm); do not use outside.
    Emits OutOfDomainWarning if aromatics is outside [8, 23] vol % or the result is outside
    [19, 25] mm (trust range); the value is still returned. Clamped to [0, 50] mm."""
    if not (8.0 <= aromatics_vol <= 23.0):
        warnings.warn(f"aromatics {aromatics_vol} vol % outside training domain 8-23", OutOfDomainWarning)
    coef, means = _bulk_model()
    d = pd.DataFrame({BULK_FEATURES[0]: [aromatics_vol],
                      BULK_FEATURES[1]: [np.nan if diaromatics_vol is None else diaromatics_vol]})
    out = float(predict_bulk(d, coef, means)[0])
    if not (19.0 <= out <= 25.0):
        warnings.warn(f"predicted smoke point {out:.1f} mm outside trusted range 19-25 mm", OutOfDomainWarning)
    return out
