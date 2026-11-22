"""Leak-free ML pure-compound constants for the v5 blend routes (v8).

Measured constants always win; this layer only fills gaps and never reads measured values.

Property keys -> Track B pool targets (best method per `fuelblend.api.best_methods`):
    rho15        density            kg/m3   (D-MPNN long)
    Tm           tm_K               degC    (RF; pool is K, converted)
    dHfus        dhfus_kJ_mol       kJ/mol  (RF)
    logvisc293   logvisc_mPas_293   log10 mPa.s at 293 K (RF)
    flash_pure   flash_C            degC    (RF)

Temperature basis of the density label: results/v6/gnn/pool.csv is built by scripts/v6/gnn/build_pool.py from unified
rows with |T - 293.15 K| <= 5 K (median per compound), so the label is a 288 to 298 K value, NOT 288.15 K. No shift is
applied (a generic expansion slope would be an invented number); the shift is part of the substitution error that gate
G1 measures. Out-of-fold values come from the committed v7 Track B OOF files (scaffold 5-fold on the same fold ids).

Status: ml_rf_rdkit2d, ml_dmpnn_long, ml_out_of_domain, invalid_smiles, ml_model_missing.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from rdkit import Chem

from fuelblend.api import best_methods as BM

PROPS = {  # prop -> (pool target, unit, K->C offset)
    "rho15": ("density", "kg/m3", 0.0),
    "Tm": ("tm_K", "degC", -273.15),
    "dHfus": ("dhfus_kJ_mol", "kJ/mol", 0.0),
    "logvisc293": ("logvisc_mPas_293", "log10 mPa.s", 0.0),
    "flash_pure": ("flash_C", "degC", 0.0),
}
_OOF_COL = {"rf_rdkit2d": ("descriptor_gbm", "rf"), "dmpnn_long": ("curves_v1_long", "long")}
_STATUS = {"rf_rdkit2d": "ml_rf_rdkit2d", "dmpnn_long": "ml_dmpnn_long"}
OOD = "ml_out_of_domain"


def _ncarbon(mol) -> int:
    return sum(a.GetSymbol() == "C" for a in mol.GetAtoms())


def _is_hc(mol) -> bool:
    return all(a.GetSymbol() in ("C", "H") for a in mol.GetAtoms())


@lru_cache(maxsize=1)
def _pool():
    from fuelblend.v7 import runner
    comp, folds, _ = runner.load_track_b_pool()
    return comp, np.asarray(folds)


@lru_cache(maxsize=None)
def _carbon_range(target: str) -> tuple:
    comp, _ = _pool()
    sub = comp[comp[target].notna()]
    n = [_ncarbon(Chem.MolFromSmiles(s)) for s in sub.smiles]
    return (min(n), max(n))


def ml_constant(smiles: str, prop: str) -> tuple:
    """(value, status) from the deployed best-method model. Hydrocarbon, carbon count inside the training range."""
    target, _, off = PROPS[prop]
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return float("nan"), "invalid_smiles"
    lo, hi = _carbon_range(target)
    if not _is_hc(mol) or not lo <= _ncarbon(mol) <= hi:
        return float("nan"), OOD
    from fuelblend.api.predict_pure import predict_pure
    r = predict_pure(smiles, [target])[smiles][target]
    method = BM.best_method(target)["method"]
    if r["status"] != "ok":
        return float("nan"), "ml_model_missing"
    return float(r["value"]) + off, _STATUS[method]


def oof_constants(prop: str) -> pd.DataFrame:
    """Out-of-fold value for every pool compound: columns inchikey, smiles, value, fold, method."""
    target, _, off = PROPS[prop]
    method = BM.best_method(target)["method"]
    sub, suf = _OOF_COL[method]
    comp, folds = _pool()
    o = pd.read_csv(_root() / "results/v7/track_b" / sub / "oof.csv")
    if not (o.inchikey.values == comp.inchikey.values).all():
        raise ValueError("OOF rows are not aligned with the Track B pool")
    df = pd.DataFrame({"inchikey": comp.inchikey.values, "smiles": comp.smiles.values,
                       "value": o[f"{target}_{suf}"].values + off, "fold": folds, "method": method})
    return df[np.isfinite(df.value)].reset_index(drop=True)


def _root():
    from pathlib import Path
    return Path(__file__).resolve().parents[3]
