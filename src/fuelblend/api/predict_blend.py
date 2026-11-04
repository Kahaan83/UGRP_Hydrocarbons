"""predict_blend.py — composition-to-properties wrapper (T2).

API
---
    from fuelblend.api.predict_blend import predict

    result = predict(
        composition,        # list of (inchikey|smiles, mass_fraction)
        T_K=298.15,         # temperature in Kelvin
        T_C=None,           # alternative: temperature in Celsius (explicit kwarg only)
        mole_fractions=None # if given: fractions in composition are mole fractions
    )

Returns a dict::

    {
      "viscosity_dyn_mPas": {"value": x, "std": s, "method": str, "warnings": [...]},
      "density_kg_m3":      {"value": x, "std": s, "method": str, "warnings": [...]},
      "flash_point_C":      {"value": x, "std": s, "interval": [lo, hi],
                             "method": str, "warnings": [...]},
      "T_K": float,
      "warnings": [global_warnings],
    }

Property paths (documented in reports/T2_WRAPPER.md)
------------------------------------------------------
viscosity: pure GBM ensemble → Grunberg–Nissan baseline → + mixture-residual GBM correction
density  : pure GBM ensemble → ideal-volume baseline    → + mixture-residual GBM correction
flash    : Le Chatelier + Raoult (Antoine, LFL) if available; else weighted pure GBM; else None
"""
from __future__ import annotations

import json
import math
import sys
import warnings as _warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import joblib

# ── path setup (scripts do the same) ─────────────────────────────────────────
_BASE = Path(__file__).resolve().parents[3]   # fuelblend_data/
sys.path.insert(0, str(_BASE / "src"))

from fuelblend.features.pure_features import (    # noqa: E402
    pure_feature_row, FEATURE_NAMES, CLASSES,
)
from fuelblend.features.mixture_features import extract_mixture_features  # noqa: E402

# ── constants ────────────────────────────────────────────────────────────────
_SEEDS = [11, 23, 37, 41, 59]
_PURE_MODEL_DIR = _BASE / "runs" / "pure_models"
_MIX_MODEL_DIR  = _BASE / "runs" / "mixture_models"
_COMPOUNDS_CSV  = _BASE / "data" / "unified" / "compounds.csv"
_MODEL_TABLE    = _BASE / "data" / "processed" / "model_table.csv"
_BASELINE_PREDS = _BASE / "data" / "runs" / "baseline_preds.csv"

T_DOMAIN_LO, T_DOMAIN_HI = 230.0, 470.0   # K
P_TOTAL = 101.325                           # kPa
FLASH_LO_K, FLASH_HI_K = -150.0 + 273.15, 400.0 + 273.15

# split-conformal α=0.10 safety quantile on S1 flash over-prediction residuals
# computed once here from baseline_preds.csv (see _load_flash_conformal())
# stored as module-level to avoid repeated file reads
_FLASH_CONF_Q: float | None = None

# ── lazy-loaded caches ────────────────────────────────────────────────────────
_compounds: pd.DataFrame | None = None
_smiles_to_ik: dict[str, str] | None = None
_pure_models: dict | None = None      # prop -> list of 5 fitted models
_mix_models:  dict | None = None      # prop -> list of 5 fitted models (or None if unavailable)
_pure_fits:   dict | None = None      # {"eta": {ik: {a,b,lo,hi}}, "rho": {ik: ...}}
_compounds_dict: dict | None = None   # {ik: {mw, n_carbon, cls}}


# ── RDKit import (optional; degrades gracefully) ──────────────────────────────
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchi
    RDLogger.DisableLog("rdApp.*")
    _RDKIT = True
except ImportError:
    _RDKIT = False


# =============================================================================
# Data loading helpers
# =============================================================================

def _load_compounds() -> pd.DataFrame:
    global _compounds, _smiles_to_ik
    if _compounds is None:
        _compounds = pd.read_csv(_COMPOUNDS_CSV, low_memory=False).set_index("inchikey")
        _smiles_to_ik = {}
        for ik, smi in zip(_compounds.index, _compounds["smiles"]):
            if pd.notna(smi):
                _smiles_to_ik[str(smi).strip()] = str(ik)
    return _compounds


def _load_compounds_dict() -> dict:
    global _compounds_dict
    if _compounds_dict is None:
        c = _load_compounds()
        _compounds_dict = {}
        for ik, row in c.iterrows():
            _compounds_dict[str(ik)] = {
                "mw": float(row["mw"]) if pd.notna(row.get("mw")) else 100.0,
                "n_carbon": float(row["n_carbon"]) if pd.notna(row.get("n_carbon")) else 7.0,
                "cls": str(row["cls"]) if pd.notna(row.get("cls")) else "other",
            }
    return _compounds_dict


def _load_pure_models() -> dict:
    global _pure_models
    if _pure_models is None:
        _pure_models = {}
        for prop in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C", "vapor_pressure_kPa"):
            models = []
            for s in _SEEDS:
                p = _PURE_MODEL_DIR / f"{prop}__seed_{s}.joblib"
                if p.exists():
                    models.append(joblib.load(p))
            _pure_models[prop] = models if len(models) == 5 else []
    return _pure_models


def _load_mix_models() -> dict:
    """Load mixture-residual models; train+save them on first use if missing."""
    global _mix_models
    if _mix_models is not None:
        return _mix_models

    _mix_models = {}
    for prop in ("viscosity_dyn_mPas", "density_kg_m3"):
        models = []
        for s in _SEEDS:
            p = _MIX_MODEL_DIR / f"{prop}__seed_{s}.joblib"
            if p.exists():
                models.append(joblib.load(p))
        if len(models) == 5:
            _mix_models[prop] = models
        else:
            _mix_models[prop] = None   # will be trained on demand

    # Train any that are missing
    _ensure_mix_models()
    return _mix_models


def _ensure_mix_models():
    """Lazily train and save mixture-residual models if joblibs are absent.

    Same hyperparams as train_mixture_residual.py (150 iters, seeds [11,23,37,41,59]).
    Join preds ↔ model_table on (system_id, property, T_K rounded to 4dp) to recover
    comp_inchikeys / w_mass / x_mole / quality_weight (not stored in baseline_preds.csv).
    """
    from sklearn.ensemble import HistGradientBoostingRegressor

    need = [p for p, m in _mix_models.items() if m is None]
    if not need:
        return

    print(f"[predict_blend] Mixture-residual models not found for {need}. Training now (< 1 min)...")

    mt = pd.read_csv(_MODEL_TABLE, low_memory=False)
    preds_df = pd.read_csv(_BASELINE_PREDS, low_memory=False)
    locked = json.loads((_BASE / "splits" / "locked.json").read_text())
    lk_sys = set(locked.get("multi_system_ids_locked", []))

    # Build a lookup table from mt keyed by (system_id, property, T_K rounded to 4dp)
    mt_clean = mt[~mt["system_id"].isin(lk_sys)].copy()
    mt_clean["T_K_r4"] = mt_clean["T_K"].round(4).astype(str)
    mt_lookup = (
        mt_clean[["system_id", "property", "T_K_r4",
                  "comp_inchikeys", "w_mass", "x_mole", "quality_weight"]]
        .drop_duplicates(subset=["system_id", "property", "T_K_r4"])
    )

    preds_df = preds_df.copy()
    preds_df["T_K_r4"] = pd.to_numeric(preds_df["T_K"], errors="coerce").round(4).astype(str)

    preds_merged = preds_df.merge(
        mt_lookup, on=["system_id", "property", "T_K_r4"], how="inner",
    )
    n_matched = len(preds_merged)
    n_total = len(preds_df)

    if n_matched < n_total * 0.8:
        print(f"[predict_blend] Low join coverage ({n_matched}/{n_total}); "
              "falling back to no mixture-residual correction.")
        for p in need:
            _mix_models[p] = []
        return
    if n_matched < n_total:
        print(f"[predict_blend] Joined {n_matched}/{n_total} pred rows to model_table.")
    preds_df = preds_merged

    cd = _load_compounds_dict()
    _MIX_MODEL_DIR.mkdir(parents=True, exist_ok=True)

    cfgs = {
        "viscosity_dyn_mPas": dict(baseline="gn_mass", residual="log",
                                   max_depth=6, target_log=True),
        "density_kg_m3":      dict(baseline="dens_idvol", residual="linear",
                                   max_depth=5, target_log=False),
    }

    for prop in need:
        cfg = cfgs[prop]
        base_col = cfg["baseline"]
        mask = (
            (preds_df["property"] == prop)
            & (preds_df["conflict"] == False)
            & preds_df[base_col].notna()
            & (preds_df[base_col] > 0 if cfg["target_log"] else True)
            & (preds_df["value"] > 0 if cfg["target_log"] else True)
            & preds_df["comp_inchikeys"].notna()
        )
        data = preds_df[mask].copy()

        feats = [extract_mixture_features(r, cd) for _, r in data.iterrows()]
        X = pd.DataFrame(feats)

        if cfg["target_log"]:
            residual = np.log(data["value"].values) - np.log(data[base_col].values)
        else:
            residual = data["value"].values - data[base_col].values

        train_mask = (data["split"] == "S1").values
        models = []
        for s in _SEEDS:
            gbr = HistGradientBoostingRegressor(
                max_iter=150, max_depth=cfg["max_depth"],
                min_samples_leaf=20, learning_rate=0.05, random_state=s,
            )
            sw = np.nan_to_num(data["quality_weight"].values[train_mask], nan=1.0)
            gbr.fit(X[train_mask], residual[train_mask], sample_weight=sw)
            joblib.dump(gbr, _MIX_MODEL_DIR / f"{prop}__seed_{s}.joblib")
            models.append(gbr)

        _mix_models[prop] = models
        s2 = (data["split"] == "S2").values
        if s2.sum() > 0:
            s2_preds = np.mean([m.predict(X[s2]) for m in models], axis=0)
            rmse_base = float(np.sqrt(np.mean(residual[s2] ** 2)))
            rmse_mod  = float(np.sqrt(np.mean((residual[s2] - s2_preds) ** 2)))
            print(f"  [{prop}] S2 residual RMSE: {rmse_base:.4f} -> {rmse_mod:.4f}")
        print(f"  [{prop}] {len(models)} models saved to {_MIX_MODEL_DIR}")


def _load_pure_fits() -> dict:
    """Build per-compound ln-eta and rho-linear fits from model_table pure rows.

    Replicates baselines_eval.build_fits exactly (same _fit_loglin/_fit_linear).
    """
    global _pure_fits
    if _pure_fits is not None:
        return _pure_fits

    mt = pd.read_csv(_MODEL_TABLE, low_memory=False)
    pure = mt[mt["kind"] == "pure"].copy()

    def _jload(s):
        try:
            return json.loads(s) if isinstance(s, str) else []
        except Exception:
            return []

    pure["inchikey"] = pure["comp_inchikeys"].map(lambda s: (_jload(s) or [""])[0])

    def _fit_loglin(T, y):
        T, ly = np.asarray(T, float), np.log(np.asarray(y, float))
        lt = np.log(T)
        if len(T) == 1:
            return dict(a=float(ly[0]), b=0.0, lo=float(T[0]), hi=float(T[0]))
        if len(T) >= 3:
            i, j = np.triu_indices(len(T), 1)
            dlt = lt[j] - lt[i]
            ok = np.abs(dlt) > 1e-9
            b = float(np.median((ly[j] - ly[i])[ok] / dlt[ok])) if ok.any() else 0.0
            a = float(np.median(ly - b * lt))
            keep = np.abs(ly - (a + b * lt)) < np.log(3.0)
            if keep.sum() >= 2:
                T, ly, lt = T[keep], ly[keep], lt[keep]
        b, a = np.polyfit(lt, ly, 1)
        return dict(a=float(a), b=float(b), lo=float(T.min()), hi=float(T.max()))

    def _fit_linear(T, y):
        T, y = np.asarray(T, float), np.asarray(y, float)
        if len(T) == 1:
            return dict(a=float(y[0]), b=0.0, lo=float(T[0]), hi=float(T[0]))
        b, a = np.polyfit(T, y, 1)
        return dict(a=float(a), b=float(b), lo=float(T.min()), hi=float(T.max()))

    fe, fr = {}, {}
    p = pure[(pure["property"] == "viscosity_dyn_mPas") & (pure["value"] > 0)]
    for ik, g in p.groupby("inchikey"):
        g = g.groupby(g["T_K"].round(3))["value"].median()
        fe[str(ik)] = _fit_loglin(g.index.values, g.values)

    p2 = pure[pure["property"] == "density_kg_m3"]
    for ik, g in p2.groupby("inchikey"):
        g = g.groupby(g["T_K"].round(3))["value"].median()
        fr[str(ik)] = _fit_linear(g.index.values, g.values)

    _pure_fits = {"eta": fe, "rho": fr}
    return _pure_fits


def _load_flash_conformal() -> float:
    """Split-conformal safety quantile: 90th percentile of S1 over-prediction residuals.

    over-prediction residual = chat_lflT - measured  (positive = model says safer than reality)
    α = 0.10 → 1-α quantile = 90th percentile.
    Returns the quantile value in °C (positive = conservative upper bound offset).
    """
    global _FLASH_CONF_Q
    if _FLASH_CONF_Q is not None:
        return _FLASH_CONF_Q

    df = pd.read_csv(_BASELINE_PREDS, low_memory=False)
    s1 = df[(df["property"] == "flash_point_C") & (df["split"] == "S1") & df["chat_lflT"].notna()]
    if len(s1) == 0:
        _FLASH_CONF_Q = float("nan")
        return _FLASH_CONF_Q
    # residual = pred - measured; positive = over-prediction (dangerous direction)
    res = s1["chat_lflT"].values - s1["value"].values
    _FLASH_CONF_Q = float(np.quantile(res, 0.90))  # 90th pct → 90% coverage
    return _FLASH_CONF_Q


# =============================================================================
# Input resolution: smiles / inchikey → inchikey + compounds.csv lookup
# =============================================================================

def _smiles_to_inchikey(smiles: str) -> str | None:
    if not _RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    inchi = MolToInchi(mol)
    if inchi is None:
        return None
    from rdkit.Chem.inchi import InchiToInchiKey
    return InchiToInchiKey(inchi)


def _resolve_components(composition, mole_fractions) -> tuple[list, list, list[str]]:
    """Returns (inchikeys, w_mass_norm, global_warnings)."""
    comp_df = _load_compounds()
    smiles_to_ik = _smiles_to_ik if _smiles_to_ik is not None else {}

    global_warn: list[str] = []
    inchikeys, fracs = [], []

    for token, frac in composition:
        token = str(token).strip()
        # 1. direct inchikey lookup
        if token in comp_df.index:
            ik = token
        # 2. SMILES in compounds.csv
        elif token in smiles_to_ik:
            ik = smiles_to_ik[token]
        # 3. Try RDKit canonicalise
        else:
            ik = _smiles_to_inchikey(token) if _RDKIT else None
            if ik and ik in comp_df.index:
                pass  # found via inchikey from smiles
            else:
                if ik:  # we got an inchikey but it's not in our registry
                    global_warn.append(f"unknown_compound: {token!r} resolved to {ik} but not in compounds registry; using partial data")
                else:
                    global_warn.append(f"unknown_compound: {token!r} could not be resolved; predictions may degrade")
                    ik = token   # keep as-is; downstream code will handle None lookups
        inchikeys.append(ik)
        fracs.append(float(frac))

    # convert mole → mass if requested
    if mole_fractions is not None:
        mole_fracs = np.array(mole_fractions, dtype=float)
        mws = []
        for ik in inchikeys:
            mw = float(comp_df.loc[ik, "mw"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "mw"]) else 100.0
            mws.append(mw)
        mws = np.array(mws, dtype=float)
        masses = mole_fracs * mws
        total = masses.sum()
        if total <= 0:
            global_warn.append("fraction_conversion_failed: could not convert mole→mass")
            w_norm = np.ones(len(inchikeys)) / len(inchikeys)
        else:
            w_norm = masses / total
    else:
        w_arr = np.array(fracs, dtype=float)
        total = w_arr.sum()
        if abs(total - 1.0) > 1e-4:
            global_warn.append(f"fraction_normalized: sum was {total:.6f}")
        if total > 0:
            w_norm = w_arr / total
        else:
            w_norm = np.ones(len(inchikeys)) / len(inchikeys)
            global_warn.append("fraction_error: all fractions zero")

    return inchikeys, w_norm.tolist(), global_warn


# =============================================================================
# Pure-GBM ensemble prediction for a single component
# =============================================================================

def _pure_gbm_predict(prop: str, ik: str, smiles: str, cls_val: str, T_K: float) -> tuple[float | None, float | None]:
    """Returns (mean, std) over seeds in the model's native space (ln for viscosity, etc.)."""
    models = _load_pure_models().get(prop, [])
    if not models:
        return None, None
    feat = pure_feature_row(ik, smiles, T_K, cls_val)
    preds = []
    for m in models:
        try:
            p = float(m.predict(feat.reshape(1, -1))[0])
            if math.isfinite(p):
                preds.append(p)
        except Exception:
            pass
    if not preds:
        return None, None
    return float(np.mean(preds)), float(np.std(preds))


def _pred_value(prop: str, raw: float) -> float:
    """Invert the training transform (matches train_pure_models.py)."""
    if prop == "viscosity_dyn_mPas":
        return float(np.exp(raw))
    if prop == "vapor_pressure_kPa":
        return float(10.0 ** raw)
    return float(raw)


# =============================================================================
# Domain checks
# =============================================================================

def _check_T_domain(T_K: float) -> list[str]:
    if not (T_DOMAIN_LO <= T_K <= T_DOMAIN_HI):
        return [f"T_out_of_domain: T={T_K:.1f} K outside [{T_DOMAIN_LO}, {T_DOMAIN_HI}] K training range"]
    return []


def _check_pure_extrapolation(ik: str, T_K: float, prop: str) -> list[str]:
    """Check if the component has pure-model data within ±20 K of T."""
    fits = _load_pure_fits()
    fmap = fits["eta"] if "viscosity" in prop else fits["rho"]
    f = fmap.get(ik)
    if f is None:
        return [f"pure_extrapolation: no pure-data curve for {ik}"]
    if not (f["lo"] - 20 <= T_K <= f["hi"] + 20):
        return [f"pure_extrapolation: {ik} pure data [{f['lo']:.0f},{f['hi']:.0f}] K; T={T_K:.0f} K outside ±20 K"]
    return []


# =============================================================================
# Physics baselines (replicating baselines_eval.py formulas exactly)
# =============================================================================

_T_PAD = 2.0  # K of allowed extrapolation outside pure-data range (matches baselines_eval.py)


def _eta_at(ik: str, T: float) -> float | None:
    f = _load_pure_fits()["eta"].get(ik)
    if f is None or not (f["lo"] - _T_PAD <= T <= f["hi"] + _T_PAD):
        return None
    return float(np.exp(f["a"] + f["b"] * np.log(T)))


def _rho_at(ik: str, T: float) -> float | None:
    f = _load_pure_fits()["rho"].get(ik)
    if f is None or not (f["lo"] - _T_PAD <= T <= f["hi"] + _T_PAD):
        return None
    v = float(f["a"] + f["b"] * T)
    return v if v > 0 else None


def _gn_mass(inchikeys, w_mass, T_K) -> float | None:
    """Grunberg–Nissan mass-weighted: exp(sum w_i ln eta_i)."""
    vals = []
    for ik, wi in zip(inchikeys, w_mass):
        e = _eta_at(ik, T_K)
        if e is None or e <= 0:
            return None
        vals.append(wi * math.log(e))
    return float(math.exp(sum(vals)))


def _refutas(inchikeys, w_mass, T_K) -> tuple[float | None, float | None, float | None]:
    """Returns (eta_dyn_mPas, nu_cSt, rho_mix_kg_m3)."""
    vb, wr, rhos = 0.0, 0.0, []
    for ik, wi in zip(inchikeys, w_mass):
        e = _eta_at(ik, T_K)
        if e is None or e <= 0:
            return None, None, None
        r = _rho_at(ik, T_K)
        if r is None:
            return None, None, None
        nu = e / (r / 1000.0)  # mPa.s / (kg/m3 /1000) = mm2/s = cSt
        if nu < 0.2:
            return None, None, None
        vb += wi * (14.534 * math.log(math.log(nu + 0.8)) + 10.975)
        wr += wi
        rhos.append(r)
    if wr <= 0:
        return None, None, None
    vbn = vb / wr
    nu_mix = float(math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8)
    rho_mix = 1.0 / sum(wi / r for wi, r in zip(w_mass, rhos))
    eta_dyn = nu_mix * rho_mix / 1000.0
    return eta_dyn, nu_mix, rho_mix


def _dens_idvol(inchikeys, w_mass, T_K) -> float | None:
    """Ideal-volume density: 1 / sum(w_i/rho_i)."""
    inv = 0.0
    for ik, wi in zip(inchikeys, w_mass):
        r = _rho_at(ik, T_K)
        if r is None:
            return None
        inv += wi / r
    return (1.0 / inv) if inv > 0 else None


# ── Antoine + Le Chatelier flash ─────────────────────────────────────────────

def _psat_kpa(A, B, C, T_K: float) -> float:
    d = T_K + C
    if not math.isfinite(d) or d <= 1.0:
        return math.nan
    return float(10.0 ** (A - B / d))


def _mole_from_mass(inchikeys, w_mass) -> list[float] | None:
    comp_df = _load_compounds()
    num = []
    for ik, wi in zip(inchikeys, w_mass):
        if ik not in comp_df.index or pd.isna(comp_df.loc[ik, "mw"]):
            return None
        num.append(float(wi) / float(comp_df.loc[ik, "mw"]))
    s = sum(num)
    return [v / s for v in num] if s > 0 else None


def _flash_lechatelier(inchikeys, x_mole, warnings_out: list[str], lfl_temp=True) -> float | None:
    """Le Chatelier + Raoult. Returns flash point °C or None."""
    from scipy.optimize import brentq

    comp_df = _load_compounds()
    A, B, C, lfl = [], [], [], []

    for ik in inchikeys:
        if ik not in comp_df.index:
            warnings_out.append(f"no_antoine: {ik} not in registry")
            return None
        row = comp_df.loc[ik]
        if pd.isna(row.get("antoine_A")) or pd.isna(row.get("antoine_B")) or pd.isna(row.get("antoine_C")):
            warnings_out.append(f"no_antoine: {ik} missing Antoine constants")
            return None
        if pd.isna(row.get("lfl_vol_pct")):
            warnings_out.append(f"no_lfl: {ik} missing LFL")
            return None
        # Antoine window check
        Tlo = float(row["antoine_Tlo"]) if pd.notna(row.get("antoine_Tlo")) else 0.0
        Thi = float(row["antoine_Thi"]) if pd.notna(row.get("antoine_Thi")) else 9999.0
        if FLASH_LO_K > Thi + 5 or FLASH_HI_K < Tlo - 5:
            warnings_out.append(f"antoine_out_of_range: {ik} Antoine window [{Tlo:.0f},{Thi:.0f}] K "
                                 f"does not overlap flash search range")
            return None
        A.append(float(row["antoine_A"]))
        B.append(float(row["antoine_B"]))
        C.append(float(row["antoine_C"]))
        lfl.append(float(row["lfl_vol_pct"]) / 100.0)

    x = np.asarray(x_mole, float)
    if x.sum() <= 0:
        return None

    def haz(T_K: float) -> float:
        tc = T_K - 273.15
        fac = max(1.0 - 0.000721 * (tc - 25.0), 0.05) if lfl_temp else 1.0
        h = 0.0
        for xi, a, b, c, l in zip(x, A, B, C, lfl):
            p = _psat_kpa(a, b, c, T_K)
            if not math.isfinite(p):
                return math.nan
            h += xi * p / (P_TOTAL * l * fac)
        return h - 1.0

    lo, hi = FLASH_LO_K, FLASH_HI_K
    try:
        flo, fhi = haz(lo), haz(hi)
        if not math.isfinite(flo) or flo > 0:
            warnings_out.append("no_root_in_bracket: H(T_lo)>0 or NaN for flash Le Chatelier")
            return None
        if not math.isfinite(fhi) or fhi < 0:
            warnings_out.append("no_root_in_bracket: H(T_hi)<0 or NaN for flash Le Chatelier")
            return None
        T_flash = brentq(haz, lo, hi, xtol=1e-6, rtol=1e-10)
        return float(T_flash - 273.15)
    except Exception as e:
        warnings_out.append(f"solver_error: {e}")
        return None


# =============================================================================
# Mixture-residual GBM correction
# =============================================================================

def _make_mix_row(inchikeys, w_mass, x_mole, T_K, n_comp):
    """Build a pseudo-row for extract_mixture_features."""
    return {
        "comp_inchikeys": json.dumps(inchikeys),
        "w_mass": json.dumps(w_mass),
        "x_mole": json.dumps(x_mole),
        "T_K": T_K,
    }


def _mix_residual_predict(prop: str, mix_row: dict) -> tuple[float | None, float | None]:
    """Returns (mean_residual, std_residual) from mixture GBM ensemble, or (None, None)."""
    mix_models = _load_mix_models()
    models = mix_models.get(prop)
    if not models:
        return None, None
    cd = _load_compounds_dict()
    feats = extract_mixture_features(mix_row, cd)
    if not feats:
        return None, None
    X = pd.DataFrame([feats])
    preds = []
    for m in models:
        try:
            preds.append(float(m.predict(X)[0]))
        except Exception:
            pass
    if not preds:
        return None, None
    return float(np.mean(preds)), float(np.std(preds))


# =============================================================================
# Per-property prediction paths
# =============================================================================

def _predict_viscosity(inchikeys, w_mass, x_mole, T_K) -> dict:
    warns = list(_check_T_domain(T_K))
    for ik in inchikeys:
        warns.extend(_check_pure_extrapolation(ik, T_K, "viscosity_dyn_mPas"))

    # ── Physics baseline: try Refutas, fall back to GN mass ──────────────────
    eta_dyn, nu_mix, rho_mix = _refutas(inchikeys, w_mass, T_K)
    if eta_dyn is not None:
        phys_val = eta_dyn
        phys_method = "refutas"
    else:
        phys_val = _gn_mass(inchikeys, w_mass, T_K)
        phys_method = "gn_mass"

    # ── Fall back to pure-GBM ensemble average if physics fails ──────────────
    if phys_val is None:
        # Try pure GBMs per component, then log-linear mix
        comp_df = _load_compounds()
        comp_preds, comp_stds = [], []
        for ik, wi in zip(inchikeys, w_mass):
            smi = str(comp_df.loc[ik, "smiles"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "smiles"]) else ""
            cls_val = str(comp_df.loc[ik, "cls"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "cls"]) else "other"
            mean_raw, std_raw = _pure_gbm_predict("viscosity_dyn_mPas", ik, smi, cls_val, T_K)
            if mean_raw is not None:
                comp_preds.append((wi, mean_raw, std_raw or 0.0))
            else:
                warns.append(f"no_pure_model: could not predict viscosity for {ik}")
        if comp_preds:
            ln_mix = sum(wi * mn for wi, mn, _ in comp_preds)
            var_mix = sum((wi * st) ** 2 for wi, _, st in comp_preds)
            phys_val = float(np.exp(ln_mix))
            phys_std = float(phys_val * np.sqrt(var_mix))  # delta method
            phys_method = "gn_fallback"
            return {"value": phys_val, "std": phys_std, "method": phys_method, "warnings": warns}
        else:
            warns.append("mixture_size_warning") if len(inchikeys) > 5 else None
            return {"value": None, "std": None, "method": "none", "warnings": warns}

    # ── Mixture-residual GBM correction ──────────────────────────────────────
    mix_row = _make_mix_row(inchikeys, w_mass, x_mole, T_K, len(inchikeys))
    res_mean, res_std = _mix_residual_predict("viscosity_dyn_mPas", mix_row)

    if res_mean is not None:
        # residual is in ln-space; apply: ln(eta_corrected) = ln(phys) + residual
        ln_phys = math.log(phys_val)
        ln_corrected = ln_phys + res_mean
        val = float(math.exp(ln_corrected))
        # std: propagate residual ensemble std through exp
        std = float(math.exp(ln_corrected) * abs(res_std)) if res_std else None
        method = f"residual_gbm+{phys_method}"
    else:
        val = phys_val
        std = None
        method = phys_method

    if len(inchikeys) > 5:
        warns.append("mixture_size_warning: n_comp>5 is outside residual model training range")

    return {"value": val, "std": std, "method": method, "warnings": warns}


def _predict_density(inchikeys, w_mass, x_mole, T_K) -> dict:
    warns = list(_check_T_domain(T_K))
    for ik in inchikeys:
        warns.extend(_check_pure_extrapolation(ik, T_K, "density_kg_m3"))

    phys_val = _dens_idvol(inchikeys, w_mass, T_K)

    if phys_val is None:
        # Fallback: weighted average of pure GBM predictions
        comp_df = _load_compounds()
        wt_vals, wt_stds = [], []
        for ik, wi in zip(inchikeys, w_mass):
            smi = str(comp_df.loc[ik, "smiles"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "smiles"]) else ""
            cls_val = str(comp_df.loc[ik, "cls"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "cls"]) else "other"
            mean_raw, std_raw = _pure_gbm_predict("density_kg_m3", ik, smi, cls_val, T_K)
            if mean_raw is not None:
                wt_vals.append((wi, mean_raw))
                wt_stds.append((wi, std_raw or 0.0))
            else:
                warns.append(f"no_pure_model: could not predict density for {ik}")
        if wt_vals:
            val = float(sum(wi * v for wi, v in wt_vals))
            std = float(math.sqrt(sum((wi * s) ** 2 for wi, s in wt_stds)))
            return {"value": val, "std": std, "method": "gn_fallback", "warnings": warns}
        return {"value": None, "std": None, "method": "none", "warnings": warns}

    # Mixture-residual correction (density residual ≈ 0.7% S2 improvement — applied but noted)
    mix_row = _make_mix_row(inchikeys, w_mass, x_mole, T_K, len(inchikeys))
    res_mean, res_std = _mix_residual_predict("density_kg_m3", mix_row)

    if res_mean is not None:
        val = phys_val + res_mean   # residual in kg/m3 (linear)
        std = float(abs(res_std)) if res_std else None
        method = "residual_gbm+dens_idvol"
    else:
        val = phys_val
        std = None
        method = "dens_idvol"

    if len(inchikeys) > 5:
        warns.append("mixture_size_warning: n_comp>5 is outside residual model training range")

    return {"value": val, "std": std, "method": method, "warnings": warns}


def _predict_flash(inchikeys, w_mass, x_mole, T_K) -> dict:
    """Flash-point prediction path with Le Chatelier fallback to pure-GBM weighted average."""
    warns = list(_check_T_domain(T_K))
    comp_df = _load_compounds()

    # Check non-HC
    for ik in inchikeys:
        if ik in comp_df.index and not bool(comp_df.loc[ik, "is_hc"]):
            warns.append(f"non_hc_component: {ik} is not a hydrocarbon; flash prediction unreliable")

    flash_warns: list[str] = []
    val = _flash_lechatelier(inchikeys, x_mole, flash_warns, lfl_temp=True)

    if val is not None:
        method = "chat"
        warns.extend(flash_warns)
        warns.append("flash_data_thin: mixture flash training is sparse (123 binary + 3 multi rows)")
        # Conformal interval
        q = _load_flash_conformal()
        interval = [val - abs(q), val + abs(q)] if math.isfinite(q) else [None, None]
        return {"value": val, "std": None, "interval": interval, "method": method, "warnings": warns}

    warns.extend(flash_warns)

    # Fallback: weighted pure-GBM flash estimates
    wt_vals, wt_stds = [], []
    for ik, wi in zip(inchikeys, w_mass):
        smi = str(comp_df.loc[ik, "smiles"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "smiles"]) else ""
        cls_val = str(comp_df.loc[ik, "cls"]) if ik in comp_df.index and pd.notna(comp_df.loc[ik, "cls"]) else "other"
        mean_raw, std_raw = _pure_gbm_predict("flash_point_C", ik, smi, cls_val, T_K)
        if mean_raw is not None:
            wt_vals.append((wi, float(mean_raw)))
            wt_stds.append((wi, float(std_raw) if std_raw else 0.0))

    if wt_vals:
        gbm_val = float(sum(wi * v for wi, v in wt_vals))
        gbm_std = float(math.sqrt(sum((wi * s) ** 2 for wi, s in wt_stds)))
        warns.append("flash_data_thin: mixture flash training is sparse (123 binary + 3 multi rows)")
        warns.append("flash_fallback_gbm_mix: Le Chatelier unavailable (missing Antoine/LFL); "
                     "using mass-weighted average of pure-GBM flash — mixture non-ideality ignored")
        return {"value": gbm_val, "std": gbm_std, "interval": [None, None], "method": "gbm_mix", "warnings": warns}

    warns.append("flash_none: cannot predict flash point (no Antoine/LFL and no pure GBM coverage)")
    return {"value": None, "std": None, "interval": [None, None], "method": "none", "warnings": warns}


# =============================================================================
# Public API
# =============================================================================

def predict(
    composition,
    T_K: float = 298.15,
    T_C: float | None = None,
    mole_fractions=None,
) -> dict[str, Any]:
    """Predict viscosity, density, and flash point for a hydrocarbon blend.

    Parameters
    ----------
    composition : list of (str, float)
        Each element is (inchikey|smiles, mass_fraction). Fractions are auto-
        normalized if their sum deviates from 1.0 by < 1e-4. Pass exactly
        mass fractions; for mole fractions use `mole_fractions=`.
    T_K : float
        Temperature in Kelvin (default 298.15). Use T_C= for Celsius.
    T_C : float, optional
        Temperature in Celsius. If given, T_K is ignored.
    mole_fractions : list of float, optional
        If given, the second element of each composition tuple is ignored and
        these mole fractions (aligned to composition order) are converted to
        mass fractions using MW from compounds.csv.

    Returns
    -------
    dict with keys:
        viscosity_dyn_mPas, density_kg_m3, flash_point_C  — each a dict with
            value (float|None), std (float|None), method (str), warnings (list)
        flash_point_C additionally has interval ([lo,hi] or [None,None])
        T_K (float), warnings (list of global warnings)
    """
    if T_C is not None:
        T_K = float(T_C) + 273.15
    T_K = float(T_K)

    inchikeys, w_mass, global_warn = _resolve_components(composition, mole_fractions)

    # Non-HC check at global level
    comp_df = _load_compounds()
    for ik in inchikeys:
        if ik in comp_df.index and not bool(comp_df.loc[ik, "is_hc"]):
            global_warn.append(f"non_hc_input: {ik} is not a hydrocarbon — model trained on HC only")

    x_mole = _mole_from_mass(inchikeys, w_mass)
    if x_mole is None:
        x_mole = w_mass   # fallback (equal mole/mass only for equal MW — wrong but non-crashing)
        global_warn.append("mole_fraction_approx: could not compute exact mole fractions from MW; using mass fractions")

    visc  = _predict_viscosity(inchikeys, w_mass, x_mole, T_K)
    dens  = _predict_density(inchikeys, w_mass, x_mole, T_K)
    flash = _predict_flash(inchikeys, w_mass, x_mole, T_K)

    return {
        "viscosity_dyn_mPas": visc,
        "density_kg_m3": dens,
        "flash_point_C": flash,
        "T_K": T_K,
        "warnings": global_warn,
    }


# =============================================================================
# CLI demo: three realistic examples
# =============================================================================

def _demo():
    import pprint

    print("=" * 70)
    print("EXAMPLE (a): n-decane / n-dodecane 50/50 wt% at T=298.15 K")
    print("Searching model_table for matching measured row...")
    print("=" * 70)

    dec_ik = "DIOQZVSQGTUSAI-UHFFFAOYSA-N"   # n-decane
    dod_ik = "SNRUBQQJIBEYMU-UHFFFAOYSA-N"   # n-dodecane
    dec_smi = "CCCCCCCCCC"
    dod_smi = "CCCCCCCCCCCC"

    # find closest measured row at 298.15K
    try:
        mt = pd.read_csv(_MODEL_TABLE, low_memory=False)
        v = mt[(mt["kind"] == "binary") & (mt["property"] == "viscosity_dyn_mPas")]
        mask = v["comp_inchikeys"].apply(lambda s: dec_ik in str(s) and dod_ik in str(s))
        cands = v[mask & (abs(v["T_K"] - 298.15) < 1.0)]
        if len(cands):
            # find closest to 50/50 by mass
            best = None
            best_diff = 99.0
            for _, row in cands.iterrows():
                ws = json.loads(row["w_mass"])
                iks = json.loads(row["comp_inchikeys"])
                dec_idx = iks.index(dec_ik)
                diff = abs(ws[dec_idx] - 0.5)
                if diff < best_diff:
                    best_diff, best = diff, row
            if best is not None:
                ws = json.loads(best["w_mass"])
                iks = json.loads(best["comp_inchikeys"])
                dec_idx = iks.index(dec_ik)
                print(f"  Closest measured row: w_decane={ws[dec_idx]:.3f}, T={best.T_K:.1f} K, "
                      f"viscosity={best['value']:.4f} mPa.s")
    except Exception as e:
        print(f"  (could not load measured data: {e})")

    comp_a = [(dec_smi, 0.5), (dod_smi, 0.5)]
    r = predict(comp_a, T_K=298.15)
    print(f"\n  Predicted viscosity: {r['viscosity_dyn_mPas']['value']:.4f} mPa.s "
          f"(method={r['viscosity_dyn_mPas']['method']})")
    print(f"  Predicted density:   {r['density_kg_m3']['value']:.2f} kg/m3")
    fp = r["flash_point_C"]
    print(f"  Flash point:         {fp['value']:.1f} °C (method={fp['method']}, "
          f"interval={[f'{v:.1f}' if v else None for v in fp['interval']]})")
    if r["warnings"]:
        print(f"  Global warnings: {r['warnings']}")

    print()
    print("=" * 70)
    print("EXAMPLE (b): 3-component blend: n-tetradecane/n-decane/methane (wt: 57.8/41.6/0.6%)")
    print("(A real model_table row; methane triggers a low-viscosity warning)")
    print("=" * 70)
    tet_smi = "CCCCCCCCCCCCCC"     # n-tetradecane
    met_smi = "C"                  # methane
    comp_b = [(tet_smi, 0.578), (dec_smi, 0.416), (met_smi, 0.006)]
    r2 = predict(comp_b, T_K=294.4)
    print(f"  Predicted viscosity: {r2['viscosity_dyn_mPas']['value']:.4f} mPa.s "
          f"(method={r2['viscosity_dyn_mPas']['method']})")
    print(f"  Predicted density:   {r2['density_kg_m3']['value']:.2f} kg/m3")
    fp2 = r2["flash_point_C"]
    print(f"  Flash point:         {fp2['value']} °C (method={fp2['method']})")
    if r2["warnings"]:
        print(f"  Global warnings: {r2['warnings']}")
    for prop in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C"):
        if r2[prop]["warnings"]:
            print(f"  {prop} warnings: {r2[prop]['warnings'][:2]}")

    print()
    print("=" * 70)
    print("EXAMPLE (c): Unknown compound + n-decane — fallback path")
    print("=" * 70)
    comp_c = [(dec_smi, 0.7), ("CCCOOOO[Zn]CCCC", 0.3)]
    r3 = predict(comp_c, T_K=310.0)
    print(f"  Global warnings: {r3['warnings']}")
    print(f"  Viscosity method: {r3['viscosity_dyn_mPas']['method']}")
    print(f"  Flash method:     {r3['flash_point_C']['method']}")
    pprint.pprint({k: v for k, v in r3.items() if k not in ("warnings",)})


if __name__ == "__main__":
    _demo()
