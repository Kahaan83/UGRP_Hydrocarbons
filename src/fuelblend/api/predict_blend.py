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
      "flash_point_C":      {"value": x, "std": s, "method": str, "warnings": [...]},
      # every property dict also carries (v4): lo, hi, interval, model_name, in_domain, margin, domain_reasons, tuple
      "uncertified": {freezing_point_C|smoke_point_mm|ron|mon|dcn|dvpe_kpa|pour_point_C: no-margin placeholders},
      "T_K": float,
      "warnings": [global_warnings],
    }

Property paths (documented in reports/T2_WRAPPER.md)
------------------------------------------------------
v4 routing (Task 18). Intervals come ONLY from configs/uncertainty_v4.yaml `applied_margin`; 'inf' => (nan, nan).
v5: predict() delegates to ``fuelblend.api.predict_route.predict_fuel`` (the BASE scoring route): pure data from the frozen
provider grid (never rebuilt), a component with no own data for a property is replaced by a nearest same-kind representative
for that property only, non-hydrocarbons are unmapped, and margins apply only when nothing was substituted or unmapped.
composition_kind="known" (default): viscosity raw Refutas (any n_comp) with the lab_mixtures.raw_refutas margin of the
           temperature bin; density ideal volume with side-specific lab_mixtures.ideal_volume over/under margins.
           Hydrocarbons only; a bin without a valid margin is not certified.
composition_kind="reconstructed": viscosity recalibrated Refutas (-40..20 C; raw Refutas outside, uncertified);
           density margin density_kg_m3.crc647_ideal_volume.all_T; fuel_type='jet' required.
           Viscosity and density are None (status T_out_of_grid) outside the 230-470 K grid.
flash    : ideal Le Chatelier + Raoult on the anchored Antoine / LFL registry with the LFL overlay (U1), independent of T_K;
           margin = max(lab ideal, D93 ideal) applied_margin (lab 3.70 C, set by one system: warning lab_margin_set_by_one_system);
           the v4 routing value is variants['v4_designated'] on predict_fuel. There is no fallback model: if flash fails after
           substitution the value is None and the method 'none'.
The mixture-residual GBM is opt-in (use_mixture_residual=True) and has no margin. No model is ever trained at load:
a missing bundle raises FileNotFoundError naming the command that builds it.
Margins were validated on CRC 647 jet fuels only: pass fuel_type='jet' or in_domain is False.
Extra predict() kwargs: fuel_type (None|'jet'|...; 'jet' is TRUSTED from the caller, the composition is not checked
against a jet specification) and use_mixture_residual (default False; opt-in GBM correction, no margin, needs the bundle).
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

T_DOMAIN_LO, T_DOMAIN_HI = 230.0, 470.0   # K
VISC_RECAL_T_C = (-40.0, 20.0)       # CRC 647 window the Refutas recalibration + its yaml margin were validated on
DENSITY_MARGIN_T_C = (-0.2, 69.5)    # CRC 647 density residual temperature range (yaml scope)
CERTIFIED_FUEL_KIND = "jet"          # every yaml margin was validated on CRC 647 / AJFTD jet fuels only
COMPOSITION_KINDS = ("known", "reconstructed")

# v4 routing: (property, method) -> yaml margin keys combined with max(). Mirrors src/fuelblend/design/models_v4.py.
MARGIN_KEYS = {
    ("viscosity_dyn_mPas", "refutas_recalibrated"): ["viscosity_kin_cSt.crc647.recalibrated_refutas.all_T"],
    ("density_kg_m3", "dens_idvol"): ["density_kg_m3.crc647_ideal_volume.all_T"],
    # yaml recommendation: NOT the 0-floored D93 margin alone; use the larger of the lab and D93 applied margins
    ("flash_point_C", "ideal_le_chatelier_anchored"): ["flash_point_C.lab_mixtures.ideal", "flash_point_C.real_fuels_D93_reconstructed.ideal"],
}
# composition_kind="known": raw Refutas / ideal volume + LAB-mixture margins (hydrocarbon-only, not jet specific).
# Viscosity margins are binned by temperature (same bins as engine_v4.viscosity_t_bin); density margins are side-specific:
# lo = value - over, hi = value + under.
LAB_VISC_KEY = "viscosity_kin_cSt.lab_mixtures.raw_refutas.{bin}"
LAB_DENSITY_KEYS = {"over": "density_kg_m3.lab_mixtures.ideal_volume.over",
                    "under": "density_kg_m3.lab_mixtures.ideal_volume.under"}
# properties with no valid margin: reported as (nan, nan) and in_domain False; never given a made-up number
UNCERTIFIED_PROPERTIES = {
    "freezing_point_C": "yaml freezing margins are 'inf' / not carried for blends (no valid applied_margin)",
    "smoke_point_mm": "yaml smoke_point margin is 'inf' (no valid applied_margin)",
    "ron": "octane margins exist only for base gasoline + one blendstock; no multi-component model is certified",
    "mon": "octane margins exist only for base gasoline + one blendstock; no multi-component model is certified",
    "dcn": "no diesel residuals exist; no margin",
    "dvpe_kpa": "DVPE data units are unconfirmed; excluded from uncertainty_v4.yaml (no margin)",
    "pour_point_C": "no v4 task evaluated pour point (no margin)",
}

_BUILD_MIX_CMD = ".venv/Scripts/python.exe extras/fuelblend_data/scripts/train_mixture_residual.py   (run from the repo root)"
_BUILD_PURE_CMD = ".venv/Scripts/python.exe extras/fuelblend_data/scripts/train_pure_models.py   (run from the repo root)"
_UNCERTAINTY_YAML = _BASE / "configs" / "uncertainty_v4.yaml"
_VISC_PRED = _BASE / "results" / "v4" / "V" / "viscosity_predictions_v4.csv"

# ── lazy-loaded caches ────────────────────────────────────────────────────────
_compounds: pd.DataFrame | None = None
_smiles_to_ik: dict[str, str] | None = None
_pure_models: dict | None = None      # prop -> list of 5 fitted models
_mix_models:  dict | None = None      # prop -> list of 5 fitted models (or None if unavailable)
_pure_fits:   dict | None = None      # {"eta": {ik: {a,b,lo,hi}}, "rho": {ik: ...}}
_compounds_dict: dict | None = None   # {ik: {mw, n_carbon, cls}}
_recal_b: tuple | None = None
_margins_cache: dict | None = None


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
            models, missing = [], []
            for s in _SEEDS:
                p = _PURE_MODEL_DIR / f"{prop}__seed_{s}.joblib"
                if p.exists():
                    models.append(joblib.load(p))
                else:
                    missing.append(p)
            if missing:
                _pure_models = None
                raise _bundle_missing("pure-component", missing, _BUILD_PURE_CMD)
            _pure_models[prop] = models
    return _pure_models


def _bundle_missing(kind: str, paths: list, build_cmd: str) -> FileNotFoundError:
    shown = ", ".join(sorted({p.parent.as_posix() for p in paths}))
    return FileNotFoundError(
        f"predict_blend: {kind} model bundle missing ({len(paths)} file(s) not found under {shown}). "
        f"predict_blend never trains models at load time. Build the bundle with: {build_cmd}"
    )


def _load_mix_models() -> dict:
    """Load the mixture-residual joblibs. NEVER trains: a missing file raises FileNotFoundError (v4 Task 18)."""
    global _mix_models
    if _mix_models is not None:
        return _mix_models
    missing, loaded = [], {}
    for prop in ("viscosity_dyn_mPas", "density_kg_m3"):
        models = []
        for s in _SEEDS:
            p = _MIX_MODEL_DIR / f"{prop}__seed_{s}.joblib"
            if p.exists():
                models.append(joblib.load(p))
            else:
                missing.append(p)
        loaded[prop] = models
    if missing:
        raise _bundle_missing("mixture-residual", missing, _BUILD_MIX_CMD)
    _mix_models = loaded
    return _mix_models


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


def _mole_from_mass(inchikeys, w_mass) -> list[float] | None:
    comp_df = _load_compounds()
    num = []
    for ik, wi in zip(inchikeys, w_mass):
        if ik not in comp_df.index or pd.isna(comp_df.loc[ik, "mw"]):
            return None
        num.append(float(wi) / float(comp_df.loc[ik, "mw"]))
    s = sum(num)
    return [v / s for v in num] if s > 0 else None


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

def _predict_viscosity(inchikeys, w_mass, x_mole, T_K, use_residual: bool = False,
                       composition_kind: str = "known") -> dict:
    warns = list(_check_T_domain(T_K))
    for ik in inchikeys:
        warns.extend(_check_pure_extrapolation(ik, T_K, "viscosity_dyn_mPas"))

    # ── known composition: raw Refutas for any number of components (a recalibration that switched on at n=2 made the
    #    result discontinuous between a pure compound and a 99.9/0.1 blend).  reconstructed composition (CRC fuels):
    #    recalibrated Refutas inside its -40..20 C validation window; raw Refutas outside it.
    eta_dyn, nu_mix, rho_mix = _refutas(inchikeys, w_mass, T_K)
    if eta_dyn is not None:
        lo_k, hi_k = VISC_RECAL_T_C[0] + 273.15, VISC_RECAL_T_C[1] + 273.15
        if composition_kind == "known":
            phys_val = eta_dyn
            phys_method = "refutas_raw"
        elif lo_k <= T_K <= hi_k:
            b0, b1 = _recal_beta()
            phys_val = eta_dyn * math.exp(b0 + b1 * 1000.0 / T_K)
            phys_method = "refutas_recalibrated"
        else:
            phys_val = eta_dyn
            phys_method = "refutas"
            warns.append(f"recalibration_not_applied: T={T_K - 273.15:.1f} C outside the CRC 647 recalibration window "
                         f"{VISC_RECAL_T_C} C; raw Refutas returned, no certified margin")
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

    # ── Optional (opt-in) mixture-residual GBM correction: no v4 margin exists for it ─────────────
    res_mean, res_std = None, None
    if use_residual:
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

    if len(inchikeys) > 5 and use_residual:
        warns.append("mixture_size_warning: n_comp>5 is outside residual model training range")

    return {"value": val, "std": std, "method": method, "warnings": warns}


def _predict_density(inchikeys, w_mass, x_mole, T_K, use_residual: bool = False) -> dict:
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

    # Optional (opt-in) mixture-residual correction: the v4 margin is for ideal volume mixing only
    res_mean, res_std = None, None
    if use_residual:
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

    if len(inchikeys) > 5 and use_residual:
        warns.append("mixture_size_warning: n_comp>5 is outside residual model training range")

    return {"value": val, "std": std, "method": method, "warnings": warns}


def _predict_flash(inchikeys, w_mass, x_mole, T_K) -> dict:
    """Flash point on the BASE scoring route (``predict_route``); returns the legacy {value, std, method, warnings}.
    No fallback: a flash that cannot be computed after substitution is value None / method 'none' with the reason."""
    from fuelblend.api import predict_route
    r = predict_route.predict_fuel(list(zip(inchikeys, w_mass)), T_K=T_K)["flash_point_C"]
    return {"value": r["value"], "std": None, "method": r["method"], "warnings": list(r["warnings"])}


# =============================================================================
# v4 recalibration + yaml margins (applied_margin ONLY)
# =============================================================================

def _recal_beta() -> tuple[float, float]:
    """ln(nu_true/nu_Refutas) = b0 + b1*1000/T, least squares on CRC 647 DEV fuels (same fit as models_v4.Library)."""
    global _recal_b
    if _recal_b is None:
        if not _VISC_PRED.exists():
            raise FileNotFoundError(
                f"predict_blend: {_VISC_PRED} missing (needed for the Refutas recalibration). Build it with: "
                ".venv/Scripts/python.exe scripts/v4/build_viscosity_leaderboard_v4.py")
        d = pd.read_csv(_VISC_PRED)
        d = d[~d.locked.astype(bool)]
        resid = np.log(d.y.to_numpy()) - np.log(d.nu_refutas.to_numpy())
        B = np.column_stack([np.ones(len(d)), 1000.0 / d.T_K.to_numpy()])
        b = np.linalg.lstsq(B, resid, rcond=None)[0]
        _recal_b = (float(b[0]), float(b[1]))
    return _recal_b


def _load_margins() -> dict:
    global _margins_cache
    if _margins_cache is None:
        import yaml
        if not _UNCERTAINTY_YAML.exists():
            raise FileNotFoundError(
                f"predict_blend: {_UNCERTAINTY_YAML} missing. Build it with: "
                ".venv/Scripts/python.exe scripts/v4/build_uncertainty_v4.py")
        with open(_UNCERTAINTY_YAML, "r", encoding="utf-8") as fh:
            _margins_cache = yaml.safe_load(fh)["margins"]
    return _margins_cache


def _applied_margin(keys, margins: dict | None = None):
    """Combine yaml entries by max(applied_margin). -> (margin|None, unit, direction, reason).

    Reads ONLY ``applied_margin``. An entry with no usable applied_margin (value 'inf') makes the whole result None:
    no number is ever substituted.
    """
    margins = _load_margins() if margins is None else margins
    vals, units, dirs = [], set(), set()
    for k in keys:
        e = margins.get(k)
        if e is None:
            return None, "", None, f"{k}: no such entry in uncertainty_v4.yaml"
        m = e.get("applied_margin")
        if m is None or isinstance(m, str) or not math.isfinite(float(m)):
            return None, "", None, f"{k}: no valid margin (applied_margin absent / inf)"
        vals.append(float(m))
        units.add(e.get("unit"))
        dirs.add(e.get("direction"))
    if not vals:
        return None, "", None, "no margin key for this model"
    unit = "ln-ratio" if "ln-ratio" in units else next(iter(units))
    direction = "over" if "over" in dirs else next(iter(dirs))
    return max(vals), unit, direction, ""


def interval_from_margin(value, margin, unit) -> tuple[float, float]:
    """value -/+ margin (additive), or value*exp(-/+ margin) for ln-ratio margins. NaN pair when no margin."""
    if value is None or margin is None or not math.isfinite(margin):
        return math.nan, math.nan
    if unit == "ln-ratio":
        return float(value) * math.exp(-margin), float(value) * math.exp(margin)
    return float(value) - margin, float(value) + margin


def _lab_scope(prop: str, method: str, composition_kind: str) -> bool:
    """True when the lab-mixture margins (hydrocarbon-only, not jet specific, not CRC-window limited) are the ones in use."""
    return (composition_kind == "known" and (prop, method) in (("viscosity_dyn_mPas", "refutas_raw"),
                                                               ("density_kg_m3", "dens_idvol")))


def _domain_reasons(prop: str, inchikeys, T_K: float, fuel_type, lab: bool = False) -> list[str]:
    """Reasons the yaml margin does NOT apply (empty list -> in_domain True).

    ``lab`` = the lab-mixture margins are in use: they are not jet specific and not tied to the CRC temperature windows
    (viscosity is binned by temperature; the bin's own margin decides), so those two checks are skipped."""
    from fuelblend.mixture.common import CLASSES as _CRC
    why = []
    if lab:
        pass
    elif fuel_type is None:
        why.append("fuel_type not declared: margins were validated on CRC 647 jet fuels only (pass fuel_type='jet' if it is one)")
    elif str(fuel_type).lower() != CERTIFIED_FUEL_KIND:
        why.append(f"fuel_type={fuel_type!r}: margins validated on jet fuels only")
    comp = _load_compounds()
    for ik in inchikeys:
        if ik not in comp.index:
            why.append(f"{ik}: not in compounds registry")
            continue
        if not bool(comp.loc[ik, "is_hc"]):
            why.append(f"{ik}: not a hydrocarbon")
        c = comp.loc[ik, "cls"]
        if pd.isna(c) or str(c) not in _CRC:
            why.append(f"{ik}: class {c!r} outside the validated classes {_CRC}")
    tc = T_K - 273.15
    if prop == "viscosity_dyn_mPas" and not lab and not (VISC_RECAL_T_C[0] <= tc <= VISC_RECAL_T_C[1]):
        why.append(f"T={tc:.1f} C outside viscosity validation window {VISC_RECAL_T_C} C")
    if prop == "density_kg_m3" and not lab and not (DENSITY_MARGIN_T_C[0] <= tc <= DENSITY_MARGIN_T_C[1]):
        why.append(f"T={tc:.1f} C outside density validation range {DENSITY_MARGIN_T_C} C")
    return why


def _lab_non_hc(inchikeys) -> list[str]:
    comp = _load_compounds()
    return [ik for ik in inchikeys if ik not in comp.index or not bool(comp.loc[ik, "is_hc"])]


def _attach_uncertainty(prop: str, res: dict, inchikeys, T_K: float, fuel_type, composition_kind: str = "reconstructed") -> dict:
    """Add lo/hi/model_name/in_domain/tuple to a property result using uncertainty_v4.yaml applied_margin only."""
    lab = _lab_scope(prop, res["method"], composition_kind)
    sided = None
    if lab:
        bad = _lab_non_hc(inchikeys)
        if bad:
            margin, unit, direction = None, "", None
            reason = f"lab-mixture margins are hydrocarbon-only; non-hydrocarbon or unregistered component(s) {bad}"
        elif prop == "viscosity_dyn_mPas":
            from fuelblend.design.engine_v4 import viscosity_t_bin
            margin, unit, direction, reason = _applied_margin([LAB_VISC_KEY.format(bin=viscosity_t_bin(T_K - 273.15))])
        else:                                   # density: side-specific
            m_over, unit, _, r1 = _applied_margin([LAB_DENSITY_KEYS["over"]])
            m_under, _, _, r2 = _applied_margin([LAB_DENSITY_KEYS["under"]])
            if m_over is None or m_under is None:
                margin, unit, direction, reason = None, "", None, "; ".join(r for r in (r1, r2) if r)
            else:
                margin, direction, reason, sided = max(m_over, m_under), "sided", "", (m_over, m_under)
    else:
        keys = MARGIN_KEYS.get((prop, res["method"]))
        if keys is None:
            margin, unit, direction = None, "", None
            reason = f"model {res['method']!r} has no uncertainty_v4.yaml margin"
        else:
            margin, unit, direction, reason = _applied_margin(keys)
    dom = _domain_reasons(prop, inchikeys, T_K, fuel_type, lab=lab) if res.get("value") is not None else ["no value"]
    if margin is None:
        dom = dom + [reason]            # uncertified: in_domain must be False
    if sided is not None:                # lo = value - over-prediction margin, hi = value + under-prediction margin
        lo, hi = float(res["value"]) - sided[0], float(res["value"]) + sided[1]
        res.update(margin_over=sided[0], margin_under=sided[1])
    else:
        lo, hi = interval_from_margin(res.get("value"), margin, unit)
    res.update(lo=lo, hi=hi, interval=[lo, hi], model_name=res["method"], in_domain=(not dom),
               margin=margin, margin_unit=unit, margin_direction=direction, domain_reasons=dom)
    if direction == "over":
        res["safe_bound"] = lo          # over-prediction is the dangerous direction for flash: lo is the safe side
    res["tuple"] = (res.get("value"), lo, hi, res["method"], res["in_domain"])
    return res


# =============================================================================
# Public API
# =============================================================================

def predict(
    composition,
    T_K: float = 298.15,
    T_C: float | None = None,
    mole_fractions=None,
    fuel_type: str | None = None,
    use_mixture_residual: bool = False,
    composition_kind: str = "known",
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
    fuel_type : str, optional
        Declared fuel kind. Only 'jet' can yield in_domain=True (margins validated on CRC 647 jet fuels). The value
        is trusted from the caller; the composition is not checked to be a jet fuel. None -> in_domain False.
    use_mixture_residual : bool
        Default False. True adds the mixture-residual GBM correction (needs the joblib bundle, else FileNotFoundError);
        it has no uncertainty_v4.yaml margin, so the result is not certified.
    composition_kind : {"known", "reconstructed"}
        "known" (default): the composition is a list of known compounds (a lab blend or a design). Viscosity is raw
        Refutas for any number of components and density is ideal volume mixing; both are judged by the hydrocarbon
        lab-mixture margins (viscosity binned by temperature, density side-specific), which are not jet specific.
        A temperature bin without a valid margin (below 10 C at present) is not certified.
        Cold limit: viscosity and density come from the frozen 230-470 K pure-property grid, so the cold path is
        raw Refutas as well; a bin without a valid margin (below 283 K at present) is in_domain=False. Outside the
        grid viscosity and density are None with status T_out_of_grid (flash does not depend on T).
        "reconstructed": the composition was reconstructed from bulk fuel properties (CRC 647 style); the CRC
        recalibrated-Refutas path and its CRC margins are used and fuel_type='jet' is still required for in_domain.
    mole_fractions : list of float, optional
        If given, the second element of each composition tuple is ignored and
        these mole fractions (aligned to composition order) are converted to
        mass fractions using MW from compounds.csv.

    Returns
    -------
    dict with keys:
        viscosity_dyn_mPas, density_kg_m3, flash_point_C  — each a dict with
            value (float|None), std (float|None), method (str), warnings (list)
        every property also has lo, hi (NaN if no valid margin), model_name, in_domain, tuple
        T_K (float), warnings (list of global warnings)
    """
    if composition_kind not in COMPOSITION_KINDS:
        raise ValueError(f"composition_kind must be one of {COMPOSITION_KINDS}, got {composition_kind!r}")
    if T_C is not None:
        T_K = float(T_C) + 273.15
    T_K = float(T_K)

    from fuelblend.api import predict_route
    if use_mixture_residual:
        _load_mix_models()                      # FileNotFoundError naming the build command when the bundle is missing
    route = predict_route.predict_fuel(list(composition), T_K=T_K, fuel_type=fuel_type,
                                       mole_fractions=mole_fractions, composition_kind=composition_kind)
    global_warn = list(route["warnings"])
    visc, dens, flash = route["viscosity_dyn_mPas"], route["density_kg_m3"], route["flash_point_C"]
    for res in (visc, dens, flash):
        res["warnings"] = list(_check_T_domain(T_K)) + res["warnings"]

    if use_mixture_residual:                    # opt-in GBM correction on top of the route value; no yaml margin exists for it
        inchikeys = [c["ik"] for c in route["components"]]
        w_mass = [c["w"] for c in route["components"]]
        x_mole = _mole_from_mass(inchikeys, w_mass) or w_mass
        mix_row = _make_mix_row(inchikeys, w_mass, x_mole, T_K, len(inchikeys))
        for prop, res in (("viscosity_dyn_mPas", visc), ("density_kg_m3", dens)):
            if res["value"] is None:
                continue
            mean, std = _mix_residual_predict(prop, mix_row)
            if mean is None:
                continue
            if prop == "viscosity_dyn_mPas":
                res["value"] = float(math.exp(math.log(res["value"]) + mean))
                res["std"] = float(res["value"] * abs(std)) if std else None
            else:
                res["value"] = float(res["value"] + mean)
                res["std"] = float(abs(std)) if std else None
            res["method"] = f"residual_gbm+{res['method']}"
            _attach_uncertainty(prop, res, inchikeys, T_K, fuel_type, composition_kind)

    unc = {p: {"value": None, "lo": math.nan, "hi": math.nan, "model_name": "none", "in_domain": False,
               "reason": r, "tuple": (None, math.nan, math.nan, "none", False)}
           for p, r in UNCERTIFIED_PROPERTIES.items()}

    return {
        "viscosity_dyn_mPas": visc,
        "density_kg_m3": dens,
        "flash_point_C": flash,
        "uncertified": unc,
        "T_K": T_K,
        "warnings": global_warn,
    }


def predict_tuples(*args, **kwargs) -> dict[str, tuple]:
    """Same inputs as ``predict``; returns {property: (value, lo, hi, model_name, in_domain)}.

    lo/hi are NaN when no valid yaml margin exists (never invented); in_domain is False for non-jet / undeclared fuel
    types, out-of-window temperature, unresolved or non-validated compounds, and uncertified models/properties.
    """
    r = predict(*args, **kwargs)
    out = {p: r[p]["tuple"] for p in ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C")}
    out.update({p: v["tuple"] for p, v in r["uncertified"].items()})
    return out


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
          f"interval=({fp['lo']:.1f}, {fp['hi']:.1f}), in_domain={fp['in_domain']})")
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
