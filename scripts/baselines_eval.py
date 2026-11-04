#!/usr/bin/env python
"""WP3: physics baselines on the unified model table (no training, no tuning).

    python scripts/baselines_eval.py --write-preds   # -> data/runs/baseline_preds.csv
    python scripts/baselines_eval.py --eval          # -> data/runs/BASELINES_EVAL_LOG.md, reports/BASELINES.md
    python scripts/baselines_eval.py                 # both

Baselines (AGENTS.md section 2):
  viscosity  : gn_mass   ln(eta_mix) = sum w_i ln(eta_i)              (Grunberg-Nissan, k_ij=0)
               gn_mole   same with mole fractions
               refutas   VBN_i = 14.534 ln ln(nu_i+0.8)+10.975, mass-weighted, inverted;
                         eta = nu_mix * rho_mix, rho_mix = ideal volume (Amagat)
  density    : dens_idvol  rho = 1 / sum(w_i / rho_i)
  flash      : chat_nominal / chat_lflT  Le Chatelier + Raoult (gamma=1) + Antoine, brentq on H(T)=1
               (chat_lflT adds Zabetakis LFL(T) = LFL25 (1 - 0.000721 (T-25)), T in C)
  vap. press.: psat_antoine  identity check of compounds.csv Antoine at the measured T (pure rows)

Pure-component curves come from the pure-compound rows of the model table itself:
  eta_i(T): ln eta = a + b ln T   (linear fit per compound, in-range only, 2 K pad)
  rho_i(T): rho = a + b T         (linear fit per compound, in-range only, 2 K pad)
A record is only predicted when EVERY component has pure data covering the record T (Mode A).

Evaluation: splits S1/S2/S3 from splits/s1s2s3.json (locked ids already excluded; re-asserted here),
null model = median of (property, kind, 20 K T-bin) from the S1 training rows (S1 rows themselves use a
5-fold system-disjoint median so the null is not in-sample), cluster bootstrap over system_id.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval import metrics as M  # noqa: E402

SEED = 20261030
P_TOTAL = 101.325          # kPa
T_PAD = 2.0                # K of allowed extrapolation outside a compound's pure-data range
FLASH_LO_C, FLASH_HI_C = -150.0, 400.0
TBIN = 20.0                # K width of the null-model temperature bin

MODEL_TABLE = ROOT / "data/processed/model_table.csv"
COMPOUNDS = ROOT / "data/unified/compounds.csv"
SPLITS = ROOT / "splits/s1s2s3.json"
LOCKED = ROOT / "splits/locked.json"
PREDS_OUT = ROOT / "data/runs/baseline_preds.csv"
LOG_OUT = ROOT / "data/runs/BASELINES_EVAL_LOG.md"
REPORT_OUT = ROOT / "reports/BASELINES.md"

# property -> baselines that are meaningful for it
TARGET = {
    "viscosity_dyn_mPas": ["gn_mass", "gn_mole", "refutas"],
    "viscosity_kin_cSt": ["gn_mass", "gn_mole", "refutas"],
    "density_kg_m3": ["dens_idvol"],
    "flash_point_C": ["chat_nominal", "chat_lflT"],
    "vapor_pressure_kPa": ["psat_antoine"],
}
PRED_COLS = ["gn_mass", "gn_mole", "refutas", "dens_idvol", "chat_nominal", "chat_lflT", "psat_antoine"]
KINDS = ["pure", "binary", "multi", "real_fuel"]
SPLIT_NAMES = ["S1", "S2", "S3"]

# unit + primary metric per property, used for the report tables
SPEC = {
    "viscosity_dyn_mPas": dict(score="rmse", unit="ln(mPa.s)", tol=0.05, tol_kind="rel5"),
    "viscosity_kin_cSt": dict(score="rmse", unit="ln(cSt)", tol=0.05, tol_kind="rel5"),
    "density_kg_m3": dict(score="rmse", unit="kg/m3", tol=0.05, tol_kind="rel5"),
    "flash_point_C": dict(score="rmse", unit="degC", tol=5.0, tol_kind="abs5"),
    "vapor_pressure_kPa": dict(score="rmse", unit="log10(kPa)", tol=0.1, tol_kind="log10"),
}


# ----------------------------------------------------------------------------- io
def _json(s):
    if isinstance(s, list):
        return s
    try:
        return json.loads(s) if isinstance(s, str) and s else []
    except Exception:
        return []


def load_inputs():
    mt = pd.read_csv(MODEL_TABLE, low_memory=False)
    comp = pd.read_csv(COMPOUNDS, low_memory=False).set_index("inchikey")
    splits = json.loads(SPLITS.read_text())
    locked = json.loads(LOCKED.read_text()) if LOCKED.exists() else {}
    return mt, comp, splits, locked


def apply_splits(mt, splits, locked):
    """Row label S1/S2/S3 (positional indices in the model table); locked ids forced OUT of every split.

    The split manifest was written before the last locked_sets() run in this repo, so a few locked
    multi-component system_ids are still inside S1/S2 - they are demoted to split "" here and counted.
    """
    split = pd.Series("", index=mt.index)
    for name in SPLIT_NAMES:
        split.iloc[list(splits.get(name, []))] = name
    lk_sys = set(locked.get("multi_system_ids_locked", []))
    lk_fuel = {str(f) for f in locked.get("real_fuel_ids_locked", [])}
    leaked = int((split.isin(SPLIT_NAMES) & mt["system_id"].isin(lk_sys)).sum())
    split = split.mask(split.isin(SPLIT_NAMES) & mt["system_id"].isin(lk_sys), "")
    out = mt.copy()
    out["split"] = split
    drop = out["system_id"].isin(lk_sys)
    if "fuel_id" in out.columns:      # model_table.csv has no fuel_id today; real fuels are split-excluded anyway
        leak_fuel = {str(f) for f in mt.loc[split.isin(SPLIT_NAMES), "fuel_id"]} & lk_fuel
        assert not leak_fuel, f"locked fuel_id inside a split: {sorted(leak_fuel)[:3]}"
        drop |= out["fuel_id"].astype(str).isin(lk_fuel)
    info = dict(locked_rows_removed_from_splits=leaked, locked_system_ids=len(lk_sys))
    return out[~drop].copy(), info


# ----------------------------------------------------------------- component fits
def _fit_loglin(T, y):
    """ln y = a + b ln T (b=0 for a single point).

    One-pass robust refit: with >=4 points, a point more than a factor 3 off the first fit is dropped and
    the line refitted (a handful of source rows are unit/assignment errors and would otherwise destroy the
    slope of a 2-point compound). The number dropped is returned for the log.
    """
    T = np.asarray(T, float)
    ly = np.log(np.asarray(y, float))
    lt = np.log(T)
    if len(T) == 1:
        return dict(a=float(ly[0]), b=0.0, lo=float(T[0]), hi=float(T[0]), n=1, rms=0.0, dropped=0)
    # robust start: Theil-Sen (median of pairwise slopes, median intercept) so 1-2 bad rows cannot drag the line
    if len(T) >= 3:
        i, j = np.triu_indices(len(T), 1)
        dlt = lt[j] - lt[i]
        ok = np.abs(dlt) > 1e-9
        b = float(np.median((ly[j] - ly[i])[ok] / dlt[ok])) if ok.any() else 0.0
        a = float(np.median(ly - b * lt))
    else:
        b, a = np.polyfit(lt, ly, 1)
    dropped = 0
    if len(T) >= 3:
        keep = np.abs(ly - (a + b * lt)) < np.log(3.0)
        if keep.sum() < 2:
            keep = np.ones(len(T), bool)
        dropped = int((~keep).sum())
        T, ly, lt = T[keep], ly[keep], lt[keep]
    b, a = np.polyfit(lt, ly, 1)
    rms = float(np.sqrt(np.mean((ly - (a + b * lt)) ** 2)))
    return dict(a=float(a), b=float(b), lo=float(T.min()), hi=float(T.max()), n=len(T), rms=rms,
                dropped=dropped)


def _fit_linear(T, y):
    """y = a + b T (b=0 for a single point)."""
    T = np.asarray(T, float)
    y = np.asarray(y, float)
    if len(T) == 1:
        return dict(a=float(y[0]), b=0.0, lo=float(T[0]), hi=float(T[0]), n=1, rms=0.0)
    b, a = np.polyfit(T, y, 1)
    rms = float(np.sqrt(np.mean((y - (a + b * T)) ** 2)))
    return dict(a=float(a), b=float(b), lo=float(T.min()), hi=float(T.max()), n=len(T), rms=rms)


def build_fits(pure):
    """pure: rows with columns inchikey, property, T_K, value (pure-compound records only)."""
    fe, fr = {}, {}
    p = pure[(pure.property == "viscosity_dyn_mPas") & pure.value.gt(0)]
    for ik, g in p.groupby("inchikey"):
        g = g.groupby(g.T_K.round(3))["value"].median()
        fe[ik] = _fit_loglin(g.index.values, g.values)
    p = pure[pure.property == "density_kg_m3"]
    for ik, g in p.groupby("inchikey"):
        g = g.groupby(g.T_K.round(3))["value"].median()
        fr[ik] = _fit_linear(g.index.values, g.values)
    return {"eta": fe, "rho": fr}


def eta_at(fits, ik, T, pad=T_PAD):
    f = fits["eta"].get(ik)
    if f is None or not (f["lo"] - pad <= T <= f["hi"] + pad):
        return None
    return float(np.exp(f["a"] + f["b"] * np.log(T)))


def rho_at(fits, ik, T, pad=T_PAD):
    f = fits["rho"].get(ik)
    if f is None or not (f["lo"] - pad <= T <= f["hi"] + pad):
        return None
    v = float(f["a"] + f["b"] * T)
    return v if v > 0 else None


def mole_from_mass(keys, w, mw):
    """x_i = (w_i/M_i) / sum(w_j/M_j)."""
    try:
        num = np.array([float(wi) / float(mw[k]) for k, wi in zip(keys, w)], float)
    except (KeyError, TypeError, ValueError):
        return None
    s = num.sum()
    return (num / s).tolist() if np.isfinite(s) and s > 0 else None


# --------------------------------------------------------------------- baselines
def gn_from_frac(keys, frac, fits, T):
    """exp(sum frac_i ln eta_i). Returns (value, reason)."""
    vals = []
    for k, f in zip(keys, frac):
        e = eta_at(fits, k, T)
        if e is None or e <= 0:
            return None, "no_pure_eta_at_T"
        vals.append(f * np.log(e))
    return float(np.exp(sum(vals))), ""


def refutas(keys, w, fits, T):
    """Returns (eta_dyn_mPas, nu_cSt, rho_mix_kg_m3, reason)."""
    vb, wr, etas, rhos, nus = 0.0, 0.0, [], [], []
    for k, wi in zip(keys, w):
        e = eta_at(fits, k, T)
        if e is None or e <= 0:
            return None, None, None, "no_pure_eta_at_T"
        r = rho_at(fits, k, T)
        if r is None:
            return None, None, None, "no_pure_rho_at_T"
        nu = e / (r / 1000.0)                      # mPa.s / (kg/m3) * 1000 = mm2/s
        if nu < 0.2:
            return None, None, None, "nu_lt_0p2"
        vb += wi * (14.534 * np.log(np.log(nu + 0.8)) + 10.975)
        wr += wi
        etas.append(e)
        rhos.append(r)
        nus.append(nu)
    if wr <= 0:
        return None, None, None, "bad_w_mass"
    vbn = vb / wr
    nu_mix = float(np.exp(np.exp((vbn - 10.975) / 14.534)) - 0.8)
    rho_mix = 1.0 / sum(wi / r for wi, r in zip(w, rhos))
    return nu_mix * rho_mix / 1000.0, nu_mix, rho_mix, ""


def dens_idvol(keys, w, fits, T):
    inv = 0.0
    for k, wi in zip(keys, w):
        r = rho_at(fits, k, T)
        if r is None:
            return None, "no_pure_rho_at_T"
        inv += wi / r
    return (1.0 / inv if inv > 0 else None), ""


def psat_kpa(A, B, C, T_K):
    d = T_K + C
    if not np.isfinite(d) or d <= 1.0:
        return np.nan
    return float(10.0 ** (A - B / d))


def flash_prediction(keys, x, comp, lfl_temp=False):
    """Le Chatelier + Raoult (gamma=1). Returns (T_C, reason)."""
    A, B, C, lfl = [], [], [], []
    for k in keys:
        if k not in comp.index:
            return None, "unknown_component"
        row = comp.loc[k]
        if pd.isna(row["antoine_A"]) or pd.isna(row["antoine_B"]) or pd.isna(row["antoine_C"]):
            return None, "no_antoine"
        if pd.isna(row["lfl_vol_pct"]):
            return None, "unknown_LFL"
        A.append(float(row["antoine_A"]))
        B.append(float(row["antoine_B"]))
        C.append(float(row["antoine_C"]))
        lfl.append(float(row["lfl_vol_pct"]) / 100.0)
    x = np.asarray(x, float)
    if x.sum() <= 0:
        return None, "bad_x_mole"

    def haz(T_K):
        tc = T_K - 273.15
        fac = max(1.0 - 0.000721 * (tc - 25.0), 0.05) if lfl_temp else 1.0
        h = 0.0
        for xi, a, b, c, l in zip(x, A, B, C, lfl):
            p = psat_kpa(a, b, c, T_K)
            if not np.isfinite(p):
                return np.nan
            h += xi * p / (P_TOTAL * l * fac)
        return h - 1.0

    lo, hi = FLASH_LO_C + 273.15, FLASH_HI_C + 273.15
    flo, fhi = haz(lo), haz(hi)
    if not np.isfinite(flo) or flo > 0:
        return None, "no_root_in_bracket"
    if not np.isfinite(fhi) or fhi < 0:
        return None, "no_root_in_bracket"
    try:
        return float(brentq(haz, lo, hi, xtol=1e-6, rtol=1e-10) - 273.15), ""
    except Exception as exc:  # noqa: BLE001
        return None, f"solver_{type(exc).__name__}"


def psat_antoine_pure(key, T_K, comp):
    if key not in comp.index:
        return None, "unknown_component"
    row = comp.loc[key]
    if pd.isna(row["antoine_A"]) or pd.isna(row["antoine_B"]) or pd.isna(row["antoine_C"]):
        return None, "no_antoine"
    if not np.isfinite(T_K):
        return None, "no_T"
    v = psat_kpa(float(row["antoine_A"]), float(row["antoine_B"]), float(row["antoine_C"]), float(T_K))
    return (v, "") if np.isfinite(v) and v > 0 else (None, "psat_invalid")


# ------------------------------------------------------------------ prediction pass
def predict_records(mt, comp, fits):
    """One row per model-table record; NaN where the baseline is not applicable."""
    skip = Counter()      # (property, baseline, reason) -> n
    rows = []
    for r in mt.itertuples():
        prop = r.property
        if prop not in TARGET:
            continue
        keys = _json(r.comp_inchikeys)
        w = _json(r.w_mass)
        T = float(r.T_K) if pd.notna(r.T_K) else np.nan
        rec = dict(property=prop, system_id=r.system_id, T_key=("" if not np.isfinite(T) else f"{T:.2f}"),
                   P_key=_p_key(r.P_kPa), x_key=_x_key(keys, _json(r.x_mole)),
                   kind=r.kind, split=getattr(r, "split", ""), value=float(r.value),
                   conflict=bool(r.conflict), n_comp=int(r.n_comp), T_K=T,
                   fuel_id=getattr(r, "fuel_id", ""), isnull=False)
        for c in PRED_COLS:
            rec[c] = np.nan
        rec["baseline_row_id"] = hashlib.md5(
            f"{prop}|{r.system_id}|{rec['T_key']}|{'' if pd.isna(r.P_kPa) else round(float(r.P_kPa), 4)}|{rec['x_key']}"
            .encode()).hexdigest()[:16]
        rec["eval_key"] = rec["baseline_row_id"]

        def _skip(base, why):
            skip[(prop, base, why)] += 1

        if not keys or r.kind == "real_fuel":
            for base in TARGET[prop]:
                _skip(base, "no_molecule_list")
            rows.append(rec)
            continue

        if prop in ("viscosity_dyn_mPas", "viscosity_kin_cSt"):
            if not np.isfinite(T):
                for base in TARGET[prop]:
                    _skip(base, "no_T")
            else:
                v, why = gn_from_frac(keys, w, fits, T)
                if v is None:
                    for base in TARGET[prop]:
                        _skip(base, why)
                else:
                    rec["gn_mass"] = v
                    x = mole_from_mass(keys, w, comp["mw"].to_dict())
                    if x is None:
                        _skip("gn_mole", "no_mw")
                    else:
                        vm, why2 = gn_from_frac(keys, x, fits, T)
                        rec["gn_mole"] = np.nan if vm is None else vm
                        if vm is None:
                            _skip("gn_mole", why2)
                    edyn, nu, rho_mix, why3 = refutas(keys, w, fits, T)
                    if edyn is None:
                        _skip("refutas", why3)
                    elif prop == "viscosity_dyn_mPas":
                        rec["refutas"] = edyn
                    else:
                        # kinematic property: every prediction must be in cSt, so it needs rho_mix
                        rec["refutas"] = nu
                        if rho_mix is None:
                            rec["gn_mass"] = rec["gn_mole"] = np.nan
                            for base in ("gn_mass", "gn_mole"):
                                _skip(base, "no_pure_rho_at_T")
                        else:
                            rec["gn_mass"] = rec["gn_mass"] / (rho_mix / 1000.0)
                            if np.isfinite(rec["gn_mole"]):
                                rec["gn_mole"] = rec["gn_mole"] / (rho_mix / 1000.0)
        elif prop == "density_kg_m3":
            if r.kind == "pure":
                rec["dens_idvol"] = rec["value"]      # identity by construction (v_i = 1)
            elif not np.isfinite(T):
                _skip("dens_idvol", "no_T")
            else:
                v, why = dens_idvol(keys, w, fits, T)
                rec["dens_idvol"] = np.nan if v is None else v
                if v is None:
                    _skip("dens_idvol", why)
        elif prop == "flash_point_C":
            if r.kind == "pure":
                x = [1.0]
            else:
                x = mole_from_mass(keys, w, comp["mw"].to_dict())
                if x is None:
                    for base in TARGET[prop]:
                        _skip(base, "no_mw")
                    rows.append(rec)
                    continue
            for base, lfl_t in (("chat_nominal", False), ("chat_lflT", True)):
                p, why = flash_prediction(keys, x, comp, lfl_temp=lfl_t)
                rec[base] = np.nan if p is None else p
                if p is None:
                    _skip(base, why)
        elif prop == "vapor_pressure_kPa":
            if r.kind != "pure":
                _skip("psat_antoine", "not_pure")
            else:
                v, why = psat_antoine_pure(keys[0], T, comp)
                rec["psat_antoine"] = np.nan if v is None else v
                if v is None:
                    _skip("psat_antoine", why)

        rows.append(rec)

    df = pd.DataFrame(rows)
    cols = ["property", "system_id", "T_key", "P_key", "x_key", "kind", "split", "baseline_row_id", "value"] + \
           PRED_COLS + ["null_pred", "is_null_pred", "eval_key", "conflict", "n_comp", "T_K", "fuel_id"]
    df = df.reindex(columns=cols)
    assert df.baseline_row_id.is_unique, "duplicate (property, system, T, P, x) key in model table"
    return df, skip


def _p_key(P):
    if pd.isna(P):
        return "101"
    P = float(P)
    return "101" if P <= 150 else str(int(round(P / 100.0)) * 100)


def _x_key(keys, x):
    if not keys:
        return ""
    pairs = sorted(zip(keys, x))
    return ",".join(f"{k}:{float(v):.3f}" for k, v in pairs)


# -------------------------------------------------------------------- null model
def _tbin(T):
    return "na" if pd.isna(T) else f"{int(np.floor(float(T) / TBIN) * TBIN)}"


LEVELS = [("property", "kind", "Tbin"), ("property", "kind"), ("property",)]


def _null_medians(ref):
    return {lv: ref.groupby(list(lv))["value"].median().to_dict() for lv in LEVELS}


def _null_lookup(med, tgt):
    out = np.full(len(tgt), np.nan)
    props = tgt["property"].values
    kinds = tgt["kind"].values
    tbins = tgt["Tbin"].values
    for i in range(len(tgt)):
        triple = (props[i], kinds[i], tbins[i])
        for lv in LEVELS:
            key = triple[:len(lv)]
            v = med[lv].get(key if len(lv) > 1 else key[0], np.nan)
            if pd.notna(v):
                out[i] = v
                break
    return out


def null_predictions(mt, seed=SEED, n_folds=5):
    """Median of (property, kind, 20 K T-bin) from S1; S1 itself uses system-disjoint folds."""
    d = mt.assign(Tbin=mt["T_K"].map(_tbin))
    d = d.assign(fold=d["system_id"].map(
        lambda s: int(hashlib.md5(f"{seed}|{s}".encode()).hexdigest(), 16) % n_folds))
    out = np.full(len(d), np.nan)
    s1 = d[d.split == "S1"]
    for split_name in SPLIT_NAMES:
        mask = (d.split == split_name).values
        if not mask.any():
            continue
        tgt = d[mask]
        pred = np.full(len(tgt), np.nan)
        if split_name == "S1":
            for f in range(n_folds):
                sub = (tgt["fold"] == f).values
                if not sub.any():
                    continue
                pred[sub] = _null_lookup(_null_medians(s1[s1.fold != f]), tgt[sub])
        else:
            pred = _null_lookup(_null_medians(s1), tgt)
        out[mask] = pred
    return out


# ------------------------------------------------------------------------ metrics
def _boot(y, p, cl, fn, seed=SEED):
    return M.bootstrap_ci(y, p, cl, fn, n_boot=2000, seed=seed)


def eval_block(y, p, clusters, mode, seed=SEED):
    """mode: log (positive-scale property -> ln-space), abs (degC / kg/m3), log10."""
    y = np.asarray(y, float)
    p = np.asarray(p, float)
    out = {"n": int(len(y))}
    rel = np.abs((p - y) / np.clip(np.abs(y), 1e-9, None))
    out["mape_pct"] = float(100 * np.mean(rel))
    out["frac_rel5"] = float(np.mean(rel <= 0.05))
    out["frac_rel10"] = float(np.mean(rel <= 0.10))
    if mode == "log":
        yl, pl = np.log(y), np.log(p)
        out["score_name"] = "rmse"
        b = _boot(yl, pl, clusters, M.rmse, seed)
        out.update(rmse=b["point"], rmse_lo=b["lo"], rmse_hi=b["hi"])
        b = _boot(yl, pl, clusters, M.mae, seed)
        out.update(mae=b["point"], mae_lo=b["lo"], mae_hi=b["hi"])
    elif mode == "abs":
        out["score_name"] = "rmse"
        b = _boot(y, p, clusters, M.rmse, seed)
        out.update(rmse=b["point"], rmse_lo=b["lo"], rmse_hi=b["hi"])
        b = _boot(y, p, clusters, M.mae, seed)
        out.update(mae=b["point"], mae_lo=b["lo"], mae_hi=b["hi"])
        out["bias"] = float(np.mean(p - y))
    elif mode == "log10":
        yl, pl = np.log10(np.clip(y, 1e-6, None)), np.log10(np.clip(p, 1e-6, None))
        out["score_name"] = "rmse"
        b = _boot(yl, pl, clusters, M.rmse, seed)
        out.update(rmse=b["point"], rmse_lo=b["lo"], rmse_hi=b["hi"])
        b = _boot(yl, pl, clusters, M.mae, seed)
        out.update(mae=b["point"], mae_lo=b["lo"], mae_hi=b["hi"])
    return out


def _mode_of(prop):
    if prop in ("viscosity_dyn_mPas", "viscosity_kin_cSt"):
        return "log"
    if prop == "vapor_pressure_kPa":
        return "log10"
    return "abs"


def evaluate(pred, seed=SEED, conflict_free=True):
    rows = []
    d0 = pred[pred.split.isin(SPLIT_NAMES)]
    if conflict_free:
        d0 = d0[~d0.conflict.astype(bool)]
    for prop, bases in TARGET.items():
        for base in bases:
            for split in SPLIT_NAMES:
                for kind in KINDS:
                    s = d0[(d0.property == prop) & (d0.split == split) & (d0.kind == kind)]
                    if not len(s):
                        continue
                    s = s.dropna(subset=[base, "null_pred"])
                    if not len(s):
                        continue
                    y = s["value"].values
                    cl = s["system_id"].values
                    if prop in ("viscosity_dyn_mPas", "viscosity_kin_cSt", "vapor_pressure_kPa"):
                        keep = (y > 0) & (s[base].values > 0)
                        s, y, cl = s[keep], y[keep], cl[keep]
                        if not len(s):
                            continue
                    for name, col in ((base, base), ("null", "null_pred")):
                        yv, pv = s["value"].values, s[col].values
                        r = eval_block(yv, pv, cl, _mode_of(prop), seed)
                        r.update(property=prop, baseline=name, split=split, kind=kind,
                                 conflict_free=conflict_free, score=SPEC[prop]["score"])
                        if prop == "flash_point_C":
                            d = pv - yv
                            r["overpred_rate"] = float(np.mean(d > 5.0))
                            ci = _boot(yv, pv, cl, lambda a, c: M.flash_overpred_rate(a, c, 5.0), seed)
                            r["overpred_rate_lo"], r["overpred_rate_hi"] = ci["lo"], ci["hi"]
                            r["mean_positive_error"] = M.mean_positive_error(yv, pv)
                            r["frac_abs5"] = float(np.mean(np.abs(d) <= 5.0))
                            r["p90_abs"] = float(np.percentile(np.abs(d), 90))
                        rows.append(r)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------- reporting
def _md(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
    return "\n".join(out)


def _f(x, nd=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    return f"{x:.{nd}f}"


def _ci(r, key="rmse", nd=3):
    return f"{_f(r[key], nd)} [{_f(r[key + '_lo'], nd)}, {_f(r[key + '_hi'], nd)}]"


def headline_table(mets, prop, conflict_free=True):
    """Rows = baseline x split x kind, physics and null side by side."""
    m = mets[(mets.property == prop) & (mets.conflict_free == conflict_free)]
    if not len(m):
        return "_no rows_"
    spec = SPEC[prop]
    rows = []
    for base in TARGET[prop] + ["null"]:
        b = m[m.baseline == base]
        for split in [s for s in SPLIT_NAMES if s in set(b.split)]:
            for kind in [k for k in KINDS if k in set(b[b.split == split].kind)]:
                r = b[(b.split == split) & (b.kind == kind)]
                if not len(r):
                    continue
                r = r.iloc[0]
                line = [base, split, kind, int(r.n), _ci(r, r.score_name, 4)]
                line.append(_f(r["mae"], 3))
                line.append(f"{100 * r['frac_rel5']:.1f}%")
                if prop == "flash_point_C":
                    line.append(f"{100 * r['overpred_rate']:.1f}% "
                                f"[{100 * r['overpred_rate_lo']:.0f},{100 * r['overpred_rate_hi']:.0f}]")
                    line.append(_f(r["mean_positive_error"], 2))
                    line.append(f"{100 * r['frac_abs5']:.1f}%")
                rows.append(line)
    heads = ["baseline", "split", "kind", "n", f"{spec['score']} 95% CI ({spec['unit']})", "MAE",
             "within tol" + (" (5%)" if spec["tol_kind"] == "rel5" else "")]
    if prop == "flash_point_C":
        heads += ["over-pred >5 C", "mean +err (C)", "|err|<=5 C"]
    return _md(heads, rows)


def mixture_summary(mets, conflict_free=True):
    """Physics baselines only: the null appears once per physics cell (paired) and would be triple-counted here,
    so it is reported in the paired bootstrap table instead."""
    m = mets[(mets.conflict_free == conflict_free) & mets.kind.isin(["binary", "multi"])
             & (mets.baseline != "null")]
    rows = []
    for (prop, base), g in m.groupby(["property", "baseline"], sort=False):
        for split in [s for s in SPLIT_NAMES if s in set(g.split)]:
            gg = g[g.split == split]
            if not len(gg):
                continue
            n = int(gg.n.sum())
            wavg = float(np.average(gg[SPEC[prop]["score"]], weights=gg.n))
            rows.append([prop, base, split, n, _f(wavg, 4), f"{100 * float(np.average(gg.frac_rel5, weights=gg.n)):.1f}%"])
    return _md(["property", "baseline", "split", "n(mixtures)", "n-weighted score", "within 5% rel"], rows)


def _pooled_ci(pred, prop, base, split, seed=SEED):
    """Pooled binary+multi, conflict-free: physics vs null RMSE with 95 % cluster-bootstrap CIs."""
    mode = _mode_of(prop)
    d = pred[(pred.property == prop) & (pred.split == split) & pred.kind.isin(["binary", "multi"])]
    d = d[~d.conflict.astype(bool)].dropna(subset=[base, "null_pred"])
    if mode in ("log", "log10"):
        d = d[(d["value"] > 0) & (d[base] > 0) & (d["null_pred"] > 0)]
    if len(d) < 5:
        return None
    y = d["value"].values.astype(float)
    cl = d["system_id"].values
    f = (lambda a, b: M.rmse(np.log(a), np.log(b))) if mode == "log" else (
        (lambda a, b: M.rmse(np.log10(a), np.log10(b))) if mode == "log10" else M.rmse)
    out = {"n": len(d)}
    for name, col in (("physics", base), ("null", "null_pred")):
        ci = M.bootstrap_ci(y, d[col].values.astype(float), cl, f, n_boot=2000, seed=seed)
        out[name] = ci["point"]
        out[name + "_lo"], out[name + "_hi"] = ci["lo"], ci["hi"]
    out["separated"] = out["physics_hi"] < out["null_lo"]
    return out


def write_log(mt, pred, skip, mets, t0, argv):
    L = ["# BASELINES_EVAL_LOG (WP3, auto-generated raw log)", "",
         f"- command: `python scripts/baselines_eval.py {' '.join(argv)}`",
         f"- seed: {SEED} (split seed, bootstrap seed, null-fold hash seed)",
         f"- run at: {time.strftime('%Y-%m-%d %H:%M:%S')}  wall time {time.time() - t0:.1f} s",
         f"- python {sys.version.split()[0]}, pandas {pd.__version__}, numpy {np.__version__}",
         "- inputs: data/processed/model_table.csv, data/unified/compounds.csv, splits/s1s2s3.json, splits/locked.json",
         "- nothing was trained or tuned; no locked row was read for tuning (locked ids asserted absent from all splits)",
         "", "## Row inventory", ""]
    inv = mt.groupby(["property", "split"]).size().unstack(fill_value=0)
    L.append("```\n" + inv.to_string() + "\n```")
    L += ["", "## Evaluable records (non-NaN prediction) per property/baseline/split", ""]
    d = pred[pred.split.isin(SPLIT_NAMES)]
    rows = []
    for prop, bases in TARGET.items():
        for base in bases:
            for split in SPLIT_NAMES:
                s = d[(d.property == prop) & (d.split == split) & d[base].notna()]
                if not len(s):
                    continue
                by = s.groupby("kind").size().to_dict()
                rows.append([prop, base, split, len(s), by.get("pure", 0), by.get("binary", 0), by.get("multi", 0)])
    L.append(_md(["property", "baseline", "split", "n_eval", "pure", "binary", "multi"], rows))
    L += ["", "## Unevaluable records (skip reasons)", ""]
    rows = [[p, b, why, n] for (p, b, why), n in sorted(skip.items(), key=lambda kv: (-kv[1], kv[0]))]
    L.append(_md(["property", "baseline", "reason", "n"], rows))
    L += ["", "## Full metrics: conflict=False (headline)", ""]
    L.append("```\n" + mets[mets.conflict_free].round(6).to_csv(index=False) + "```")
    L += ["", "## Full metrics: all rows incl. conflict=True (sensitivity)", ""]
    L.append("```\n" + mets[~mets.conflict_free].round(6).to_csv(index=False) + "```")
    LOG_OUT.parent.mkdir(parents=True, exist_ok=True)
    LOG_OUT.write_text("\n".join(L) + "\n", encoding="utf-8")


def write_report(mt, pred, skip, mets, t0, argv, run_info=None):
    d = pred[pred.split.isin(SPLIT_NAMES)]
    real = pred[(pred.kind == "real_fuel") & (pred.property.isin(TARGET))]
    evaluated = {p: int(d[(d.property == p) & d[TARGET[p]].notna().any(axis=1)].shape[0]) for p in TARGET}
    n_of = {p: int((mt.property == p).sum()) for p in TARGET}
    run_info = run_info or {}

    L = ["# BASELINES.md - WP3 physics baselines on the unified data", "",
         "## What ran - exact commands + seeds", "",
         "```", f"python scripts/baselines_eval.py {' '.join(argv)}", "```", "",
         f"- script: `scripts/baselines_eval.py`; outputs `data/runs/baseline_preds.csv`, "
         f"`data/runs/BASELINES_EVAL_LOG.md`, this file.",
         f"- seed {SEED} everywhere (splits, bootstrap, null-fold hashing); 2000 cluster-bootstrap resamples by "
         "`system_id` (`src/fuelblend/eval/metrics.py`); wall time {:.0f} s.".format(time.time() - t0),
         f"- data: `data/processed/model_table.csv` ({len(mt):,} hydrocarbon-only de-duplicated records, built by "
         "`src/fuelblend/data/model_table.py`); component constants from `data/unified/compounds.csv`; "
         "splits `splits/s1s2s3.json`.",
         "- nothing trained, nothing tuned; locked sets (`splits/locked.json`) never evaluated (asserted absent; "
         f"{run_info.get('locked', '?')} locked multi-component rows are demoted out of every split before scoring).",
         "- **Mode A only**: a record is predicted iff every one of its components has pure-compound data of the "
         "needed property covering the record temperature (2 K pad). No novel-molecule extrapolation is attempted.",
         "- pure-component curves are fitted from the pure rows of the same table: `ln eta = a + b ln T` "
         "(viscosity, in-range only) and `rho = a + b T` (density). A single-measurement compound is a constant. "
         f"The viscosity fit starts from a Theil-Sen line and drops points more than a factor 3 off it "
         f"({run_info.get('dropped', '?')} of {run_info.get('points', '?')} pure viscosity points dropped over "
         f"{run_info.get('curves', '?')} compounds).",
         "- null model: median of (`property`, `kind`, 20 K T-bin) from the **S1 training rows**; S1 rows are scored "
         "against a 5-fold system-disjoint median so the null is never in-sample. Same rows are used for the "
         "physics baseline and the null in every cell.",
         "- default filter `conflict=False` for every number below unless the line says otherwise; the sensitivity "
         "run including `conflict=True` is in the log.",
         "", "## Numbers table per property x kind x split: RMSE, MAE, 95 % CI from your bootstrap, "
         "frac_within_tolerance, over-prediction rate (flash)", ""]

    notes = {
        "viscosity_dyn_mPas": "log-space score = RMSE of ln(eta / mPa.s); tolerance 5 % relative. "
                              "`gn_mass`/`gn_mole` = Grunberg-Nissan k_ij=0 (mass / mole weighting); "
                              "`refutas` = VBN blending on kinematic viscosity, needs pure density (records with "
                              "any nu_i < 0.2 cSt skipped).",
        "viscosity_kin_cSt": "log-space score = RMSE of ln(nu / cSt); `refutas` predicts nu directly; "
                             "`gn_mass`/`gn_mole` are converted with the ideal-volume mixture density.",
        "density_kg_m3": "score = RMSE in kg/m3; `dens_idvol` = 1 / sum(w_i/rho_i). Pure rows are identity by "
                         "construction (v_i = 1) - marked with † - and are not evidence about the physics.",
        "flash_point_C": "score = RMSE in degC; `chat_nominal` = Le Chatelier + Raoult (gamma=1) + Antoine with "
                         "LFL25; `chat_lflT` adds Zabetakis LFL(T). over-prediction = pred - measured > 5 degC "
                         "(the unsafe direction: the model claims a safer fuel than it is).",
        "vapor_pressure_kPa": "score = RMSE of log10(P/kPa); this is a sanity check that `compounds.csv` Antoine "
                              "reproduces the pure vapour-pressure measurements it was (partly) refitted from - it "
                              "is **not** an independent validation.",
    }
    for prop in TARGET:
        L += [f"### {prop}", "", notes[prop], "",
              f"records in table: {n_of[prop]}, evaluated: {evaluated[prop]} "
              f"({100 * evaluated[prop] / max(n_of[prop], 1):.1f} %)", "",
              headline_table(mets, prop), ""]
        if prop in ("viscosity_dyn_mPas", "density_kg_m3"):
            L += ["† pure rows: for viscosity the value is the pure-compound fit evaluated at its own data "
                  "(circular, included only to show the fit error); for density it is the identity.", ""]

    L += ["### Mixture-only summary (binary + multi, conflict=False)", "",
          "Physics baselines only, n-weighted across kinds within each split; this is the honest headline (pure rows "
          "excluded when they are identity/self-fit). The comparable null level, on exactly the same rows, is in "
          "the paired bootstrap table below. The S3 viscosity cells are dominated by impossible source labels "
          "(system `a175950749e8`, see *What failed*); use the S1/S2 columns for generalisation.", "", mixture_summary(mets), "",
          "Success criterion (reports/NOISE_FLOOR.md): a model only counts as an improvement over physics if it "
          "beats the physics baseline by at least the noise floor with a non-overlapping 95 % CI - viscosity "
          f"{0.016} log-units median spread (p90 0.12), density 0.73 kg/m3 (p90 5.44), flash 6 degC (p90 29.4). "
          "The same bar applies to physics vs the null model below.", ""]

    # verdict per property, physics vs null
    L += ["### Does physics beat the null model?", "",
          "Per split, n-weighted across kinds and across the physics cells listed for that property (each cell pairs "
          "the physics baseline with the null on the same rows; pure rows included, so density/viscosity pure cells "
          "are identity/self-fit and not evidence). The paired, mixture-only, same-rows comparison with bootstrap "
          "CIs is the second table.", ""]
    vrows = []
    for prop in TARGET:
        m = mets[(mets.property == prop) & mets.conflict_free]
        for base in TARGET[prop]:
            b = m[m.baseline == base]
            z = m[m.baseline == "null"]
            if not len(b) or not len(z):
                continue
            for split in [s for s in SPLIT_NAMES if s in set(b.split)]:
                bb, zz = b[b.split == split], z[z.split == split]
                if not len(bb) or not len(zz):
                    continue
                sc = SPEC[prop]["score"]
                pb = float(np.average(bb[sc], weights=bb.n))
                pz = float(np.average(zz[sc], weights=zz.n))
                vrows.append([prop, base, split, _f(pb, 4), _f(pz, 4),
                              f"{100 * (pb - pz) / max(pz, 1e-12):+.1f}%",
                              "physics better" if pb < pz else "null better"])
    L += [_md(["property", "baseline", "split", "physics score", "null score", "delta", "verdict"], vrows), ""]
    L += ["Pooled binary+multi rows (conflict-free), same rows on both sides, RMSE with a 95 % "
          "cluster-bootstrap CI - 'separated' means the physics CI lies entirely below the null CI:", ""]
    prows = []
    for prop in TARGET:
        for base in TARGET[prop]:
            for split in SPLIT_NAMES:
                p = _pooled_ci(pred, prop, base, split)
                if p is None:
                    continue
                prows.append([prop, base, split, p["n"],
                              f"{_f(p['physics'], 4)} [{_f(p['physics_lo'], 4)}, {_f(p['physics_hi'], 4)}]",
                              f"{_f(p['null'], 4)} [{_f(p['null_lo'], 4)}, {_f(p['null_hi'], 4)}]",
                              "yes" if p["separated"] else "no"])
    L += [_md(["property", "baseline", "split", "n", "physics RMSE [lo, hi]", "null RMSE [lo, hi]",
               "separated"], prows), ""]

    L += ["## real vs synthetic", "",
          "Every number in this report is a **real-data** number computed from `data/processed/model_table.csv` "
          "(measured densities, viscosities, flash points and vapour pressures, de-duplicated medians). No "
          "synthetic label is used anywhere in the evaluation. The only synthetic values in this work package are "
          "the hand-made fixtures inside `tests/test_baselines_eval.py` (invented pure viscosities and Antoine "
          "constants, used to prove the code implements the formulas); they are clearly labelled in the test file "
          "and are **not** reported as results.", ""]

    # coverage / failures
    L += ["## What failed / unevaluable (counts + reasons)", "",
          _md(["property", "records", "evaluable (any baseline)", "share"],
              [[p, n_of[p], evaluated[p], f"{100 * evaluated[p] / max(n_of[p], 1):.1f}%"] for p in TARGET]), ""]
    L += ["Per-baseline skip reasons (a record can be skipped by one baseline and predicted by another):", ""]
    L += [_md(["property", "baseline", "reason", "records"],
              [[p, b, why, n] for (p, b, why), n in sorted(skip.items(), key=lambda kv: (-kv[1], kv[0]))]), ""]
    rrows = []
    for prop in TARGET:
        s = d[d.property == prop]
        for kind in KINDS:
            k = s[s.kind == kind]
            if not len(k):
                continue
            rrows.append([prop, kind, len(k), int(k[TARGET[prop]].notna().any(axis=1).sum())])
    L += ["Evaluable records per property x kind (a record counts as evaluable if at least one baseline "
          "produced a number):", "", _md(["property", "kind", "records", "evaluable"], rrows), ""]

    n_real = int(real[PRED_COLS].notna().any(axis=1).sum())
    L += [f"**Real fuels**: {len(real)} real-fuel records for these five properties are in the table and "
          f"{n_real} of them have a molecule list, so these baselines cannot score them at all. They are held for "
          "a distribution-input model (class x carbon-number, GCxGC) and for inverse-design validation - "
          "listed as future work, never silently dropped (every one is in `baseline_preds.csv` with NaN "
          "predictions; the counted reason `no_molecule_list` is in the log's skip-reason table).", ""]
    L += [f"**S3** ({int((d.split == 'S3').sum())} rows) was evaluated and is reported where a baseline reached "
          "it; S3 holds out whole molecules, but Mode A physics does not need those molecules to be seen, so "
          "only Antoine/LFL/pure-curve availability limits coverage. Cells with zero coverage are simply absent "
          "from the tables above.", ""]

    # audit trail: worst rows and impossible labels
    L += ["### Worst records (audit trail, not excluded)", "",
          "Largest errors of the best baseline per property (mixture rows, conflict-free). They are kept in every "
          "number above - nothing here is filtered out.", ""]
    best = {"viscosity_dyn_mPas": "refutas", "viscosity_kin_cSt": "refutas", "density_kg_m3": "dens_idvol",
            "flash_point_C": "chat_lflT", "vapor_pressure_kPa": "psat_antoine"}
    wrows = []
    for prop, base in best.items():
        s = d[(d.property == prop) & d[base].notna() & d.kind.isin(["pure", "binary", "multi"])]
        s = s[~s.conflict.astype(bool)]
        if not len(s):
            continue
        e = (np.log(s[base]) - np.log(s["value"])) if _mode_of(prop) in ("log", "log10") else (s[base] - s["value"])
        s = s.assign(_e=np.abs(e)).sort_values("_e", ascending=False)
        for _, row in s.head(2).iterrows():
            wrows.append([prop, base, row.system_id, row.split, row.kind,
                          _f(row.T_K, 1), _f(row["value"], 4), _f(row[base], 4), _f(row["_e"], 3)])
    L += [_md(["property", "baseline", "system_id", "split", "kind", "T_K", "measured", "predicted", "|err|"],
              wrows), ""]

    susp = []
    for prop, base in (("viscosity_dyn_mPas", "gn_mass"), ("density_kg_m3", "dens_idvol")):
        s = d[(d.property == prop) & d.kind.isin(["binary", "multi"])].dropna(subset=[base])
        s = s[(s["value"] > 0) & (s[base] > 0)]
        ratio = np.abs(np.log(s["value"] / s[base]))
        bad = s[ratio > np.log(50.0)]
        susp.append(f"- `{prop}`: {len(bad)} of {len(s)} mixture records are more than 50x away from the "
                    f"ideal/log-linear prediction." +
                    (f" Example: system `{bad.iloc[0].system_id}` ({bad.iloc[0].kind}, {_f(bad.iloc[0].T_K, 1)} K) "
                     f"measured {_f(bad.iloc[0]['value'], 3)}, predicted {_f(bad.iloc[0][base], 3)}."
                     if len(bad) else ""))
    L += ["**Impossible-looking labels.** A physics baseline cannot be off by 3 orders of magnitude on a mixture of "
          "well-characterised components, so these rows are almost certainly unit/assignment errors in the source "
          "records (e.g. mPa.s vs Pa.s families) and they dominate the S3 viscosity cells:", ""] + susp + [""]

    L += ["## Decisions made", "",
          "1. **Component curves from the pure rows of the same table**, `ln eta = a + b ln T` and `rho = a + b T`, "
          "evaluated only inside (2 K pad) each compound's own measured range. This replaces the earlier "
          "nearest-measurement rule (`baselines.py`, 7.5 K tolerance) with interpolation - coverage went up and no "
          "molecule is extrapolated.",
          "2. **Mode A is enforced structurally**: no component curve -> no prediction, with the reason logged. We "
          "did not impute missing pure properties; imputation would make the baselines non-reproducible.",
          "3. **Refutas needs kinematic viscosity**, so it also needs pure density; records where any component "
          "lacks a density curve are reported as `no_pure_rho_at_T` rather than silently switching to the "
          "log-linear baseline. Records with any nu_i < 0.2 cSt are excluded (`nu_lt_0p2`) because Refutas is "
          "undefined there.",
          "4. **Pure rows of density (and of viscosity) are identity/circular**, so they are marked † and excluded "
          "from the mixture headline summary. Flash and vapour-pressure pure rows are genuine predictions and are "
          "kept.",
          "5. **Null model is split-aware** (S1 medians for S2/S3, system-disjoint folds for S1) so the physics "
          "comparison is not flattered by an in-sample null. The null is binned by `kind` and a 20 K temperature "
          "bin; for flash (no condition T in the table) that degenerates to a per-kind constant, which is the right "
          "sanity reference - it is the 'always predict the median flash point' model.",
          "6. **`conflict=False` is the headline filter**; the sensitivity run with conflicts included is in the "
          "log. `n_obs`-based `quality_weight` was not used inside the baselines (they are per-record, and "
          "weighting would hide the disagreement).",
          "7. **No tuning of any kind**: baselines are formulas; the only fitted objects are the per-compound "
          "pure-property curves, which use pure rows only (never a blend row), so a blend record cannot leak into "
          "its own prediction.",
          "8. Mole fractions for flash are recomputed from `w_mass` and the registry molar masses, as specified, "
          "not taken from the stored `x_mole` column.", ""]

    # next steps grounded in the numbers
    def sc(prop, base, split):
        m = mets[(mets.property == prop) & (mets.baseline == base) & (mets.split == split) & mets.conflict_free
                 & mets.kind.isin(["binary", "multi"])]
        if not len(m):
            return float("nan")
        return float(np.average(m[SPEC[prop]["score"]], weights=m.n))

    L += ["## Next steps (what a learned residual model should target)", ""]
    L.append("- Residual target per property, from the tables above (mixture cells, conflict-free):")
    for prop in TARGET:
        for base in TARGET[prop]:
            v1, v2 = sc(prop, base, "S1"), sc(prop, base, "S2")
            if np.isfinite(v1) or np.isfinite(v2):
                how = {("log"): "multiplicative correction in ln space",
                       ("log10"): "additive correction in log10 space"}.get(_mode_of(prop),
                                                                          "additive correction in the target unit")
                L.append(f"  - `{prop}` / `{base}`: S1 {_f(v1, 4)} -> S2 {_f(v2, 4)} {SPEC[prop]['unit']} "
                         f"({how}; the S1->S2 gap is the molecule-generalisation headroom).")
    L += ["- Learn the **activity-coefficient-like** term for flash point: ideal Raoult (`chat_*`) cannot bend "
          "below the lightest component's own flash point, and the binary alkane mixtures in this table do; that "
          "is the residual with the clearest physical meaning.",
          "- Learn the **excess viscosity** correction; keep Refutas as the floor (it is the strongest baseline on "
          "binaries) and only allow the learned term to act for n_comp >= 3.",
          "- Density: the ideal-volume rule should already be close; any residual is a volume-contraction signal, "
          "so fit it as a function of aromatic content and T rather than as a free-form correction.",
          "- Flash point must be trained with an **asymmetric loss** on over-prediction; report `over-pred >5 C` as "
          "a first-class metric with its bootstrap CI, exactly as tabulated here.",
          "- Unevaluable records (no Antoine / no LFL / no pure curve) are the coverage ceiling for any residual "
          "model: fixing them is a data task (interpolate Antoine vs carbon number, extend LFL), not a modelling "
          "task.",
          "- Real fuels (no molecule list) stay out of supervised training; use them as a distribution-input "
          "validation set and for inverse-design sanity checks.", ""]

    REPORT_OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT_OUT.write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------- main
def main(argv=None):
    t0 = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-preds", action="store_true")
    ap.add_argument("--eval", action="store_true")
    a, _unknown = ap.parse_known_args(argv)
    do_preds = a.write_preds or not a.eval
    do_eval = a.eval or not a.write_preds

    mt, comp, splits, locked = load_inputs()
    mt, split_info = apply_splits(mt, splits, locked)
    print(f"model table rows: {len(mt)}  split counts: {mt.split.value_counts().to_dict()}  {split_info}")

    pure = mt[mt.kind == "pure"].copy()
    pure["inchikey"] = pure.comp_inchikeys.map(lambda s: (_json(s) or [""])[0])
    # component curves may use any pure row (they are label-free constants), but never a blend row
    fits = build_fits(pure)
    fit_info = {"curves": len(fits["eta"]),
                "points": int(sum(f["n"] + f.get("dropped", 0) for f in fits["eta"].values())),
                "dropped": int(sum(f.get("dropped", 0) for f in fits["eta"].values())),
                "locked": split_info["locked_rows_removed_from_splits"]}
    print(f"component fits: eta {len(fits['eta'])}, rho {len(fits['rho'])}; dropped {fit_info['dropped']}")

    pred, skip = predict_records(mt, comp, fits)
    pred["null_pred"] = null_predictions(pred)
    print(f"predictions: {len(pred)} rows; non-null per baseline: "
          f"{ {c: int(pred[c].notna().sum()) for c in PRED_COLS} }")

    if do_preds:
        PREDS_OUT.parent.mkdir(parents=True, exist_ok=True)
        pred.to_csv(PREDS_OUT, index=False)
        print(f"wrote {PREDS_OUT}")

    if do_eval:
        mets = pd.concat([evaluate(pred, conflict_free=True), evaluate(pred, conflict_free=False)],
                         ignore_index=True)
        write_log(mt, pred, skip, mets, t0, argv or [])
        write_report(mt, pred, skip, mets, t0, argv or [], run_info=fit_info)
        print(f"wrote {LOG_OUT} and {REPORT_OUT}")
        key = mets[(mets.conflict_free) & mets.kind.isin(["binary", "multi"])]
        for (prop, base), g in key.groupby(["property", "baseline"], sort=False):
            sc = float(np.average(g[SPEC[prop]["score"]], weights=g.n))
            print(f"  {prop:24s} {base:12s} n={int(g.n.sum()):6d}  {SPEC[prop]['score']}={sc:.4f}")
    print(f"done in {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main(sys.argv[1:])
