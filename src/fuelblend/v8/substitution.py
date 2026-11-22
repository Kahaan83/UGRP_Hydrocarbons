"""G1: replace measured pure constants by out-of-fold ML values inside the v5 blend routes and compare accuracy.

One constant kind per property; all components of a system are substituted together. A system whose components do not
all have an ML value is skipped (status `skipped_no_ml_value`); a flash system with a component lacking a measured
Antoine is skipped (status `skipped_no_measured_antoine`) and is never scored with a substitute Antoine (trap 3).

Substitution definitions (documented limits):
  density    rho_i(T) shifted by a constant so that rho_i(293.15 K) equals the ML value (the Track B label basis is 293.15 +/- 5 K,
             so the grid slope converts it to any T, including 288.15 K); ideal volume mixing.
  viscosity  ln eta_i(T) shifted by a constant so that eta_i(293.15 K) equals the ML value (log10 mPa.s). The measured
             T-curve shape is kept: the v7 Andrade slope head is not deployable (it failed against baselines), so this
             arm is a lower bound on the substitution cost.
  flash      LFL_i = 100 * Psat_i(Tf_ML) / P_atm with the measured Antoine; ideal Le Chatelier (designated lab route).
  freezing   ideal Schroeder-van Laar with ML Tm (degC) and dHfus (kJ/mol).
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
TOL = {"density": 2.0, "viscosity": 0.02, "flash": 2.0, "freezing": 2.0}
R_GAS = 8.314462618
P_ATM = 101.325
COLS = ["system_id", "status", "n_rows", "mae_measured", "mae_ml", "bias_ml", "n_substituted", "split_role"]


@lru_cache(maxsize=None)
def _ml(prop: str) -> dict:
    from fuelblend.v8.ml_constants import oof_constants
    d = oof_constants(prop)
    return dict(zip(d.inchikey, d.value))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@lru_cache(maxsize=1)
def _mix():
    from fuelblend.mixture.cold_flow import compound_meta
    from fuelblend.mixture.common import CLASSES, CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    _, cls, nc = compound_meta(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    return reg, prov, cls, nc, CLASSES


class ShiftedProvider:
    """Wraps a PurePropertyProvider; adds constant per-compound shifts to ln eta and rho."""

    def __init__(self, base, d_ln_eta: dict, d_rho: dict):
        self.base, self.d_le, self.d_rho = base, d_ln_eta, d_rho
        self.valid, self.T = base.valid, base.T

    def lookup(self, idx, T):
        le, rho = self.base.lookup(idx, T)
        le, rho = le.copy(), rho.copy()
        a = np.asarray(idx)
        for j in np.unique(a[a >= 0]):
            if int(j) in self.d_le:
                le[a == j] += self.d_le[int(j)]
            if int(j) in self.d_rho:
                rho[a == j] += self.d_rho[int(j)]
        return le, rho


class _RefutasShim:
    """Same formula as scripts/v4/build_viscosity_leaderboard_v4.py::refutas_nu (that module cannot be imported: it needs the
    removed scripts/v3 file). Mass-fraction VBN blending, ideal-volume density."""

    @staticmethod
    def refutas_nu(provider, idx, w, T_K):
        ln_eta, rho = provider.lookup(idx[None, :], T_K)
        ln_eta, rho = ln_eta[0], rho[0]
        nu_i = 1000.0 * np.exp(ln_eta) / rho
        vbn_i = 14.534 * np.log(np.log(nu_i + 0.8)) + 10.975
        vbn = float((w * vbn_i).sum())
        nu = float(np.exp(np.exp((vbn - 10.975) / 14.534)) - 0.8)
        return nu, float(1.0 / (w / rho).sum()), ln_eta, rho


def _shifts(prop: str, override):
    reg, prov, *_ = _mix()
    d_le, d_rho = {}, {}
    ml_rho, ml_v = _ml("rho15"), _ml("logvisc293")
    for ik, j in reg.ik2i.items():
        if not prov.valid[j]:
            continue
        if prop == "density":
            g = prov.lookup(np.array([[j]]), 293.15)[1][0][0]
            m = g if override == "measured" else ml_rho.get(ik, np.nan)
            if np.isfinite(m):
                d_rho[j] = float(m - g)
        elif prop == "viscosity":
            g = prov.lookup(np.array([[j]]), 293.15)[0][0][0]
            m = g if override == "measured" else (math.log(10.0) * ml_v[ik] if ik in ml_v else np.nan)
            if np.isfinite(m):
                d_le[j] = float(m - g)
    return d_le, d_rho


def _frame(rows: list) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=COLS)


# ----------------------------------------------------------------------------------------------- density
def _density(set_name: str, override) -> pd.DataFrame:
    reg, prov, *_ = _mix()
    roles = pd.read_csv(ROOT / "results/v4/L/density_lab_residuals.csv", usecols=["system_id", "split_role"]) \
        .drop_duplicates().set_index("system_id").split_role.to_dict()
    m = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", low_memory=False)
    m = m[m.kind.isin(["binary", "multi"]) & m.is_hc.astype(bool) & (m.property == "density_kg_m3") & ~m.conflict.astype(bool)]
    _, d_rho = _shifts("density", override)
    sp = ShiftedProvider(prov, {}, d_rho)
    rows = []
    for sid, g in m.groupby("system_id"):
        if sid not in roles or (set_name and roles[sid] != set_name):
            continue
        ks = json.loads(g.iloc[0].comp_inchikeys)
        ii = [reg.ik2i.get(k) for k in ks]
        if any(i is None or not prov.valid[i] for i in ii):
            continue
        ns = sum(i in d_rho for i in ii)
        if ns < len(ii):
            rows.append(dict(system_id=sid, status="skipped_no_ml_value", n_rows=len(g), n_substituted=ns, split_role=roles[sid]))
            continue
        e0, e1 = [], []
        ia = np.array(ii)
        for r in g.itertuples():
            w = np.array(json.loads(r.w_mass), float)
            if not np.all(np.isfinite(w)) or w.sum() <= 0:
                continue
            w = w / w.sum()
            for P, acc in ((prov, e0), (sp, e1)):
                rho = P.lookup(ia[None, :], float(r.T_K))[1][0]
                acc.append(1.0 / (w / rho).sum() - float(r.value))
        if e0:
            rows.append(dict(system_id=sid, status="ok", n_rows=len(e0), mae_measured=np.mean(np.abs(e0)),
                             mae_ml=np.mean(np.abs(e1)), bias_ml=np.mean(e1), n_substituted=ns, split_role=roles[sid]))
    return _frame(rows)


# ----------------------------------------------------------------------------------------------- viscosity
def _viscosity(set_name: str, override) -> pd.DataFrame:
    reg, prov, cls, nc, CLASSES = _mix()
    s1 = _load("lab_v5", "scripts/v5/eval_lab_known_composition_v5.py")
    lb = _RefutasShim()
    meas = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", low_memory=False)
    roles = pd.read_csv(ROOT / "results/v5/S1/viscosity_lab_residuals_v5.csv", usecols=["system_id", "split_role"]) \
        .drop_duplicates().set_index("system_id").split_role.to_dict()
    sel = s1.lab_selection(meas)
    sel = sel[sel.system_id.isin([s for s, r in roles.items() if not set_name or r == set_name])]
    d_le, _ = _shifts("viscosity", override)
    sp = ShiftedProvider(prov, d_le, {})
    base = s1.score_rows(sel, reg, prov, cls, nc, lb, CLASSES)
    base = base[base.method == "refutas"]
    ml = s1.score_rows(sel, reg, sp, cls, nc, lb, CLASSES)
    ml = ml[ml.method == "refutas"]
    sys_iks = sel.groupby("system_id").comp_inchikeys.first().map(json.loads).to_dict()
    rows = []
    for sid, g in base.groupby("system_id"):
        ns = sum(reg.ik2i.get(k) in d_le for k in sys_iks[sid])
        if ns < len(sys_iks[sid]):
            rows.append(dict(system_id=sid, status="skipped_no_ml_value", n_rows=len(g), n_substituted=ns, split_role=roles[sid]))
            continue
        h = ml[ml.system_id == sid]
        if len(h) != len(g):                              # substituted Refutas undefined (VBN needs nu + 0.8 > 1) for some rows
            rows.append(dict(system_id=sid, status="ml_route_failed", n_rows=len(g), n_substituted=ns, split_role=roles[sid]))
            continue
        rows.append(dict(system_id=sid, status="ok", n_rows=len(g), mae_measured=g.resid.abs().mean(),
                         mae_ml=h.resid.abs().mean(), bias_ml=h.resid.mean(), n_substituted=ns, split_role=roles[sid]))
    return _frame(rows)


# ----------------------------------------------------------------------------------------------- flash
def _lfl_from_flash(row, tf_c: float) -> float:
    T = tf_c + 273.15
    p = 10 ** (float(row.antoine_A) - float(row.antoine_B) / (T + float(row.antoine_C)))
    return 100.0 * p / P_ATM


def _flash(set_name: str, override) -> pd.DataFrame:
    from rdkit import Chem
    from fuelblend.baselines import flash_unifac as fu
    from fuelblend.v8.ml_constants import _pool
    if set_name == "locked":
        return _frame([])                                 # no locked lab flash systems exist (all 31 are dev)
    f = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    reg = fu._default_registry()
    ml = _ml("flash_pure")
    comp = _pool()[0]
    smi2ik = {Chem.MolToSmiles(Chem.MolFromSmiles(s)): ik for ik, s in zip(comp.inchikey, comp.smiles)}
    rows, acc = [], {}
    for r in f.itertuples():
        smi = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in json.loads(r.comp_smiles)]
        x = json.loads(r.x_mole)
        if any(s not in reg.index or pd.isna(reg.loc[s, "antoine_A"]) for s in smi):
            rows.append(dict(system_id=r.system_id, status="skipped_no_measured_antoine", n_rows=1, split_role="dev"))
            continue
        iks = [smi2ik.get(s) for s in smi]
        if override != "measured" and any(k is None or k not in ml for k in iks):
            rows.append(dict(system_id=r.system_id, status="skipped_no_ml_value", n_rows=1, split_role="dev"))
            continue
        reg2 = reg.copy()
        if override != "measured":
            for s, k in zip(smi, iks):
                reg2.loc[s, "lfl_vol_pct"] = _lfl_from_flash(reg2.loc[s], float(ml[k]))
        a = fu.flash_point_unifac_detail(smi, x, registry=reg, use_unifac=False)["flash_C"]
        b = fu.flash_point_unifac_detail(smi, x, registry=reg2, use_unifac=False)["flash_C"]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        y = float(r.measured_C)
        acc.setdefault(r.system_id, []).append((a - y, b - y, len(smi)))
    for sid, v in acc.items():
        e0, e1, n = zip(*v)
        rows.append(dict(system_id=sid, status="ok", n_rows=len(v), mae_measured=np.mean(np.abs(e0)),
                         mae_ml=np.mean(np.abs(e1)), bias_ml=np.mean(e1), n_substituted=n[0], split_role="dev"))
    out = _frame(rows)
    ok = set(out[out.status == "ok"].system_id)
    return out[(out.status == "ok") | ~out.system_id.isin(ok)].drop_duplicates("system_id").reset_index(drop=True)


# ----------------------------------------------------------------------------------------------- freezing
def _freezing(set_name: str, override) -> pd.DataFrame:
    if set_name == "locked":
        return _frame([])                                 # SLE lab systems have no locked role
    d = pd.read_csv(ROOT / "results/v5/SLE2/sle_eval_extended.csv")
    d = d[d.scored.astype(bool) & d.ideal_v4_verified_K.notna()]
    tm, dh = _ml("Tm"), _ml("dHfus")
    consts = pd.read_csv(ROOT / "results/v5/SLE2/constants_extended.csv").set_index("inchikey")
    rows = []
    for sid, g in d.groupby("system_id"):
        iks = json.loads(g.inchikeys.iloc[0])
        if override == "measured":
            tabs = {k: (consts.loc[k, "tm_K"], consts.loc[k, "dh_J_mol"]) for k in iks}
        elif all(k in tm and k in dh for k in iks):
            tabs = {k: (tm[k] + 273.15, dh[k] * 1000.0) for k in iks}
        else:
            rows.append(dict(system_id=sid, status="skipped_no_ml_value", n_rows=len(g), split_role="dev"))
            continue
        e0, e1 = [], []
        for r in g.itertuples():
            x = json.loads(r.x_mole)
            ts = [1.0 / (1.0 / tabs[k][0] - R_GAS / tabs[k][1] * math.log(xi)) for xi, k in zip(x, iks) if xi > 0]
            if not ts or any((not np.isfinite(t)) or t <= 0 for t in ts):
                continue
            e0.append(r.ideal_v4_verified_K - r.T_obs_K)
            e1.append(max(ts) - r.T_obs_K)
        if e0:
            rows.append(dict(system_id=sid, status="ok", n_rows=len(e0), mae_measured=np.mean(np.abs(e0)),
                             mae_ml=np.mean(np.abs(e1)), bias_ml=np.mean(e1), n_substituted=len(iks), split_role="dev"))
    return _frame(rows)


_FN = {"density": _density, "viscosity": _viscosity, "flash": _flash, "freezing": _freezing}


def score_substitution(prop: str, set_name: str = "dev", _ml_override=None) -> pd.DataFrame:
    """One row per system: system_id, status, n_rows, mae_measured, mae_ml, bias_ml, n_substituted, split_role."""
    return _FN[prop](set_name, _ml_override)


def summarize(df: pd.DataFrame, prop: str, n_boot: int = 2000, seed: int = 42, alpha: float = 0.10) -> dict:
    """Paired group bootstrap of dMAE = mean over systems of (mae_ml - mae_measured); upper end of the one-sided 90 % bound."""
    ok = df[df.status == "ok"]
    d = (ok.mae_ml - ok.mae_measured).to_numpy(float)
    out = dict(prop=prop, n_systems=len(ok), n_skipped=int((df.status != "ok").sum()), tol=TOL[prop])
    if len(d) == 0:
        return dict(out, d_mae=np.nan, ci_lo=np.nan, ci_hi=np.nan, verdict="no_data")
    rng = np.random.default_rng(seed)
    bs = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(n_boot)])
    hi, lo = float(np.quantile(bs, 1 - alpha)), float(np.quantile(bs, alpha))
    return dict(out, d_mae=float(d.mean()), ci_lo=lo, ci_hi=hi, verdict="pass" if hi <= TOL[prop] else "fail")
