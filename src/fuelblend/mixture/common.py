"""
fuelblend.mixture.common
========================
Shared building blocks for T8 (DeepSets mixture model) and T9 (inverse design).

  * CompoundRegistry      static descriptors + Antoine/LFL for every pure compound
  * PurePropertyProvider  pure-compound ln(eta) and rho vs T, grid-cached on disk
  * physics baselines     ideal-volume density, Refutas viscosity, Le Chatelier flash
                          (torch, vectorised, differentiable -> usable inside T9)
  * assemble_batch        ONE code path that turns raw tensors into model inputs,
                          used by both training (T8) and optimisation (T9)

Install at: fuelblend_data/src/fuelblend/mixture/common.py  (+ empty __init__.py)

ASSUMPTIONS (marked ADAPTER where you may need to edit):
  [ADAPTER-1] pure eta(T), rho(T) come from fuelblend.api.predict_blend.predict() called on a
              single-component input [(inchikey, 1.0)]. Slow on the first run, cached afterwards.
              If you have a direct pure-GBM function, pass it as PurePropertyProvider(fetch=...).
  [ADAPTER-2] baselines are re-implemented here (Refutas / ideal volume / Le Chatelier with
              Raoult + Antoine + LFL). Check that baseline numbers printed by T8 on S2 match your
              WP3 numbers; if not, swap the three functions below for your own.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch

LN10 = math.log(10.0)
P_ATM_KPA = 101.325
T_LO, T_HI = 230.0, 470.0          # API validity window (K)
PROPS = ("viscosity_dyn_mPas", "density_kg_m3", "flash_point_C")
CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]


# --------------------------------------------------------------------------- registry
class CompoundRegistry:
    """Static (T-independent) per-compound arrays, indexed 0..n-1 in compounds.csv order."""

    def __init__(self, compounds_csv):
        from fuelblend.features.pure_features import pure_feature_row

        df = pd.read_csv(compounds_csv).drop_duplicates("inchikey").reset_index(drop=True)
        n = len(df)
        self.ik = df["inchikey"].tolist()
        self.ik2i = {k: i for i, k in enumerate(self.ik)}
        self.smiles = df["smiles"].tolist()

        feats = np.full((n, 16), np.nan)          # 11 RDKit + 5 class one-hot (no T terms)
        for i in range(n):
            try:
                feats[i] = pure_feature_row(self.ik[i], self.smiles[i], 298.15,
                                            df["cls"].iat[i], CLASSES)[:16]
            except Exception:
                pass                                # e.g. cls == 'other' -> compound unusable

        mw = df["mw"].to_numpy(float)
        mw = np.where(np.isfinite(mw), mw, feats[:, 0])
        A, B, C = (df[c].to_numpy(float) for c in ("antoine_A", "antoine_B", "antoine_C"))
        lfl = df["lfl_vol_pct"].to_numpy(float) / 100.0
        has = np.isfinite(A) & np.isfinite(B) & np.isfinite(C) & np.isfinite(lfl) & (lfl > 0)

        self.has_flash = has
        self.A, self.B, self.C = (np.where(has, v, 0.0) for v in (A, B, C))
        self.lfl = np.where(has, lfl, 1.0)         # volume fraction
        self.mw = np.where(np.isfinite(mw), mw, 100.0)

        def lnp(T):                                # ln(P_sat / kPa)
            return np.clip(np.where(has, LN10 * (self.A - self.B / (T + self.C)), 0.0), -12, 8)

        raw = np.column_stack([feats, lnp(300.0), lnp(400.0),
                               np.where(has, np.log(self.lfl), 0.0), has.astype(float)])
        self.usable = np.isfinite(feats).all(1) & np.isfinite(mw)
        mu = raw[self.usable].mean(0)
        sd = raw[self.usable].std(0)
        sd[sd < 1e-6] = 1.0
        S = (raw - mu) / sd
        S[~self.usable] = 0.0
        self.S = S.astype(np.float32)
        self.d_static = S.shape[1]                 # 20


class RegistryTensors:
    def __init__(self, reg: CompoundRegistry, device):
        t = lambda a, dt=torch.float32: torch.as_tensor(np.asarray(a), dtype=dt, device=device)
        self.S = t(reg.S)
        self.mw = t(np.nan_to_num(reg.mw, nan=100.0))
        self.A, self.B, self.C, self.lfl = t(reg.A), t(reg.B), t(reg.C), t(reg.lfl)
        self.has_flash = t(reg.has_flash, torch.bool)
        self.device = device


# --------------------------------------------------------------------------- pure-property cache
_MEASURED_FITS = None


def _get_measured_fits():
    global _MEASURED_FITS
    if _MEASURED_FITS is not None:
        return _MEASURED_FITS
    mt_p = Path(__file__).resolve().parents[3] / "data/processed/model_table.csv"
    if not mt_p.exists():
        _MEASURED_FITS = ({}, {})
        return _MEASURED_FITS
    import json
    mt = pd.read_csv(mt_p, low_memory=False)
    pures = mt[mt["n_comp"] == 1].copy()
    pures["ik"] = pures["comp_inchikeys"].apply(lambda s: json.loads(s)[0] if s.startswith("[") else s)
    fr, fe = {}, {}
    for ik, g in pures[pures["property"] == "density_kg_m3"].groupby("ik"):
        T, y = g["T_K"].values, g["value"].values
        if len(T) == 1:
            fr[ik] = (float(y[0]), 0.0, float(T[0]), float(T[0]))
        else:
            b, a = np.polyfit(T, y, 1)
            fr[ik] = (float(a), float(b), float(T.min()), float(T.max()))
    for ik, g in pures[(pures["property"] == "viscosity_dyn_mPas") & (pures["value"] > 0)].groupby("ik"):
        T, y = g["T_K"].values, np.log(g["value"].values)
        if len(T) == 1:
            fe[ik] = (float(y[0]), 0.0, float(T[0]), float(T[0]))
        else:
            b, a = np.polyfit(np.log(T), y, 1)
            fe[ik] = (float(a), float(b), float(T.min()), float(T.max()))
    _MEASURED_FITS = (fr, fe)
    return _MEASURED_FITS


def _fetch_via_api(inchikey: str, T_K: float):
    fr, fe = _get_measured_fits()
    T = float(T_K)
    d = None
    if inchikey in fr:
        a, b, tlo, thi = fr[inchikey]
        if tlo - 25.0 <= T <= thi + 25.0:
            val = a + b * T
            if val > 0:
                d = float(val)
    le = None
    if inchikey in fe:
        a, b, tlo, thi = fe[inchikey]
        if tlo - 25.0 <= T <= thi + 25.0:
            le = float(a + b * math.log(T))
    if d is None or le is None:
        from fuelblend.api.predict_blend import _predict_viscosity, _predict_density          # [ADAPTER-1]
        if d is None:
            rd = _predict_density([inchikey], [1.0], [1.0], T)
            dv = rd.get("value")
            d = float(dv) if dv else float("nan")
        if le is None:
            rv = _predict_viscosity([inchikey], [1.0], [1.0], T)
            vv = rv.get("value")
            le = math.log(vv) if vv and vv > 0 else float("nan")
    return le, d


class PurePropertyProvider:
    """ln(eta_dyn / mPa.s) and rho (kg/m3) of each pure compound on a T grid, with interpolation
    (linear in 1/T for ln eta, linear in T for rho). Results are cached in an .npz file."""

    def __init__(self, reg: CompoundRegistry, cache_path, T_grid=None, fetch=None, verbose=True):
        self.reg, self.verbose = reg, verbose
        self.T = np.linspace(T_LO, T_HI, 13) if T_grid is None else np.asarray(T_grid, float)
        self.fetch = fetch or _fetch_via_api
        self.cache_path = Path(cache_path)
        n = len(reg.ik)
        self.ln_eta = np.full((n, len(self.T)), np.nan)
        self.rho = np.full((n, len(self.T)), np.nan)
        self.done = np.zeros(n, bool)
        if self.cache_path.exists():
            z = np.load(self.cache_path, allow_pickle=False)
            if z["T"].shape == self.T.shape and np.allclose(z["T"], self.T):
                for k, ik in enumerate(z["ik"]):
                    j = reg.ik2i.get(str(ik))
                    if j is not None and z["done"][k]:
                        self.ln_eta[j], self.rho[j], self.done[j] = z["ln_eta"][k], z["rho"][k], True
        self._refresh_valid()

    def _refresh_valid(self):
        self.valid = np.isfinite(self.ln_eta).all(1) & np.isfinite(self.rho).all(1)

    def _eval_one(self, i):
        le = np.full(len(self.T), np.nan)
        rh = np.full(len(self.T), np.nan)
        for k, T in enumerate(self.T):
            try:
                le[k], rh[k] = self.fetch(self.reg.ik[i], T)
            except Exception:
                pass
        for arr in (le, rh):                    # patch isolated gaps by interpolation in T
            ok = np.isfinite(arr)
            if 2 <= ok.sum() < len(arr):
                arr[~ok] = np.interp(self.T[~ok], self.T[ok], arr[ok])
        return i, le, rh

    def ensure(self, comp_idx):
        from concurrent.futures import ThreadPoolExecutor
        todo = [int(i) for i in np.unique(comp_idx) if not self.done[int(i)]]
        if not todo:
            return
        with ThreadPoolExecutor(max_workers=4) as ex:
            for n, (i, le, rh) in enumerate(ex.map(self._eval_one, todo)):
                if self.verbose and n % 25 == 0:
                    print(f"[pure-cache] {n}/{len(todo)} compounds ...", flush=True)
                self.ln_eta[i], self.rho[i], self.done[i] = le, rh, True
                if (n + 1) % 25 == 0:
                    self.cache_path.parent.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(self.cache_path, ik=np.array(self.reg.ik), T=self.T,
                                        ln_eta=self.ln_eta, rho=self.rho, done=self.done)
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(self.cache_path, ik=np.array(self.reg.ik), T=self.T,
                            ln_eta=self.ln_eta, rho=self.rho, done=self.done)
        self._refresh_valid()

    def lookup(self, idx, T):
        """idx int array (-1 = padding), T broadcastable -> (ln_eta, rho); NaN where invalid."""
        idx = np.asarray(idx)
        safe = np.where(idx < 0, 0, idx)
        T = np.clip(np.broadcast_to(np.asarray(T, float), idx.shape), self.T[0], self.T[-1])
        j = np.clip(np.searchsorted(self.T, T, side="right") - 1, 0, len(self.T) - 2)
        T0, T1 = self.T[j], self.T[j + 1]
        t = (T - T0) / (T1 - T0)
        s = (1 / T - 1 / T0) / (1 / T1 - 1 / T0)
        rho = (1 - t) * self.rho[safe, j] + t * self.rho[safe, j + 1]
        le = (1 - s) * self.ln_eta[safe, j] + s * self.ln_eta[safe, j + 1]
        pad = idx < 0
        return np.where(pad, np.nan, le), np.where(pad, np.nan, rho)


# --------------------------------------------------------------------------- physics (torch)
def density_idvol(w, rho):
    """Ideal volume mixing, w = mass fractions (B,N), rho (B,N) kg/m3 -> (B,)."""
    return 1.0 / (w / rho.clamp_min(1.0)).sum(-1).clamp_min(1e-9)


def refutas_ln_eta(w, ln_eta, rho):
    """Refutas viscosity blending. Returns ln(eta_dyn/mPa.s) of the mixture (B,)."""
    nu = 1000.0 * torch.exp(ln_eta) / rho.clamp_min(1.0)                   # cSt
    vbi = 14.534 * torch.log(torch.log(nu + 0.8).clamp_min(1e-3)) + 10.975
    vbn = (w * vbi).sum(-1)
    nu_mix = torch.exp(torch.exp((vbn - 10.975) / 14.534)) - 0.8
    return torch.log(nu_mix.clamp_min(1e-3) * density_idvol(w, rho) / 1000.0)


def lechatelier_flash_C(x, A, B, C, lfl, n_iter=30):
    """Flash point (degC) from  sum_i x_i P_i(T) / (P_atm * LFL_i) = 1  (Raoult + Antoine + Le Chatelier).
    Newton on ln F(T) with an implicit-function final step, so gradients w.r.t. x are exact."""
    lnx = torch.log(x.clamp_min(1e-30))
    lnpl = torch.log(P_ATM_KPA * lfl)

    def terms(T):
        den = (T[:, None] + C).clamp_min(5.0)
        z = lnx + LN10 * (A - B / den) - lnpl
        dz = LN10 * B / den ** 2
        return z, dz

    with torch.no_grad():
        T = torch.full(x.shape[:1], 300.0, dtype=x.dtype, device=x.device)
        for _ in range(n_iter):
            z, dz = terms(T)
            lse = torch.logsumexp(z, -1)
            d = (torch.softmax(z, -1) * dz).sum(-1).clamp_min(1e-4)
            T = (T - (lse / d).clamp(-40, 40)).clamp(120.0, 700.0)
    z, dz = terms(T)
    d = (torch.softmax(z, -1) * dz).sum(-1).detach().clamp_min(1e-4)
    return T - torch.logsumexp(z, -1) / d - 273.15


def mole_fractions(rt: RegistryTensors, idx, w):
    mask = idx >= 0
    mw = torch.where(mask, rt.mw[idx.clamp_min(0)], torch.ones_like(w))
    x = (w * mask) / mw
    return x / x.sum(-1, keepdim=True).clamp_min(1e-12)


def physics_baselines(rt: RegistryTensors, idx, w, ln_eta_i, rho_i):
    """idx (B,N) registry index, -1 = pad;  w mass fractions (B,N);  pure ln_eta_i, rho_i (B,N) at the
    row temperature (NaN allowed).  Returns base_raw (B,3) = [ln eta, rho, flash degC], valid (B,3), x."""
    mask = idx >= 0
    safe = idx.clamp_min(0)
    w = w * mask
    x = mole_fractions(rt, idx, w)
    vr_ok = ((torch.isfinite(ln_eta_i) & torch.isfinite(rho_i)) | ~mask).all(-1)
    ln_eta_i = torch.where(mask & torch.isfinite(ln_eta_i), ln_eta_i, torch.zeros_like(ln_eta_i))
    rho_i = torch.where(mask & torch.isfinite(rho_i), rho_i, torch.full_like(rho_i, 750.0))
    rho_i = torch.where(mask, rho_i, torch.ones_like(rho_i))
    rho_mix = density_idvol(w, rho_i)
    ln_eta_mix = refutas_ln_eta(w, ln_eta_i, rho_i)
    fl_ok = (rt.has_flash[safe] | ~mask).all(-1)
    flash = lechatelier_flash_C(x, rt.A[safe], rt.B[safe], rt.C[safe], rt.lfl[safe])
    base = torch.stack([ln_eta_mix, rho_mix, flash], -1)
    valid = torch.stack([vr_ok, vr_ok, fl_ok], -1) & torch.isfinite(base)
    return base, valid, x


# --------------------------------------------------------------------------- batch assembly
def ctx_np(T_K, missing):
    T_K = np.asarray(T_K, float)
    return np.stack([(1000.0 / T_K - 3.2) / 0.5, np.asarray(missing, float)], -1).astype(np.float32)


def assemble_batch(rt: RegistryTensors, idx, w, ln_eta_i, rho_i, ctx, base_raw, base_valid):
    """Single path from raw tensors to model inputs (shared by T8 training and T9 optimisation)."""
    mask = idx >= 0
    m = mask.unsqueeze(-1)
    w = w * mask
    x = mole_fractions(rt, idx, w)
    static = rt.S[idx.clamp_min(0)] * m
    pure = torch.stack([ln_eta_i, (rho_i - 750.0) / 100.0], -1) * m
    bn = torch.stack([base_raw[:, 0], (base_raw[:, 1] - 750.0) / 100.0,
                      (base_raw[:, 2] - 60.0) / 50.0], -1)
    bn = torch.where(base_valid, bn, torch.zeros_like(bn))
    return dict(static=static, pure=pure, w=w, x=x, mask=mask, ctx=ctx,
                base_in=torch.cat([bn, base_valid.float()], -1), base_raw=base_raw)
