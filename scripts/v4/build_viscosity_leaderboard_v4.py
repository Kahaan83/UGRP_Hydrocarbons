"""v4 Task 11: viscosity leaderboard on the 57 CRC 647 fuels with fuel-grouped NESTED CV.

Protocol (headline = leave-one-fuel-out):
  * dev fuels  = CRC647 fuels NOT in splits/v4/real_fuels_v4_splits.json["locked"]  (45)
  * locked     = the 12 locked CRC647 fuels: NEVER used for fitting, selection or tuning of any model; reported
                 as a separate untouched evaluation of models trained on the dev fuels only.
  * dev rows are LOO by fuel over the dev fuels; hyper-parameters (bilinear rank & lambda, ridge alpha, GP
    length-scale & noise) are chosen by INNER grouped CV inside the training fold only (R6).
  * Measured kinematic viscosity: jet_properties_v4.csv rows with provenance raw_file (CRC647, D455, 20/-20/-40 C).
  * Compositions: compositions_v4.csv, provenance 'reconstructed' (formula-derived from measured aromatics + D86):
    all errors mix reconstruction error and physics error.
Outputs: results/v4/V/viscosity_leaderboard_v4.csv, viscosity_predictions_v4.csv, selected_hyperparameters_v4.csv
"""
from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.mixture.bilinear_pairwise import BilinearPairwiseInteraction  # noqa: E402
from fuelblend.mixture.cold_flow import compound_meta  # noqa: E402
from fuelblend.mixture.common import CLASSES, CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.mixture.structured_residual import StructuredResidual, nested_loo, select_by_inner_cv  # noqa: E402
from fuelblend.thermo.walther_viscosity import fit_walther_parameters, walther_viscosity  # noqa: E402

OUT = ROOT / "results" / "v4" / "V"
COMP = ROOT / "data" / "real_fuels" / "v4" / "compositions_v4.csv"
PROPS = ROOT / "data" / "real_fuels" / "v4" / "jet_properties_v4.csv"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
MEAS = ROOT / "data" / "unified" / "measurements_dedup.csv"
SEED = 20261110
TEMPS_C = (20.0, -20.0, -40.0)
FEATS = [f"frac_{c}" for c in CLASSES] + ["nC_mean", "nC_std"]
ALPHA_GRID = [1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3]
GP_GRID = [(ls, nr) for ls in (2.0, 4.0, 8.0) for nr in (0.1, 0.3, 1.0)]
BIL_GRID = [(r, wd) for r in (2, 4, 8) for wd in (1e-4, 1e-3, 1e-2)]
LAB_T_RANGE = (253.0, 340.0)


# ------------------------------------------------------------------------ existing bootstrap helper
def _load_cluster_bootstrap_ci():
    p = ROOT / "scripts" / "v3" / "build_viscosity_leaderboard.py"
    spec = importlib.util.spec_from_file_location("v3_visc_leaderboard", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.cluster_bootstrap_ci


cluster_bootstrap_ci = _load_cluster_bootstrap_ci()


def rmse_log(e):
    return float(np.sqrt(np.mean(np.square(e))))


def mae_log(e):
    return float(np.mean(np.abs(e)))


def paired_bootstrap(err_a, err_b, clusters, metric, n_resamples=2000, seed=SEED):
    """Paired cluster bootstrap of metric(a) - metric(b) (same fuel resamples for both models)."""
    uc, inv = np.unique(clusters, return_inverse=True)
    idx = [np.where(inv == c)[0] for c in range(len(uc))]
    rng = np.random.default_rng(seed)
    d = np.empty(n_resamples)
    for b in range(n_resamples):
        s = np.concatenate([idx[c] for c in rng.integers(0, len(uc), len(uc))])
        d[b] = metric(err_a[s]) - metric(err_b[s])
    return float(metric(err_a) - metric(err_b)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


# ------------------------------------------------------------------------ composition -> representatives
class Mapper:
    """Maps a composition component to a registry compound with valid pure eta(T), rho(T) of the same
    (class, carbon number) bin (nearest within 3 carbons), like fuel_hindcast.Reps.  Perdeuterated entries use
    the non-deuterated parent's class / carbon number."""

    def __init__(self, reg, cls, nc, provider):
        self.reg, self.cls, self.nc = reg, cls, nc
        self.by = {}
        self.valid_idx = set()
        for i in range(len(reg.ik)):
            if reg.usable[i] and provider.valid[i] and cls[i] in CLASSES and np.isfinite(nc[i]):
                self.valid_idx.add(i)
                self.by.setdefault((cls[i], int(nc[i])), []).append(i)
        self.ns = {c: sorted(n for (cc, n) in self.by if cc == c) for c in CLASSES}

    def class_nc(self, ik, smiles):
        j = self.reg.ik2i.get(ik)
        if "[2H]" in str(smiles):                             # perdeuterated -> parent connectivity
            j = self.reg.ik2i.get(ik.split("-")[0] + "-UHFFFAOYSA-N", j)
        if j is None or self.cls[j] not in CLASSES or not np.isfinite(self.nc[j]):
            return None
        return self.cls[j], int(self.nc[j])

    def rep(self, c, n, max_gap=3):
        ns = self.ns.get(c, [])
        if not ns:
            return None
        b = min(ns, key=lambda k: (abs(k - n), k))
        return self.by[(c, b)][0] if abs(b - n) <= max_gap else None


def direct_idx(mp, ik):
    """r2 path: the row's OWN registry compound (by InChIKey) if usable with valid pure eta/rho, else None."""
    j = mp.reg.ik2i.get(ik)
    return None if j is None or j not in mp.valid_idx else j


def fuel_composition(df_f, mp, direct=False):
    """-> (idx, w, descriptors dict, dropped_mass_fraction) for one fuel; lump removed, renormalised.
    direct=True (provenance reconstructed_r2): no re-binning to the first compound of a (class, nC) bin; invalid -> dropped+counted."""
    df_f = df_f[~df_f.is_lumped.astype(bool) & (df_f.inchikey != "LUMP_NO_REALISTIC_REP") & (df_f.fraction > 0)]
    agg, cn, dropped = {}, [], 0.0
    for _, r in df_f.iterrows():
        if direct:
            i = direct_idx(mp, r.inchikey)
            k = (mp.cls[i], int(mp.nc[i])) if i is not None else None
        else:
            k = mp.class_nc(r.inchikey, r.smiles)
            i = mp.rep(*k) if k else None
        if k is None or i is None:
            dropped += r.fraction
            continue
        agg[i] = agg.get(i, 0.0) + r.fraction
        cn.append((k[0], k[1], r.fraction))
    tot = sum(agg.values())
    idx = np.array(list(agg), dtype=np.int64)
    w = np.array([agg[i] for i in idx]) / tot
    cw = np.array([c[2] for c in cn]); cw = cw / cw.sum()
    ncs = np.array([c[1] for c in cn], float)
    d = {f"frac_{c}": float(cw[[x[0] == c for x in cn]].sum()) for c in CLASSES}
    m = float((cw * ncs).sum())
    d["nC_mean"], d["nC_std"] = m, float(np.sqrt((cw * (ncs - m) ** 2).sum()))
    return idx, w, d, dropped / (dropped + tot)


# ------------------------------------------------------------------------ physics
def refutas_nu(provider, idx, w, T_K):
    """Kinematic Refutas (mass-fraction VBN blending, cSt) + ideal-volume density; plus component data."""
    ln_eta, rho = provider.lookup(idx[None, :], T_K)
    ln_eta, rho = ln_eta[0], rho[0]
    nu_i = 1000.0 * np.exp(ln_eta) / rho
    vbn_i = 14.534 * np.log(np.log(nu_i + 0.8)) + 10.975
    vbn = float((w * vbn_i).sum())
    nu = float(np.exp(np.exp((vbn - 10.975) / 14.534)) - 0.8)
    rho_mix = float(1.0 / (w / rho).sum())
    return nu, rho_mix, ln_eta, rho


# ------------------------------------------------------------------------ bilinear
def pad_batch(rows, d):
    n = max(r["w"].shape[0] for r in rows)
    z = torch.zeros(len(rows), n, d)
    w = torch.zeros(len(rows), n)
    for b, r in enumerate(rows):
        k = r["w"].shape[0]
        z[b, :k], w[b, :k] = r["z"], r["w"]
    return z, w


def bilinear_fit_predict_factory(rows, y_res):
    d = rows[0]["z"].shape[1]

    def fit_predict(tr, te, hp):
        rank, wd = hp
        z_tr, w_tr = pad_batch([rows[i] for i in tr], d)
        z_te, w_te = pad_batch([rows[i] for i in te], d)
        y_tr = torch.tensor(y_res[tr], dtype=torch.float32).unsqueeze(1)
        model = BilinearPairwiseInteraction(in_features=d, rank=rank, num_tasks=1, seed=0)
        opt = torch.optim.Adam(model.parameters(), lr=0.03, weight_decay=wd)
        for _ in range(40):
            opt.zero_grad()
            nn.functional.mse_loss(model(z_tr, w_tr), y_tr).backward()
            opt.step()
        with torch.no_grad():
            return model(z_te, w_te).squeeze(1).numpy().astype(float)
    return fit_predict


# ------------------------------------------------------------------------ lab training data
def build_lab(reg, provider, cls, nc):
    m = pd.read_csv(MEAS, low_memory=False)
    w = m[m.kind.isin(["binary", "multi"]) & m.is_hc.astype(bool) & (m.property == "viscosity_dyn_mPas")
          & ~m.conflict.astype(bool) & m["T_K"].between(*LAB_T_RANGE)]
    rows = []
    for r in w.itertuples():
        ks = json.loads(r.comp_inchikeys); ws = np.array(json.loads(r.w_mass), float)
        ii = [reg.ik2i.get(k) for k in ks]
        if any(i is None or not provider.valid[i] or not reg.usable[i] or cls[i] not in CLASSES
               or not np.isfinite(nc[i]) for i in ii) or not np.all(np.isfinite(ws)) or ws.sum() <= 0:
            continue
        ii = np.array(ii); ws = ws / ws.sum()
        nu, rho_mix, _, _ = refutas_nu(provider, ii, ws, r.T_K)
        d = {f"frac_{c}": float(sum(wi for i, wi in zip(ii, ws) if cls[i] == c)) for c in CLASSES}
        n_ = np.array([nc[i] for i in ii]); m_ = float((ws * n_).sum())
        d["nC_mean"], d["nC_std"] = m_, float(np.sqrt((ws * (n_ - m_) ** 2).sum()))
        # n_comp / measured_mPas / raw_refutas_mPas are additive provenance columns (used by eval_lab_known_composition.py);
        # the leaderboard itself selects only FEATS, T_K, resid and system_id from this frame.
        rows.append(dict(system_id=r.system_id, T_K=r.T_K, resid=np.log(r.value) - np.log(nu * rho_mix / 1000.0),
                         n_comp=len(ks), measured_mPas=float(r.value), raw_refutas_mPas=float(nu * rho_mix / 1000.0), **d))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------ main
def main(argv=None):
    global COMP, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--compositions", type=Path, default=COMP, help="composition table (default: the v4 file)")
    ap.add_argument("--out", type=Path, default=OUT, help="output directory (default: results/v4/V)")
    a = ap.parse_args(argv)
    COMP, OUT = a.compositions, a.out
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    _, cls, nc = compound_meta(ROOT / "data" / "unified" / "compounds.csv")
    provider = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    mp = Mapper(reg, cls, nc, provider)

    props = pd.read_csv(PROPS)
    kv = props[(props.property == "viscosity_kin_cSt") & props.fuel_id.str.startswith("CRC647_")]
    assert set(kv.provenance) == {"raw_file"}
    comp = pd.read_csv(COMP)
    direct = bool((comp.provenance == "reconstructed_r2").any())
    if direct:
        assert (comp.provenance == "reconstructed_r2").all(), "mixed r2 / non-r2 composition table"
    splits = json.load(open(SPLITS))
    locked = set(splits["locked"])
    fuels = sorted(kv.fuel_id.unique())
    dev = [f for f in fuels if f not in locked]
    lock = [f for f in fuels if f in locked]
    assert len(fuels) == 57 and len(lock) == 12 and len(dev) == 45, (len(fuels), len(dev), len(lock))
    print(f"fuels {len(fuels)}  dev {len(dev)}  locked {len(lock)}", flush=True)

    rows, rec = [], []
    drop = {}
    for f in fuels:
        idx, w, desc, dropped = fuel_composition(comp[comp.system_id == f"fuel:{f}"], mp, direct=direct)
        drop[f] = dropped
        for r in kv[kv.fuel_id == f].itertuples():
            T_K = r.T_C + 273.15
            nu, rho_mix, ln_eta, rho = refutas_nu(provider, idx, w, T_K)
            z = np.column_stack([reg.S[idx], np.nan_to_num(ln_eta[:, None]), np.nan_to_num(rho[:, None]) / 1000.0,
                                 np.full((len(idx), 1), 1000.0 / T_K)])
            rows.append(dict(z=torch.tensor(z, dtype=torch.float32), w=torch.tensor(w, dtype=torch.float32)))
            rec.append(dict(fuel_id=f, T_C=r.T_C, T_K=T_K, y=float(r.value), nu_refutas=nu, **desc))
    D = pd.DataFrame(rec)
    fid = D.fuel_id.to_numpy(); T_K = D.T_K.to_numpy(); T_C = D.T_C.to_numpy()
    y_ln = np.log(D.y.to_numpy()); phys_ln = np.log(D.nu_refutas.to_numpy())
    resid = y_ln - phys_ln
    X = D[FEATS].to_numpy(float)
    inv = 1000.0 / T_K
    is_dev = np.isin(fid, dev); is_lock = np.isin(fid, lock)
    dev_idx = np.where(is_dev)[0]; lock_idx = np.where(is_lock)[0]
    print(f"mapped compositions; mean dropped (unmappable) mass {np.mean(list(drop.values())):.4f}; "
          f"max {max(drop.values()):.4f}", flush=True)

    P = {}                                     # model -> ln-nu predictions for all rows (dev: LOO, locked: dev-trained)
    meta = {}
    hp_log = []

    # --- zero-parameter physics
    P["Raw Refutas"] = phys_ln.copy()
    meta["Raw Refutas"] = ("none (zero-parameter physics)", "zero-parameter", 0)

    # --- recalibrated Refutas (2 parameters: intercept + slope in 1000/T), LOO over dev, dev-fit for locked
    def recal_fit_predict(tr, te, hp=None):
        B = np.column_stack([np.ones(len(tr)), inv[tr]])
        beta = np.linalg.lstsq(B, resid[tr], rcond=None)[0]
        return phys_ln[te] + np.column_stack([np.ones(len(te)), inv[te]]) @ beta

    recal = np.full(len(D), np.nan)
    for f in dev:
        te = np.where(fid == f)[0]; tr = np.setdiff1d(dev_idx, te)
        recal[te] = recal_fit_predict(tr, te)
    recal[lock_idx] = recal_fit_predict(dev_idx, lock_idx)
    P["Recalibrated Refutas"] = recal
    meta["Recalibrated Refutas"] = ("dev fuels, LOO by fuel", "LOO by fuel (no tuning)", 2)

    # --- Walther (ASTM D341), anchored on Refutas at +20 and -20 C, extrapolated to -40 C (as v3)
    def walther(base_ln):
        out = np.full(len(D), np.nan)
        for f in fuels:
            ix = {tc: np.where((fid == f) & (T_C == tc))[0][0] for tc in TEMPS_C}
            A, B = fit_walther_parameters([293.15, 253.15], [np.exp(base_ln[ix[20.0]]), np.exp(base_ln[ix[-20.0]])])
            out[ix[20.0]], out[ix[-20.0]] = base_ln[ix[20.0]], base_ln[ix[-20.0]]
            out[ix[-40.0]] = np.log(walther_viscosity(233.15, A, B))
        return out

    P["Walther (D341)"] = walther(phys_ln)
    meta["Walther (D341)"] = ("none (anchored on raw Refutas 20/-20 C)", "zero-parameter", 0)
    P["Walther + recalibration"] = walther(recal)
    meta["Walther + recalibration"] = ("dev fuels, LOO by fuel (anchored on recalibrated Refutas)", "LOO by fuel (no tuning)", 2)

    # --- generic nested runner
    def run_nested(name, grid, fit_predict, label, n_params):
        t1 = time.time()
        pred_res = np.full(len(D), np.nan)
        # outer LOO over dev fuels; rows of locked fuels are excluded from selection AND training
        pr, chosen = nested_loo(fid, resid, grid, fit_predict, outer_fuels=dev, universe=dev_idx)
        pred_res[dev_idx] = pr[dev_idx]
        hp_final, _ = select_by_inner_cv(dev_idx, fid, grid, fit_predict, resid)     # final model: dev rows only
        pred_res[lock_idx] = fit_predict(dev_idx, lock_idx, hp_final)
        for f, h in chosen.items():
            hp_log.append(dict(model=name, fold=f, hp=str(h)))
        hp_log.append(dict(model=name, fold="FINAL(dev->locked)", hp=str(hp_final)))
        P[name] = phys_ln + pred_res
        meta[name] = (label, "nested LOO by fuel (inner grouped CV)", n_params)
        print(f"  {name}: {time.time() - t1:.0f}s  chosen {pd.Series([str(v) for v in chosen.values()]).value_counts().to_dict()}", flush=True)
        return pred_res

    # --- bilinear pairwise
    print("bilinear...", flush=True)
    bfp = bilinear_fit_predict_factory(rows, resid)
    run_nested("Bilinear pairwise", BIL_GRID, bfp, "dev fuels (nested)", "rank*(23+1)+23")

    # --- structured ridge
    def ridge_fp(tr, te, a):
        return StructuredResidual(FEATS, "ridge", alpha=a).fit(X[tr], T_K[tr], resid[tr]).predict(X[te], T_K[te])
    print("structured ridge...", flush=True)
    run_nested("Structured ridge", ALPHA_GRID, ridge_fp, "dev fuels (nested)", 2 + 2 * len(FEATS))

    # --- structured GP (+ 90% interval)
    def gp_fit(tr, hp):
        return StructuredResidual(FEATS, "gp", length_scale=hp[0], noise_ratio=hp[1]).fit(X[tr], T_K[tr], resid[tr])

    def gp_fp(tr, te, hp):
        return gp_fit(tr, hp).predict(X[te], T_K[te])
    print("structured GP...", flush=True)
    gp_res = run_nested("Structured GP", GP_GRID, gp_fp, "dev fuels (nested)", "nonparametric")
    # intervals: refit with the hyper-parameter chosen for that fold
    lo90 = np.full(len(D), np.nan); hi90 = np.full(len(D), np.nan)
    gp_hp = {r["fold"]: ast.literal_eval(r["hp"]) for r in hp_log if r["model"] == "Structured GP"}
    for f in dev:
        te = np.where(fid == f)[0]; tr = np.setdiff1d(dev_idx, te)
        mu, sd = gp_fit(tr, gp_hp[f]).predict(X[te], T_K[te], return_std=True)
        lo90[te], hi90[te] = phys_ln[te] + mu - 1.645 * sd, phys_ln[te] + mu + 1.645 * sd
    mu, sd = gp_fit(dev_idx, gp_hp["FINAL(dev->locked)"]).predict(X[lock_idx], T_K[lock_idx], return_std=True)
    lo90[lock_idx], hi90[lock_idx] = phys_ln[lock_idx] + mu - 1.645 * sd, phys_ln[lock_idx] + mu + 1.645 * sd

    # --- lab-trained structured residual applied zero-shot (never sees any CRC fuel)
    print("lab-trained structured residual...", flush=True)
    lab = build_lab(reg, provider, cls, nc)
    Xl, Tl, rl, sl = lab[FEATS].to_numpy(float), lab["T_K"].to_numpy(), lab["resid"].to_numpy(), lab.system_id.to_numpy()

    def lab_fp(tr, te, a):
        return StructuredResidual(FEATS, "ridge", alpha=a).fit(Xl[tr], Tl[tr], rl[tr]).predict(Xl[te], Tl[te])
    a_lab, lab_scores = select_by_inner_cv(np.arange(len(lab)), sl, ALPHA_GRID, lab_fp, rl)
    lab_model = StructuredResidual(FEATS, "ridge", alpha=a_lab).fit(Xl, Tl, rl)
    P["Lab-trained structured residual (zero-shot)"] = phys_ln + lab_model.predict(X, T_K)
    meta["Lab-trained structured residual (zero-shot)"] = (
        f"lab binary/multi mixtures only ({len(lab)} records, {lab.system_id.nunique()} systems, "
        f"{LAB_T_RANGE[0]:.0f}-{LAB_T_RANGE[1]:.0f} K); no CRC fuel", "zero-shot", 2 + 2 * len(FEATS))
    hp_log.append(dict(model="Lab-trained structured residual (zero-shot)", fold="lab grouped CV", hp=str(a_lab)))

    # ------------------------------------------------------------------ leaderboard
    err = {m: p - y_ln for m, p in P.items()}
    conds = [("overall (all T)", np.ones(len(D), bool)), ("T = 20 C", T_C == 20.0),
             ("T = -20 C", T_C == -20.0), ("T = -40 C", T_C == -40.0)]
    splits_ = [("dev (LOO by fuel)", is_dev), ("locked (untouched)", is_lock)]
    recs = []
    base = "Recalibrated Refutas"
    for m in P:
        for sname, smask in splits_:
            for cname, cmask in conds:
                mk = smask & cmask
                e, c = err[m][mk], fid[mk]
                r_pt, r_lo, r_hi = cluster_bootstrap_ci(e, c, rmse_log, n_resamples=2000, seed=SEED)
                a_pt, a_lo, a_hi = cluster_bootstrap_ci(e, c, mae_log, n_resamples=2000, seed=SEED)
                row = dict(model=m, split=sname, condition=cname, n_fuels=len(np.unique(c)), n_obs=int(mk.sum()),
                           training_data=meta[m][0], protocol=meta[m][1], n_params=meta[m][2],
                           log_rmse=round(r_pt, 4), log_rmse_ci95=f"[{r_lo:.4f}, {r_hi:.4f}]",
                           log_mae=round(a_pt, 4), log_mae_ci95=f"[{a_lo:.4f}, {a_hi:.4f}]",
                           mape_pct=round(float(np.mean(np.abs(np.expm1(e))) * 100), 2),
                           bias_log=round(float(np.mean(e)), 4))
                if m != base:
                    for nm, fn in (("rmse", rmse_log), ("mae", mae_log)):
                        d, lo, hi = paired_bootstrap(e, err[base][mk], c, fn)
                        row[f"d_log_{nm}_vs_recal"] = round(d, 4)
                        row[f"d_log_{nm}_vs_recal_ci95"] = f"[{lo:.4f}, {hi:.4f}]"
                        row[f"{nm}_ci_excludes_0"] = bool(lo > 0 or hi < 0)
                if m == "Structured GP":
                    inside = (y_ln[mk] >= lo90[mk]) & (y_ln[mk] <= hi90[mk])
                    row["gp_cov90"] = round(float(inside.mean()), 3)
                    row["gp_mean_width90_log"] = round(float(np.mean(hi90[mk] - lo90[mk])), 3)
                recs.append(row)
    lb = pd.DataFrame(recs)
    lb.to_csv(OUT / "viscosity_leaderboard_v4.csv", index=False)
    pr = D.copy()
    pr["locked"] = is_lock
    for m, p in P.items():
        pr[f"pred_{m}"] = np.exp(p)
    pr["gp90_lo"], pr["gp90_hi"] = np.exp(lo90), np.exp(hi90)
    pr["unmapped_mass"] = [drop[f] for f in fid]
    pr.to_csv(OUT / "viscosity_predictions_v4.csv", index=False)
    pd.DataFrame(hp_log).to_csv(OUT / "selected_hyperparameters_v4.csv", index=False)
    print(f"done in {time.time() - t0:.0f}s; lab records {len(lab)}; lab alpha {a_lab}", flush=True)


if __name__ == "__main__":
    main()
