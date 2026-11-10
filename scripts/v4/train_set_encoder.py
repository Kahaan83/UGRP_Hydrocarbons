"""v4 Task 16: train / evaluate the residual set encoder (CPU or GPU) against Task 11's best viscosity model.

    python scripts/v4/train_set_encoder.py --device cpu --seeds 5 --epochs 60 --out results/v4/E

Data
  * lab mixtures: data/unified/measurements_dedup.csv hydrocarbon binary/multi dynamic viscosity (same filter as
    Task 11's lab row; GroupKFold key = system_id) and the 126 Task 9 lab flash records (results/v4/F/flash_lab_loso.csv).
  * v4 real fuels: the 57 CRC 647 jet fuels (kinematic viscosity at +20/-20/-40 C, raw_file; D93 flash, raw_file) with
    their Task 3 reconstructed compositions; GroupKFold key = fuel_id.  The 12 locked fuels in
    splits/v4/real_fuels_v4_splits.json are NEVER in a training, early-stopping (none is used), selection or
    calibration index; they are predicted once by the dev-trained ensemble and scored as a separate count.
  * hyper-parameters are FIXED a priori (no selection => nothing can leak); see README.

Heads: viscosity (residual on ln nu vs raw Refutas, structured a0 + a1*1000/T) and flash (residual in degC on the Task 9
ideal Le Chatelier prediction).  Loss: masked multi-task Huber; every row carries one label.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.mixture.set_encoder_v4 import assert_no_locked, dev_group_folds, make_encoder  # noqa: E402
from fuelblend.mixture.structured_residual import StructuredResidual, group_kfold, select_by_inner_cv  # noqa: E402
from fuelblend.uq.conformal import conformal_rank, coverage, two_sided_margin, one_sided_margin  # noqa: E402

SEED = 20261110
PRED_T11 = ROOT / "results" / "v4" / "V" / "viscosity_predictions_v4.csv"
FLASH_LAB = ROOT / "results" / "v4" / "F" / "flash_lab_loso.csv"
FLASH_REAL = ROOT / "results" / "v4" / "F" / "flash_real_fuels.csv"
CACHE = ROOT / "runs" / "v4" / "cache"


def _load_lb():
    p = ROOT / "scripts" / "v4" / "build_viscosity_leaderboard_v4.py"
    spec = importlib.util.spec_from_file_location("v4_visc_lb", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules["v4_visc_lb"] = m
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------------------------------- data
class Task:
    """One (domain, head) table: padded components [R,N,d], weights [R,N], T_K [R], ctx [R,1], y_res [R], group [R]."""

    def __init__(self, head, domain, comp, w, T, ctx, y, group, extra=None):
        self.head, self.domain = head, domain
        self.comp, self.w, self.T, self.ctx, self.y = comp, w, T, ctx, y
        self.group = np.asarray(group).astype(str)
        self.extra = extra or {}

    def __len__(self):
        return len(self.y)


def _pad(items, S, d):
    n = max(len(i) for i, _ in items)
    comp = np.zeros((len(items), n, d), np.float32)
    w = np.zeros((len(items), n), np.float32)
    for r, (idx, ww) in enumerate(items):
        comp[r, :len(idx)] = S[idx]
        w[r, :len(idx)] = ww
    return comp, w


def prepare(lb):
    """Build the four task tables + side info.  Deterministic; cached to runs/v4/cache/prep.pkl."""
    from fuelblend.mixture.cold_flow import compound_meta
    from fuelblend.mixture.common import CLASSES, CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    _, cls, nc = compound_meta(ROOT / "data" / "unified" / "compounds.csv")
    provider = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    mp = lb.Mapper(reg, cls, nc, provider)
    S, d = reg.S, reg.S.shape[1]

    # ---- lab viscosity (same filter as Task 11 lab row)
    m = pd.read_csv(lb.MEAS, low_memory=False)
    wv = m[m.kind.isin(["binary", "multi"]) & m.is_hc.astype(bool) & (m.property == "viscosity_dyn_mPas")
           & ~m.conflict.astype(bool) & m["T_K"].between(*lb.LAB_T_RANGE)]
    items, T, y, grp, feats = [], [], [], [], []
    for r in wv.itertuples():
        ks = json.loads(r.comp_inchikeys); ws = np.array(json.loads(r.w_mass), float)
        ii = [reg.ik2i.get(k) for k in ks]
        if any(i is None or not provider.valid[i] or not reg.usable[i] or cls[i] not in CLASSES
               or not np.isfinite(nc[i]) for i in ii) or not np.all(np.isfinite(ws)) or ws.sum() <= 0:
            continue
        ii = np.array(ii); ws = ws / ws.sum()
        nu, rho_mix, _, _ = lb.refutas_nu(provider, ii, ws, r.T_K)
        items.append((ii, ws)); T.append(r.T_K); grp.append(r.system_id)
        y.append(np.log(r.value) - np.log(nu * rho_mix / 1000.0))
        dd = {f"frac_{c}": float(sum(wi for i, wi in zip(ii, ws) if cls[i] == c)) for c in CLASSES}
        n_ = np.array([nc[i] for i in ii]); m_ = float((ws * n_).sum())
        dd["nC_mean"], dd["nC_std"] = m_, float(np.sqrt((ws * (n_ - m_) ** 2).sum()))
        feats.append([dd[f] for f in lb.FEATS])
    comp, w = _pad(items, S, d)
    lab_v = Task("viscosity", "lab", comp, w, np.array(T, np.float32), np.zeros((len(y), 1), np.float32),
                 np.array(y, np.float32), grp, dict(X=np.array(feats)))

    # ---- lab flash (Task 9 ideal baseline)
    from rdkit import Chem
    fl = pd.read_csv(FLASH_LAB)
    ik2smi = {}
    items, y, grp, T = [], [], [], []
    for r in fl.itertuples():
        smis = json.loads(r.comp_smiles); xs = np.array(json.loads(r.x_mole), float)
        ii = []
        for s in smis:
            ik = Chem.MolToInchiKey(Chem.MolFromSmiles(s))
            ii.append(reg.ik2i.get(ik))
        if any(i is None or not reg.usable[i] for i in ii):
            continue
        ii = np.array(ii); mass = xs * reg.mw[ii]; mass = mass / mass.sum()
        items.append((ii, mass)); y.append(r.measured_C - r.ideal_C); grp.append(r.system_id); T.append(293.15)
    comp, w = _pad(items, S, d)
    lab_f = Task("flash", "lab", comp, w, np.array(T, np.float32), np.zeros((len(y), 1), np.float32),
                 np.array(y, np.float32), grp)

    # ---- real fuels (CRC 647)
    props = pd.read_csv(lb.PROPS)
    kv = props[(props.property == "viscosity_kin_cSt") & props.fuel_id.str.startswith("CRC647_")]
    comp_df = pd.read_csv(lb.COMP)
    fr = pd.read_csv(FLASH_REAL)
    fr = fr[(fr.method == "D93") & fr.fuel_id.str.startswith("CRC647_")]
    fuels = sorted(kv.fuel_id.unique())
    fcomp = {}
    for f in fuels:
        idx, ww, desc, dropped = lb.fuel_composition(comp_df[comp_df.system_id == f"fuel:{f}"], mp)
        fcomp[f] = (idx, ww, desc)
    items, T, y, grp, feats, meta = [], [], [], [], [], []
    for f in fuels:
        idx, ww, desc = fcomp[f]
        for r in kv[kv.fuel_id == f].itertuples():
            T_K = r.T_C + 273.15
            nu, _, _, _ = lb.refutas_nu(provider, idx, ww, T_K)
            items.append((idx, ww)); T.append(T_K); grp.append(f)
            y.append(np.log(r.value) - np.log(nu)); feats.append([desc[k] for k in lb.FEATS])
            meta.append((r.T_C, nu, float(r.value)))
    comp, w = _pad(items, S, d)
    real_v = Task("viscosity", "real", comp, w, np.array(T, np.float32), np.ones((len(y), 1), np.float32),
                  np.array(y, np.float32), grp,
                  dict(X=np.array(feats), T_C=np.array([m_[0] for m_ in meta]), nu_ref=np.array([m_[1] for m_ in meta]),
                       nu_meas=np.array([m_[2] for m_ in meta])))
    items, y, grp, T, ideal = [], [], [], [], []
    for r in fr.itertuples():
        if r.fuel_id not in fcomp or not np.isfinite(r.ideal_C):
            continue
        idx, ww, _ = fcomp[r.fuel_id]
        items.append((idx, ww)); y.append(r.measured_C - r.ideal_C); grp.append(r.fuel_id); T.append(293.15)
        ideal.append(r.ideal_C)
    comp, w = _pad(items, S, d)
    real_f = Task("flash", "real", comp, w, np.array(T, np.float32), np.ones((len(y), 1), np.float32),
                  np.array(y, np.float32), grp, dict(ideal=np.array(ideal)))
    return dict(lab_v=lab_v, lab_f=lab_f, real_v=real_v, real_f=real_f, d_comp=d)


def get_data(lb, refresh=False):
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / "prep.pkl"
    if p.exists() and not refresh:
        return pickle.load(open(p, "rb"))
    D = prepare(lb)
    pickle.dump(D, open(p, "wb"))
    return D


# ------------------------------------------------------------------------------------------- training
def sub(task: Task, idx):
    idx = np.asarray(idx, int)
    return Task(task.head, task.domain, task.comp[idx], task.w[idx], task.T[idx], task.ctx[idx], task.y[idx],
                task.group[idx], {k: v[idx] for k, v in task.extra.items()})


def _tt(task, dev):
    return {k: torch.as_tensor(getattr(task, k)).to(dev) for k in ("comp", "w", "T", "ctx", "y")}


def train_one(tasks: list, d_comp: int, pool: str, seed: int, epochs: int, device, lab_bs=512, real_bs=32,
              lr=3e-3, wd=1e-3, hidden=64):
    """Masked multi-task training: one minibatch from every non-empty task table per step, summed Huber losses.
    Fixed schedule (no early stopping / no selection) -> no held-out index is ever consulted."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = make_encoder(d_comp, pool, hidden, d_context=1).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1))
    tt = [(t, _tt(t, device)) for t in tasks if len(t) > 0]
    n_big = max(len(t) for t, _ in tt if t.domain == "lab" and t.head == "viscosity") if any(
        t.domain == "lab" and t.head == "viscosity" for t, _ in tt) else max(len(t) for t, _ in tt)
    steps = max(1, int(np.ceil(n_big / lab_bs)))
    scale = {"viscosity": 0.3, "flash": 5.0}
    for ep in range(epochs):
        model.train()
        for _ in range(steps):
            loss = 0.0
            for t, g in tt:
                bs = lab_bs if t.domain == "lab" and t.head == "viscosity" else real_bs
                b = torch.as_tensor(rng.integers(0, len(t), min(bs, len(t))), device=device)
                pred = model.residual(t.head, g["comp"][b], g["w"][b], g["T"][b], g["ctx"][b])
                loss = loss + torch.nn.functional.huber_loss(pred / scale[t.head], g["y"][b] / scale[t.head], delta=1.0)
            opt.zero_grad(); loss.backward(); opt.step()
        sched.step()
    return model.eval()


@torch.no_grad()
def predict_residual(models, task: Task, device):
    """Per-seed residual predictions [n_seeds, n_rows] for task rows."""
    g = _tt(task, device)
    return np.stack([m.residual(task.head, g["comp"], g["w"], g["T"], g["ctx"]).cpu().numpy() for m in models])


def train_ensemble(tasks, d_comp, pool, seeds, epochs, device, **kw):
    return [train_one(tasks, d_comp, pool, s, epochs, device, **kw) for s in seeds]


# ------------------------------------------------------------------------------------------- metrics
def rmse(e): return float(np.sqrt(np.mean(np.square(e))))
def mae(e): return float(np.mean(np.abs(e)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--out", default=str(ROOT / "results" / "v4" / "E"))
    ap.add_argument("--folds", type=int, default=9, help="outer grouped folds over the dev real fuels")
    ap.add_argument("--lab-folds", type=int, default=5, help="GroupKFold(system_id) folds for the lab evaluation (0=skip)")
    ap.add_argument("--pools", default="attention,sum")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--save-models", default="", help="directory for .pt checkpoints of the final dev-trained ensemble")
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args(argv)
    if a.threads:
        torch.set_num_threads(a.threads)
    dev_ = torch.device(a.device)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    lb = _load_lb()
    D = get_data(lb, a.refresh)
    d = D["d_comp"]
    rv, rf, lv, lf = D["real_v"], D["real_f"], D["lab_v"], D["lab_f"]
    locked = set(json.load(open(lb.SPLITS))["locked"])
    seeds = [SEED + i for i in range(a.seeds)]
    print(f"data: lab visc {len(lv)} ({len(set(lv.group))} sys) lab flash {len(lf)} real visc {len(rv)} "
          f"({len(set(rv.group))} fuels) real flash {len(rf)}; threads {torch.get_num_threads()}; prep {time.time()-t0:.0f}s", flush=True)

    # sanity: our Refutas baseline reproduces Task 11's stored predictions
    t11 = pd.read_csv(PRED_T11)
    key = {(f, round(tc, 3)): i for i, (f, tc) in enumerate(zip(t11.fuel_id, t11.T_C))}
    pos = np.array([key[(f, round(tc, 3))] for f, tc in zip(rv.group, rv.extra["T_C"])])
    assert np.allclose(t11.nu_refutas.to_numpy()[pos], rv.extra["nu_ref"], rtol=1e-6), "baseline differs from Task 11"

    fid = rv.group
    is_lock = np.isin(fid, list(locked))
    dev_rows = np.where(~is_lock)[0]; lock_rows = np.where(is_lock)[0]
    fid_f = rf.group
    flock_rows = np.where(np.isin(fid_f, list(locked)))[0]
    audit = dict(n_locked_fuels=len(locked), n_real_visc_rows_locked=int(is_lock.sum()))

    folds = dev_group_folds(fid, locked, a.folds, seed=SEED)
    assert sorted(np.concatenate([te for _, te in folds]).tolist()) == dev_rows.tolist()
    n_train_checks = 0

    def real_flash_for(groups):
        return np.where(np.isin(fid_f, list(groups)))[0]

    results = {}
    for pool in a.pools.split(","):
        pool = pool.strip()
        oof_v = np.full((a.seeds, len(rv)), np.nan); oof_f = np.full((a.seeds, len(rf)), np.nan)
        fold_of_row = np.full(len(rv), -1)
        for k, (tr, te) in enumerate(folds):
            test_f = set(fid[te])
            trf = np.where(~np.isin(fid_f, list(test_f | locked)))[0]
            tef = np.where(np.isin(fid_f, list(test_f)))[0]
            assert_no_locked(dict(train_visc=tr, test_visc=te), fid, locked)
            assert_no_locked(dict(train_flash=trf, test_flash=tef), fid_f, locked)
            assert not (set(fid[tr]) & set(fid[te])), "fold leak"
            n_train_checks += 1
            tasks = [lv, lf, sub(rv, tr), sub(rf, trf)]
            ms = train_ensemble(tasks, d, pool, seeds, a.epochs, dev_)
            oof_v[:, te] = predict_residual(ms, sub(rv, te), dev_)
            if len(tef):
                oof_f[:, tef] = predict_residual(ms, sub(rf, tef), dev_)
            fold_of_row[te] = k
            print(f"  [{pool}] fold {k+1}/{len(folds)} done ({time.time()-t0:.0f}s)", flush=True)
        # final dev -> locked
        tr_all = dev_rows; trf_all = np.where(~np.isin(fid_f, list(locked)))[0]
        assert_no_locked(dict(train_visc=tr_all), fid, locked)
        assert_no_locked(dict(train_flash=trf_all), fid_f, locked)
        ms = train_ensemble([lv, lf, sub(rv, tr_all), sub(rf, trf_all)], d, pool, seeds, a.epochs, dev_)
        lock_v = predict_residual(ms, sub(rv, lock_rows), dev_)
        lock_f = predict_residual(ms, sub(rf, flock_rows), dev_) if len(flock_rows) else np.zeros((a.seeds, 0))
        if a.save_models:
            Path(a.save_models).mkdir(parents=True, exist_ok=True)
            for s, m_ in zip(seeds, ms):
                torch.save(m_.state_dict(), Path(a.save_models) / f"set_encoder_{pool}_seed{s}.pt")
        # lab evaluation
        lab_res = None
        if a.lab_folds and pool == "attention":
            sysids = np.array(sorted(set(lv.group) | set(lf.group)))
            perm = np.random.default_rng(SEED).permutation(len(sysids))
            parts = np.array_split(sysids[perm], a.lab_folds)
            pv = np.full((a.seeds, len(lv)), np.nan); pf = np.full((a.seeds, len(lf)), np.nan)
            for k, p in enumerate(parts):
                trv = np.where(~np.isin(lv.group, p))[0]; tev = np.where(np.isin(lv.group, p))[0]
                trf_ = np.where(~np.isin(lf.group, p))[0]; tef_ = np.where(np.isin(lf.group, p))[0]
                assert not (set(lv.group[trv]) & set(lv.group[tev]))
                ms_l = train_ensemble([sub(lv, trv), sub(lf, trf_), sub(rv, dev_rows), sub(rf, trf_all)], d, pool,
                                      seeds, a.epochs, dev_)
                pv[:, tev] = predict_residual(ms_l, sub(lv, tev), dev_)
                if len(tef_):
                    pf[:, tef_] = predict_residual(ms_l, sub(lf, tef_), dev_)
                print(f"  [lab] fold {k+1}/{len(parts)} ({time.time()-t0:.0f}s)", flush=True)
            lab_res = (pv, pf)
        results[pool] = dict(oof_v=oof_v, oof_f=oof_f, fold_of_row=fold_of_row, lock_v=lock_v, lock_f=lock_f, lab=lab_res)

    # ------------------------------------------------------------------ reference models on the SAME rows
    y = rv.y                                     # ln(meas) - ln(Refutas)
    Xr, Tr = rv.extra["X"], rv.T.astype(float)
    ridge_oof = np.full(len(rv), np.nan)
    for k, (tr, te) in enumerate(folds):
        def fp(a_tr, a_te, al):
            return StructuredResidual(lb.FEATS, "ridge", alpha=al).fit(Xr[a_tr], Tr[a_tr], y[a_tr]).predict(Xr[a_te], Tr[a_te])
        al, _ = select_by_inner_cv(tr, fid, lb.ALPHA_GRID, fp, y)
        ridge_oof[te] = fp(tr, te, al)
    al_f, _ = select_by_inner_cv(dev_rows, fid, lb.ALPHA_GRID, fp, y)
    ridge_lock = fp(dev_rows, lock_rows, al_f)
    t11_ridge = np.log(t11["pred_Structured ridge"].to_numpy()[pos]) - np.log(t11.nu_refutas.to_numpy()[pos])
    t11_recal = np.log(t11["pred_Recalibrated Refutas"].to_numpy()[pos]) - np.log(t11.nu_refutas.to_numpy()[pos])

    # error vectors (pred_ln - meas_ln = pred_resid - y)
    E = {"Raw Refutas": np.zeros(len(rv)) - y,
         "Recalibrated Refutas (Task 11 LOO)": t11_recal - y,
         "Structured ridge (Task 11 LOO, published best)": t11_ridge - y,
         f"Structured ridge (same {a.folds}-fold, inner-CV alpha)": np.where(is_lock, np.nan, ridge_oof) - y}
    E["Structured ridge (same fold)"] = E.pop(f"Structured ridge (same {a.folds}-fold, inner-CV alpha)")
    lock_extra = {"Structured ridge (same fold)": ridge_lock}
    for pool, R in results.items():
        e = np.full(len(rv), np.nan)
        e[dev_rows] = R["oof_v"].mean(0)[dev_rows] - y[dev_rows]
        e[lock_rows] = R["lock_v"].mean(0) - y[lock_rows]
        E[f"Set encoder ({pool}, {a.seeds}-seed ensemble)"] = e
    for nm, v in lock_extra.items():
        E[nm][lock_rows] = v - y[lock_rows]
    T_C = rv.extra["T_C"]
    conds = [("overall (all T)", np.ones(len(rv), bool)), ("T = 20 C", T_C == 20.0), ("T = -20 C", T_C == -20.0),
             ("T = -40 C", T_C == -40.0)]
    rows = []
    primary = f"Set encoder (attention, {a.seeds}-seed ensemble)"
    refs = ["Structured ridge (Task 11 LOO, published best)", "Structured ridge (same fold)"]
    for m, e in E.items():
        for sname, smask in (("dev (outer CV by fuel)", ~is_lock), ("locked (untouched)", is_lock)):
            for cname, cm in conds:
                mk = smask & cm
                if not np.isfinite(e[mk]).all():
                    continue
                pt, lo, hi = lb.cluster_bootstrap_ci(e[mk], fid[mk], lb.rmse_log, n_resamples=2000, seed=SEED)
                row = dict(model=m, split=sname, condition=cname, n_fuels=len(set(fid[mk])), n_obs=int(mk.sum()),
                           log_rmse=round(pt, 4), log_rmse_ci95=f"[{lo:.4f}, {hi:.4f}]", log_mae=round(mae(e[mk]), 4),
                           mape_pct=round(float(np.mean(np.abs(np.expm1(e[mk]))) * 100), 2), bias_log=round(float(e[mk].mean()), 4))
                if m.startswith("Set encoder"):
                    for rname, short in zip(refs, ("t11loo", "samefold")):
                        if np.isfinite(E[rname][mk]).all():
                            dd, dlo, dhi = lb.paired_bootstrap(e[mk], E[rname][mk], fid[mk], lb.rmse_log)
                            row[f"d_rmse_vs_{short}"] = round(dd, 4)
                            row[f"d_rmse_vs_{short}_ci95"] = f"[{dlo:.4f}, {dhi:.4f}]"
                            row[f"ci_excludes_0_vs_{short}"] = bool(dlo > 0 or dhi < 0)
                rows.append(row)
    pd.DataFrame(rows).to_csv(out / "viscosity_encoder_vs_task11.csv", index=False)

    # ------------------------------------------------------------------ conformal (Task 10 functions, dev groups only)
    cal = []
    for pool, R in results.items():
        e = np.full(len(rv), np.nan)
        e[dev_rows] = R["oof_v"].mean(0)[dev_rows] - y[dev_rows]
        el = R["lock_v"].mean(0) - y[lock_rows]
        assert_no_locked(dict(calibration=dev_rows), fid, locked)
        m_all = two_sided_margin(e[dev_rows], fid[dev_rows], 0.10)
        # cross-conformal on dev: margin for fold f from the OTHER folds' OOF residuals, evaluated on fold f
        cov_rows = cov_grp = n_rows = n_grp = 0
        for k, (tr, te) in enumerate(folds):
            o = np.setdiff1d(dev_rows, te)
            mk_ = two_sided_margin(e[o], fid[o], 0.10)
            c, n = coverage(e[te], fid[te], mk_, two_sided=True)
            cov_grp += c; n_grp += n
            cov_rows += int(np.sum(np.abs(e[te]) <= mk_)); n_rows += len(te)
        cl, nl = coverage(el, fid[lock_rows], m_all, two_sided=True)
        cal.append(dict(model=f"Set encoder ({pool})", head="viscosity", alpha=0.10, kind="two-sided |ln residual|, group=fuel",
                        n_dev_groups=len(set(fid[dev_rows])), conformal_rank_k=conformal_rank(len(set(fid[dev_rows])), 0.10),
                        margin_ln=m_all, margin_multiplicative=float(np.exp(m_all)) if np.isfinite(m_all) else np.inf,
                        dev_cross_conformal_row_cov=cov_rows / n_rows, dev_cross_conformal_rows=f"{cov_rows}/{n_rows}",
                        dev_cross_conformal_group_cov=cov_grp / n_grp, dev_cross_conformal_groups=f"{cov_grp}/{n_grp}",
                        locked_groups_covered=f"{cl}/{nl}", locked_rows_covered=f"{int(np.sum(np.abs(el) <= m_all))}/{len(el)}",
                        ensemble_sd_mean_dev=float(results[pool]["oof_v"].std(0, ddof=0)[dev_rows].mean())))
    # flash (one-sided over-prediction, D93 fuels)
    for pool, R in results.items():
        dev_f = np.where(~np.isin(fid_f, list(locked)))[0]
        resid_dev = (R["oof_f"].mean(0) - rf.y)[dev_f]            # pred - measured = pred_res - y_res
        q = one_sided_margin(resid_dev, fid_f[dev_f], 0.10, "over")
        applied = max(q, 0.0) if np.isfinite(q) else np.inf
        rl = R["lock_f"].mean(0) - rf.y[flock_rows]
        cl, nl = coverage(rl, fid_f[flock_rows], applied, direction="over")
        cal.append(dict(model=f"Set encoder ({pool})", head="flash (D93 CRC fuels)", alpha=0.10, kind="one-sided over-prediction, group=fuel",
                        n_dev_groups=len(set(fid_f[dev_f])), conformal_rank_k=conformal_rank(len(set(fid_f[dev_f])), 0.10),
                        margin_ln=q, applied_margin_C=applied, locked_groups_covered=f"{cl}/{nl}"))
    pd.DataFrame(cal).to_csv(out / "conformal_calibration.csv", index=False)

    # ------------------------------------------------------------------ flash + lab tables
    frows = []
    dev_f = np.where(~np.isin(fid_f, list(locked)))[0]
    for sname, ix in (("dev (outer CV by fuel)", dev_f), ("locked (untouched)", flock_rows)):
        base = -rf.y[ix]                                      # ideal pred - measured (= -(measured - ideal))
        frows.append(dict(model="Task 9 ideal (zero parameter)", split=sname, n=len(ix), mae_C=mae(base), bias_C=float(base.mean()),
                          max_over_C=float(base.max()), n_over_5C=int((base > 5).sum())))
        if sname.startswith("dev"):
            offs = []
            for tr, te in folds:
                trf = np.where(~np.isin(fid_f, list(set(fid[te]) | locked)))[0]
                off = rf.y[trf].mean()
                tef = np.where(np.isin(fid_f, list(set(fid[te]))))[0]
                offs.append((tef, off))
            c = np.full(len(rf), np.nan)
            for tef, off in offs:
                c[tef] = off - rf.y[tef]
            e = c[ix]
        else:
            off = rf.y[dev_f].mean(); e = off - rf.y[ix]
        frows.append(dict(model="ideal + constant offset (fold-trained mean)", split=sname, n=len(ix), mae_C=mae(e), bias_C=float(e.mean()),
                          max_over_C=float(e.max()), n_over_5C=int((e > 5).sum())))
        for pool, R in results.items():
            e = (R["oof_f"].mean(0) - rf.y)[ix] if sname.startswith("dev") else (R["lock_f"].mean(0) - rf.y[ix])
            frows.append(dict(model=f"Set encoder ({pool})", split=sname, n=len(ix), mae_C=mae(e), bias_C=float(e.mean()),
                              max_over_C=float(e.max()), n_over_5C=int((e > 5).sum())))
    pd.DataFrame(frows).round(4).to_csv(out / "flash_encoder_real_d93.csv", index=False)

    lrows = []
    R = results.get("attention")
    if R is not None and R["lab"] is not None:
        pv, pf = R["lab"]
        ev = pv.mean(0) - lv.y
        # lab-trained structured ridge, same system folds
        sysids = np.array(sorted(set(lv.group) | set(lf.group)))
        perm = np.random.default_rng(SEED).permutation(len(sysids))
        parts = np.array_split(sysids[perm], a.lab_folds)
        er = np.full(len(lv), np.nan)
        for p in parts:
            tr_ = np.where(~np.isin(lv.group, p))[0]; te_ = np.where(np.isin(lv.group, p))[0]
            Xl, Tl, yl, gl = lv.extra["X"], lv.T.astype(float), lv.y, lv.group
            def fpl(a_tr, a_te, al):
                return StructuredResidual(lb.FEATS, "ridge", alpha=al).fit(Xl[a_tr], Tl[a_tr], yl[a_tr]).predict(Xl[a_te], Tl[a_te])
            al, _ = select_by_inner_cv(tr_, gl, lb.ALPHA_GRID, fpl, yl)
            er[te_] = fpl(tr_, te_, al) - yl[te_]
        for nm, e in (("Raw Refutas", -lv.y), ("Lab structured ridge (system-grouped CV)", er), ("Set encoder (attention)", ev)):
            pt, lo, hi = lb.cluster_bootstrap_ci(e, lv.group, lb.rmse_log, n_resamples=1000, seed=SEED)
            lrows.append(dict(property="viscosity (ln dyn)", model=nm, n=len(e), n_systems=len(set(lv.group)),
                              rmse=round(pt, 4), rmse_ci95=f"[{lo:.4f}, {hi:.4f}]", mae=round(mae(e), 4), bias=round(float(e.mean()), 4)))
        ef = pf.mean(0) - lf.y
        for nm, e in (("Task 9 ideal", -lf.y), ("Set encoder (attention)", ef)):
            lrows.append(dict(property="flash (degC)", model=nm, n=len(e), n_systems=len(set(lf.group)), rmse=round(rmse(e), 4),
                              rmse_ci95="", mae=round(mae(e), 4), bias=round(float(e.mean()), 4)))
        pd.DataFrame(lrows).to_csv(out / "lab_groupkfold_system.csv", index=False)

    # ------------------------------------------------------------------ predictions + audit + decision
    pr = pd.DataFrame(dict(fuel_id=fid, T_C=T_C, locked=is_lock, ln_resid_target=y, nu_meas=rv.extra["nu_meas"], nu_refutas=rv.extra["nu_ref"],
                           ridge_same_fold_resid=np.where(is_lock, np.nan, ridge_oof)))
    for pool, R in results.items():
        r = np.full(len(rv), np.nan)
        r[dev_rows] = R["oof_v"].mean(0)[dev_rows]; r[lock_rows] = R["lock_v"].mean(0)
        pr[f"encoder_{pool}_resid"] = r
        sd = np.full(len(rv), np.nan); sd[dev_rows] = R["oof_v"].std(0)[dev_rows]; sd[lock_rows] = R["lock_v"].std(0)
        pr[f"encoder_{pool}_seed_sd"] = sd
    pr.to_csv(out / "viscosity_encoder_predictions.csv", index=False)
    audit.update(n_outer_folds=len(folds), n_fold_leak_checks=n_train_checks,
                 locked_in_any_train_or_calibration_index=0,
                 seeds=a.seeds, epochs=a.epochs, threads=int(torch.get_num_threads()), device=a.device,
                 n_lab_visc=len(lv), n_lab_flash=len(lf), n_real_visc=len(rv), n_real_flash=len(rf),
                 elapsed_s=round(time.time() - t0, 1), torch=torch.__version__)
    json.dump(audit, open(out / "split_audit.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))

    df = pd.DataFrame(rows)
    pm = f"Set encoder (attention, {a.seeds}-seed ensemble)"
    g = df[(df.model == pm) & (df.split.str.startswith("dev")) & (df.condition == "overall (all T)")]
    ca = [c for c in cal if c["model"] == "Set encoder (attention)" and c["head"] == "viscosity"]
    dec = {}
    if len(g) and ca:
        g = g.iloc[0]
        wins = bool(g.get("ci_excludes_0_vs_t11loo")) and g.get("d_rmse_vs_t11loo", 1) < 0 and \
            bool(g.get("ci_excludes_0_vs_samefold")) and g.get("d_rmse_vs_samefold", 1) < 0
        covok = ca[0]["dev_cross_conformal_row_cov"] >= 0.85
        dec = dict(beats_task11_both_protocols_ci_excl0=wins, dev_cross_conformal_row_coverage=ca[0]["dev_cross_conformal_row_cov"],
                   coverage_ge_085=covok, PROMOTE_TO_MODEL_DECISION=bool(wins and covok))
    json.dump(dec, open(out / "promotion_decision.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("decision:", dec, f"\ntotal {time.time()-t0:.0f}s", flush=True)
    return dec


if __name__ == "__main__":
    main()
