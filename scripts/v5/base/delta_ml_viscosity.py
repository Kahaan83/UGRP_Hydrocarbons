"""v5 delta_ml_viscosity: learned pair-interaction residual on Refutas mixture viscosity.
Declared text: results/v5/BASE/delta_ml_viscosity/README.md.

Usage:  python scripts/v5/base/delta_ml_viscosity.py --dev-only   (dev out-of-fold + model fit; locked systems are NOT scored)
        python scripts/v5/base/delta_ml_viscosity.py              (full run; locked systems scored once by the final model)
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "v5" / "BASE" / "delta_ml_viscosity"
MODEL_DIR = ROOT / "models" / "v5" / "ml_viscosity_pair"
RES_CSV = ROOT / "results" / "v5" / "S1" / "viscosity_lab_residuals_v5.csv"
SPLIT = ROOT / "results" / "v5" / "S1" / "lab_systems_splits_v5.json"
GRID_SHA = "b931e969926698dff52d69ff51dd99a990c24ff12eca545aa98e970afd3e7daf"
SEED, N_BOOT, WINS, XMIN = 42, 2000, 0.5, 0.01
T_WARM = 283.15
MODELS = ("ml", "const", "classpair")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_by_path(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------------------ data
def regenerate_rows(ctx) -> pd.DataFrame:
    """Same rows as S1 score_rows (refutas), keeping compositions.  Asserts positional agreement with the committed table."""
    from fuelblend.mixture.common import CLASSES
    s1 = load_by_path("s1_lab", "scripts/v5/eval_lab_known_composition_v5.py")
    lb = load_by_path("build_viscosity_leaderboard_v4", "scripts/v4/build_viscosity_leaderboard_v4.py")
    meas4 = pd.read_csv(s1.MEAS, low_memory=False)
    boehm = pd.read_csv(s1.BOEHM, low_memory=False)
    boehm["sources"] = boehm["sources"].fillna("boehm2022")
    meas = pd.concat([meas4, boehm], ignore_index=True)
    reg, prov, cls, nc = ctx.reg, ctx.prov, ctx.cls, ctx.nc
    rows = []
    for r in s1.lab_selection(meas).itertuples():
        ks = json.loads(r.comp_inchikeys)
        ws = np.array(json.loads(r.w_mass), float)
        xs = np.array(json.loads(r.x_mole), float)
        ii = [reg.ik2i.get(k) for k in ks]
        if any(i is None or not prov.valid[i] or not reg.usable[i] or cls[i] not in CLASSES or not np.isfinite(nc[i]) for i in ii) \
                or not np.all(np.isfinite(ws)) or ws.sum() <= 0 or not np.all(np.isfinite(xs)):
            continue
        ii = np.array(ii)
        ws = ws / ws.sum()
        T = float(r.T_K)
        nu, rho_mix, _, _ = lb.refutas_nu(prov, ii, ws, T)
        p = nu * rho_mix / 1000.0
        if np.isfinite(p) and p > 0:
            rows.append(dict(system_id=r.system_id, source=getattr(r, "sources", ""), T_K=T, n_comp=len(ks), ii=tuple(int(i) for i in ii),
                             x=tuple(xs), ln_pred=math.log(p), ln_meas=math.log(float(r.value))))
    R = pd.DataFrame(rows)
    com = pd.read_csv(RES_CSV)
    com = com[com.method == "refutas"].reset_index(drop=True)
    assert len(com) == len(R), f"row count {len(com)} vs regenerated {len(R)}"
    assert (com.system_id.to_numpy() == R.system_id.to_numpy()).all(), "system order differs"
    assert np.abs(com.T_K.to_numpy() - R.T_K.to_numpy()).max() < 1e-5
    assert np.abs(np.log(com.predicted.to_numpy()) - R.ln_pred.to_numpy()).max() < 1e-6, "predicted differs"
    assert np.abs(np.log(com.measured.to_numpy()) - R.ln_meas.to_numpy()).max() < 1e-6, "measured differs"
    sp = json.loads(SPLIT.read_text())
    role = {**{s: "dev" for s in sp["dev"]}, **{s: "locked" for s in sp["locked"]}}
    R["role"] = R.system_id.map(role)
    assert R.role.notna().all()
    assert (com.split_role.to_numpy() == R.role.to_numpy()).all()
    R["t_bin"] = com.t_bin.to_numpy()
    return R


def build_pairs(R: pd.DataFrame, ctx) -> pd.DataFrame:
    """One entry per (row, i<j) pair with xx = x_i x_j and the feature frame."""
    from fuelblend.ml.viscosity_pair import pair_feature_frame
    rid, ia, ib, xx = [], [], [], []
    for k, (ii, x) in enumerate(zip(R.ii, R.x)):
        for a in range(len(ii)):
            for b in range(a + 1, len(ii)):
                rid.append(k)
                ia.append(ii[a])
                ib.append(ii[b])
                xx.append(x[a] * x[b])
    P = pair_feature_frame(np.array(ia), np.array(ib), R.T_K.to_numpy()[np.array(rid)], ctx)
    P["row"] = np.array(rid)
    P["xx"] = np.array(xx)
    P["ia"], P["ib"] = np.array(ia), np.array(ib)
    return P


# ------------------------------------------------------------------------------------------------ models
def train_set(R, P, systems):
    sysarr = R.system_id.to_numpy()[P.row.to_numpy()]
    m = np.isin(sysarr, list(systems)) & (R.n_comp.to_numpy()[P.row.to_numpy()] == 2) & (P.xx.to_numpy() >= XMIN)
    return P[m]


def targets(R, Ptr):
    r = (R.ln_meas.to_numpy() - R.ln_pred.to_numpy())[Ptr.row.to_numpy()]
    g = np.clip(r, -WINS, WINS) / Ptr.xx.to_numpy()
    w = Ptr.xx.to_numpy() ** 2
    sys_ = R.system_id.to_numpy()[Ptr.row.to_numpy()]
    tot = pd.Series(w).groupby(sys_).transform("sum").to_numpy()
    return g, w / tot


def fit_lgb(Ptr, g, w):
    import lightgbm as lgb
    from fuelblend.ml.viscosity_pair import FEATURES, LGB_PARAMS
    m = lgb.LGBMRegressor(**LGB_PARAMS)
    m.fit(Ptr[FEATURES], g, sample_weight=w, categorical_feature=["class_pair"])
    return m


def fit_baselines(Ptr, g, w):
    g0 = float(np.average(g, weights=w))
    cp = pd.DataFrame({"c": Ptr.class_pair.to_numpy(), "g": g, "w": w})
    tab = {int(c): float(np.average(d.g, weights=d.w)) for c, d in cp.groupby("c")}
    return g0, tab


def predict_all(model, base, P):
    from fuelblend.ml.viscosity_pair import FEATURES
    g0, tab = base
    return dict(ml=model.booster_.predict(P[FEATURES]), const=np.full(len(P), g0),
                classpair=np.array([tab.get(int(c), g0) for c in P.class_pair]))


def corr_by_row(Pp, gv, n):
    return np.bincount(Pp.row.to_numpy(), weights=Pp.xx.to_numpy() * gv, minlength=n)


# ------------------------------------------------------------------------------------------------ metrics
def macro(R, err, mask) -> pd.Series:
    d = pd.DataFrame({"s": R.system_id.to_numpy()[mask], "e": np.abs(err[mask])})
    return d.groupby("s").e.mean()


def paired_boot(a: pd.Series, b: pd.Series, rng) -> dict:
    """a, b per-system MAE; gain = mean(a - b); 95 % percentile interval over system resamples."""
    d = (a - b.reindex(a.index)).to_numpy()
    n = len(d)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(N_BOOT)])
    return dict(n_systems=int(n), mae_a=float(a.mean()), mae_b=float(b.mean()), gain=float(d.mean()),
                lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)))


def conformal(scores, alpha=0.10):
    s = np.sort(np.asarray(scores, float))
    n = len(s)
    k = math.ceil((n + 1) * (1 - alpha))
    return (float(s[k - 1]) if k <= n else float("inf")), n, k


def gate(R, E, mask, a="refutas", b="ml"):
    return paired_boot(macro(R, E[a], mask), macro(R, E[b], mask), np.random.default_rng(SEED))


def pairs_of(ii):
    return [tuple(sorted((ii[a], ii[b]))) for a in range(len(ii)) for b in range(a + 1, len(ii))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dev-only", action="store_true")
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    assert sha(ROOT / "data" / "cache" / "pure_grid.npz") == GRID_SHA, "pure grid SHA-256 changed"
    from fuelblend.ml import viscosity_pair as vp
    ctx = vp.context()
    R = regenerate_rows(ctx)
    print(f"regenerated {len(R)} refutas rows, {R.system_id.nunique()} systems; positional agreement with committed table OK", flush=True)
    P = build_pairs(R, ctx)
    n = len(R)
    resid_ref = R.ln_pred.to_numpy() - R.ln_meas.to_numpy()        # pred - meas, S1 convention
    sysarr = R.system_id.to_numpy()
    dev_sys = sorted(R.system_id[R.role == "dev"].unique())
    lock_sys = sorted(R.system_id[R.role == "locked"].unique())

    # ---- dev GroupKFold(5), out-of-fold
    from sklearn.model_selection import GroupKFold
    dev_rows = np.where(R.role.to_numpy() == "dev")[0]
    oof = {k: np.zeros(n) for k in MODELS}
    fold_pairs = {}
    for f, (tr, te) in enumerate(GroupKFold(n_splits=5).split(dev_rows, groups=sysarr[dev_rows])):
        tr_sys = set(sysarr[dev_rows[tr]])
        Ptr = train_set(R, P, tr_sys)
        g, w = targets(R, Ptr)
        model, base = fit_lgb(Ptr, g, w), fit_baselines(Ptr, g, w)
        te_rows = dev_rows[te]
        Pte = P[P.row.isin(te_rows)]
        for k, gv in predict_all(model, base, Pte).items():
            oof[k][te_rows] = corr_by_row(Pte, gv, n)[te_rows]
        pairs = {tuple(sorted((a, b))) for a, b in zip(Ptr.ia, Ptr.ib)}
        for r in te_rows:
            fold_pairs[r] = pairs
        print(f"fold {f}: train systems {len(tr_sys)}, train pair rows {len(Ptr)}, test systems {len(set(sysarr[te_rows]))}", flush=True)

    # ---- final model on all dev systems
    Ptr = train_set(R, P, set(dev_sys))
    g, w = targets(R, Ptr)
    model, base = fit_lgb(Ptr, g, w), fit_baselines(Ptr, g, w)
    model.booster_.save_model(str(MODEL_DIR / "model.txt"))
    final_pairs = {tuple(sorted((a, b))) for a, b in zip(Ptr.ia, Ptr.ib)}
    cps = pd.DataFrame({"s": sysarr[Ptr.row.to_numpy()], "c": Ptr.class_pair.to_numpy()}).drop_duplicates()
    cp_counts = {vp.CLASS_PAIRS[int(c)]: int(k) for c, k in cps.groupby("c").size().items()}
    n_train_sys = int(len(set(sysarr[Ptr.row.to_numpy()])))
    print("training systems per class pair:", cp_counts, flush=True)

    # resid_ML = (ln_pred + corr) - ln_meas = resid_ref + corr, since corr approximates r = ln_meas - ln_pred ... see sign check below
    # NOTE: corr approximates r = ln_meas - ln_pred = -resid_ref, so the corrected residual is resid_ref + corr.
    T = R.T_K.to_numpy()
    warm, cold = T >= T_WARM, T < T_WARM
    dev, lock = R.role.to_numpy() == "dev", R.role.to_numpy() == "locked"
    binary, tern = R.n_comp.to_numpy() == 2, R.n_comp.to_numpy() == 3
    E_dev = {"refutas": resid_ref, **{k: resid_ref + oof[k] for k in MODELS}}
    readings = {"seed": SEED, "n_boot": N_BOOT, "n_rows": n, "n_systems": int(R.system_id.nunique()), "dev_systems": len(dev_sys),
                "locked_systems": len(lock_sys), "train_pair_rows": int(len(Ptr)), "train_systems": n_train_sys, "class_pair_training_systems": cp_counts}
    xx_row = np.array([x[0] * x[1] if len(x) == 2 else np.nan for x in R.x])
    ep = binary & (xx_row < XMIN)
    readings["endpoint"] = dict(rows=int(ep.sum()), systems=int(R.system_id[ep].nunique()),
                                share_of_binary_systems=float(R.system_id[ep].nunique() / R.system_id[binary].nunique()),
                                mean_resid_ln=float(resid_ref[ep].mean()) if ep.any() else None,
                                mean_abs_resid_ln=float(np.abs(resid_ref[ep]).mean()) if ep.any() else None)

    r1 = {m: gate(R, E_dev, dev & warm, "refutas", m) for m in MODELS}
    r1["binary_only_ml"] = gate(R, E_dev, dev & warm & binary)
    r1["pass"] = bool(r1["ml"]["mae_b"] < r1["ml"]["mae_a"] and r1["ml"]["lo"] > 0)
    readings["R-MLV-1"] = r1
    print("R-MLV-1 ml", json.dumps(r1["ml"]), "pass", r1["pass"], flush=True)

    def margins(E, mask_role):
        out = {}
        for b in ("T_lt_263K", "T_263_283K", "T_ge_283K"):
            msk = mask_role & (R.t_bin.to_numpy() == b)
            d = pd.DataFrame({"s": sysarr[msk], "e": np.abs(E[msk])}).groupby("s").e.max()
            v, ng, k = conformal(d.to_numpy()) if len(d) >= 9 else (float("inf"), len(d), None)
            out[b] = dict(margin=v, n_dev_systems=int(ng), k=k)
        return out
    readings["R-MLV-6"] = dict(ml=margins(E_dev["ml"], dev), refutas=margins(E_dev["refutas"], dev))
    readings["R-MLV-6"]["ml_ge283_le_0.186"] = bool(readings["R-MLV-6"]["ml"]["T_ge_283K"]["margin"] <= 0.186)
    (OUT / "readings_dev.json").write_text(json.dumps(readings, indent=1, default=float))

    if args.dev_only:
        print("dev-only: cold dev(oof)", json.dumps(gate(R, E_dev, dev & cold)))
        print("margins", json.dumps(readings["R-MLV-6"], default=float))
        print("endpoint", readings["endpoint"])
        return

    # ---- locked, scored once
    Pl = P[P.row.isin(np.where(lock)[0])]
    pl = predict_all(model, base, Pl)
    E_all = {"refutas": resid_ref}
    for k in MODELS:
        E_all[k] = resid_ref + np.where(lock, corr_by_row(Pl, pl[k], n), oof[k])
    rl = {m: gate(R, E_all, lock & warm, "refutas", m) for m in MODELS}
    per_ref, per_ml = macro(R, E_all["refutas"], lock & warm), macro(R, E_all["ml"], lock & warm)
    worse = float(((per_ml - per_ref) > 0.02).mean())
    rl["share_systems_worse_gt_0.02"] = worse
    rl["binary_only_ml"] = gate(R, E_all, lock & warm & binary)
    rl["pass"] = bool(rl["ml"]["mae_b"] < rl["ml"]["mae_a"] and rl["ml"]["lo"] > 0 and worse <= 0.25)
    readings["R-MLV-2"] = rl
    print("R-MLV-2 ml", json.dumps(rl["ml"]), "worse", worse, "pass", rl["pass"], flush=True)

    tern_rows = np.where(tern)[0]
    allin = np.array([all(p in (final_pairs if lock[r] else fold_pairs[r]) for p in pairs_of(R.ii.iat[r])) for r in tern_rows])
    ge = {}

    def part(mm):
        return dict(rows=int(mm.sum()), systems=int(R.system_id[mm].nunique()), mae_refutas=float(np.abs(E_all["refutas"][mm]).mean()),
                    mae_ml=float(np.abs(E_all["ml"][mm]).mean()), mae_const=float(np.abs(E_all["const"][mm]).mean()),
                    mae_classpair=float(np.abs(E_all["classpair"][mm]).mean()))
    tm = tern.copy()
    ge["all"] = part(tm)
    ge["all"]["macro"] = gate(R, E_all, tm)
    for lab, sel in (("all_pairs_in_training", allin), ("some_pair_unseen", ~allin)):
        mm = np.zeros(n, bool)
        mm[tern_rows[sel]] = True
        if mm.any():
            ge[lab] = part(mm)
    for lab, rm in (("dev_oof", dev & tern), ("locked", lock & tern)):
        if rm.any():
            ge[lab] = part(rm)
    ge["pass"] = bool(ge["all"]["mae_ml"] <= ge["all"]["mae_refutas"])
    readings["R-MLV-3"] = ge

    cm = cold & (dev | lock)
    c4 = gate(R, E_all, cm)
    ur_ref, ur_ml = float((E_all["refutas"][cm] < 0).mean()), float((E_all["ml"][cm] < 0).mean())
    c4.update(under_pred_rate_refutas=ur_ref, under_pred_rate_ml=ur_ml, rows=int(cm.sum()), systems_dev=int(R.system_id[cm & dev].nunique()),
              systems_locked=int(R.system_id[cm & lock].nunique()), sources=sorted(set(R.source[cm])))
    c4["pass"] = bool(c4["mae_b"] < c4["mae_a"] and c4["lo"] > 0 and ur_ml <= ur_ref)
    for b in ("T_lt_263K", "T_263_283K"):
        mb = cm & (R.t_bin.to_numpy() == b)
        if mb.any():
            c4[b] = gate(R, E_all, mb)
    readings["R-MLV-4"] = c4

    cov = {}
    for k in ("ml", "refutas"):
        cov[k] = {}
        for b, e in readings["R-MLV-6"][k].items():
            msk = lock & (R.t_bin.to_numpy() == b)
            d = pd.DataFrame({"s": sysarr[msk], "e": np.abs(E_all[k][msk])}).groupby("s").e.max()
            cov[k][b] = dict(locked_systems=int(len(d)), covered=int((d <= e["margin"]).sum()) if np.isfinite(e["margin"]) else None)
    readings["R-MLV-6"]["locked_coverage"] = cov
    adopt = bool(r1["pass"] and rl["pass"])
    readings["adoption"] = dict(ships_variant_ml_pair_residual=adopt, cold_gate_passed=bool(c4["pass"]),
                                promotion_offered=bool(adopt and ge["pass"] and readings["R-MLV-6"]["ml_ge283_le_0.186"]))

    rows_out = R[["system_id", "source", "T_K", "t_bin", "n_comp", "role"]].copy()
    for k in ("refutas",) + MODELS:
        rows_out[f"resid_{k}"] = E_all[k]
    rows_out.to_csv(OUT / "rows_scored.csv", index=False, float_format="%.6g")
    inputs = (RES_CSV, SPLIT, ROOT / "data/unified/measurements_dedup.csv", ROOT / "data/real_fuels/v5/boehm2022_viscosity.csv",
              ROOT / "data/unified/compounds.csv")
    meta = dict(model="LightGBM pair g(d_i,d_j,T); ln eta_ML = ln eta_Refutas + sum_{i<j} x_i x_j g", features=vp.FEATURES, class_pairs=vp.CLASS_PAIRS,
                lgb_params=vp.LGB_PARAMS, winsor_ln=WINS, xmin=XMIN, seed=SEED, T_min_K=(100.0 if c4["pass"] else T_WARM),
                class_pair_training_systems=cp_counts, grid_sha256=GRID_SHA,
                input_sha256={str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in inputs},
                cv_scores=dict(dev_oof_warm_macro_mae_refutas=r1["ml"]["mae_a"], dev_oof_warm_macro_mae_ml=r1["ml"]["mae_b"],
                               locked_warm_macro_mae_refutas=rl["ml"]["mae_a"], locked_warm_macro_mae_ml=rl["ml"]["mae_b"]),
                training=dict(systems=n_train_sys, pair_rows=int(len(Ptr))), training_scope="lab binaries, 264 dev systems",
                margin_ln=readings["R-MLV-6"]["ml"], adoption=readings["adoption"])
    (MODEL_DIR / "meta.json").write_text(json.dumps(meta, indent=1, default=float))
    (OUT / "readings.json").write_text(json.dumps(readings, indent=1, default=float))
    print(json.dumps({k: readings[k] for k in ("R-MLV-1", "R-MLV-2", "R-MLV-3", "R-MLV-4", "R-MLV-6", "adoption")}, indent=1, default=float))


if __name__ == "__main__":
    main()
