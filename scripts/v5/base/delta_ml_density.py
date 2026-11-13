"""v5 delta_ml_density: learned pair excess-volume correction on ideal-volume mixture density.
Declared text: results/v5/BASE/delta_ml_density/README.md.
Usage: python scripts/v5/base/delta_ml_density.py --dev-only   (dev OOF + endpoint share; locked NOT scored)
       python scripts/v5/base/delta_ml_density.py              (full run; locked scored once)
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "v5" / "BASE" / "delta_ml_density"
MODEL_DIR = ROOT / "models" / "v5" / "ml_density_pair"
COMMITTED = ROOT / "results" / "v4" / "L" / "density_lab_residuals.csv"
THERMOML = ROOT / "data" / "processed" / "v2" / "thermoml_clean_hc.csv"
SEED, N_BOOT, WINS, XMIN = 42, 2000, 2.0, 0.01
MODELS = ("ml", "const", "classpair")


def load_by_path(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def regenerate():
    from rdkit import Chem
    from fuelblend.ml import density_pair as dp
    lab = load_by_path("eval_lab_known_composition", "scripts/v4/eval_lab_known_composition.py")
    ev = load_by_path("eval_bulk_props_v4", "scripts/v4/eval_bulk_props_v4.py")
    thermo = pd.read_csv(THERMOML, low_memory=False)
    pdm = ev.PureDensity(pd.DataFrame(columns=["inchikey_norm", "cls", "n_carbon"]))
    d, _ = lab.prepare_density_records(thermo)
    bad = {s for s, iks in zip(d.system_id, d.iks) if any(k not in pdm.t1 for k in iks)}
    desc = {}
    t = thermo[thermo.prop.str.contains("Mass density") & (thermo.prop_phase == "Liquid") & (thermo.n_comp >= 2)]
    for s in t.inchis.unique():
        for inchi in json.loads(s):
            info = lab._inchi_info(inchi)
            if info and info[0] not in desc:
                desc[info[0]] = dp.descriptors_from_mol(Chem.MolFromInchi(inchi))
    rows = []
    for r in d[~d.system_id.isin(bad)].itertuples():
        rho = []
        for k in r.iks:
            val, tier = pdm.rho(k, "unused", 10, r.T_K)
            rho.append(val if tier == 1 else np.nan)
        if not np.all(np.isfinite(rho)):
            continue
        w = np.array(r.w_full)
        rho = np.array(rho)
        M = np.array([desc[k][0][1] for k in r.iks])
        x = (w / M) / np.sum(w / M)
        mm = float(np.sum(x * M))
        pred = 1.0 / float(np.sum(w / rho))
        rows.append(dict(system_id=r.system_id, iks=r.iks, T_K=r.T_K, n_comp=r.n_comp, x=tuple(x), rho_i=tuple(rho), mm=mm,
                         pred=pred, meas=r.measured, v_ideal=mm / pred * 1000.0, v_meas=mm / r.measured * 1000.0))
    R = pd.DataFrame(rows)
    com = pd.read_csv(COMMITTED)
    assert len(com) == len(R), f"row count {len(com)} vs {len(R)}"
    assert (com.system_id.to_numpy() == R.system_id.to_numpy()).all()
    dp_ = np.abs(com.predicted.to_numpy() - R.pred.to_numpy()); print("max abs diff predicted", dp_.max(), "max rel", (dp_ / R.pred.to_numpy()).max())
    assert (dp_ / R.pred.to_numpy()).max() < 1e-6, "predicted differs"
    assert (np.abs(com.measured.to_numpy() - R.meas.to_numpy()) / R.meas.to_numpy()).max() < 1e-6, "measured differs"
    assert np.abs(com.T_K.to_numpy() - R["T_K"].to_numpy()).max() < 1e-5
    R["role"] = com.split_role.to_numpy()
    R["t_bin"] = com.t_bin.to_numpy()
    return R, desc


def build_pairs(R, desc):
    from fuelblend.ml import density_pair as dp
    rid, a, b = [], [], []
    for k, iks in enumerate(R.iks):
        for i in range(len(iks)):
            for j in range(i + 1, len(iks)):
                rid.append(k)
                a.append(i)
                b.append(j)
    rid, a, b = np.array(rid), np.array(a), np.array(b)
    ika = [R.iks.iat[k][i] for k, i in zip(rid, a)]
    ikb = [R.iks.iat[k][j] for k, j in zip(rid, b)]
    sa = np.array([desc[k][0] for k in ika])
    sb = np.array([desc[k][0] for k in ikb])
    ca = [desc[k][1] for k in ika]
    cb = [desc[k][1] for k in ikb]
    ra = np.array([R.rho_i.iat[k][i] for k, i in zip(rid, a)])
    rb = np.array([R.rho_i.iat[k][j] for k, j in zip(rid, b)])
    P = dp.pair_frame(sa, sb, ca, cb, ra, rb, R["T_K"].to_numpy()[rid])
    P["row"] = rid
    P["xx"] = np.array([R.x.iat[k][i] * R.x.iat[k][j] for k, i, j in zip(rid, a, b)])
    return P


def train_set(R, P, systems):
    ridx = P.row.to_numpy()
    m = np.isin(R.system_id.to_numpy()[ridx], list(systems)) & (R.n_comp.to_numpy()[ridx] == 2) & (P.xx.to_numpy() >= XMIN)
    return P[m]


def targets(R, Ptr):
    ve = (R.v_meas.to_numpy() - R.v_ideal.to_numpy())[Ptr.row.to_numpy()]
    g = np.clip(ve, -WINS, WINS) / Ptr.xx.to_numpy()
    w = Ptr.xx.to_numpy() ** 2
    tot = pd.Series(w).groupby(R.system_id.to_numpy()[Ptr.row.to_numpy()]).transform("sum").to_numpy()
    return g, w / tot


def fit_lgb(Ptr, g, w):
    import lightgbm as lgb
    from fuelblend.ml.density_pair import FEATURES, LGB_PARAMS
    m = lgb.LGBMRegressor(**LGB_PARAMS)
    m.fit(Ptr[FEATURES], g, sample_weight=w, categorical_feature=["class_pair"])
    return m


def fit_baselines(Ptr, g, w):
    g0 = float(np.average(g, weights=w))
    cp = pd.DataFrame({"c": Ptr.class_pair.to_numpy(), "g": g, "w": w})
    return g0, {int(c): float(np.average(d.g, weights=d.w)) for c, d in cp.groupby("c")}


def predict_all(model, base, P):
    from fuelblend.ml.density_pair import FEATURES
    g0, tab = base
    return dict(ml=model.booster_.predict(P[FEATURES]), const=np.full(len(P), g0),
                classpair=np.array([tab.get(int(c), g0) for c in P.class_pair]))


def corr_by_row(Pp, gv, n):
    return np.bincount(Pp.row.to_numpy(), weights=Pp.xx.to_numpy() * gv, minlength=n)


def macro(R, err, mask):
    return pd.DataFrame({"s": R.system_id.to_numpy()[mask], "e": np.abs(err[mask])}).groupby("s").e.mean()


def paired_boot(a, b, rng):
    d = (a - b.reindex(a.index)).to_numpy()
    n = len(d)
    bs = np.array([d[rng.integers(0, n, n)].mean() for _ in range(N_BOOT)])
    return dict(n_systems=int(n), mae_a=float(a.mean()), mae_b=float(b.mean()), gain=float(d.mean()),
                lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)))


def gate(R, E, mask, a="ideal", b="ml"):
    return paired_boot(macro(R, E[a], mask), macro(R, E[b], mask), np.random.default_rng(SEED))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dev-only", action="store_true")
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    from fuelblend.ml import density_pair as dp
    from fuelblend.uq import conformal as cf
    R, desc = regenerate()
    print(f"regenerated {len(R)} rows, {R.system_id.nunique()} systems; positional agreement OK", flush=True)
    P = build_pairs(R, desc)
    n = len(R)
    sysarr = R.system_id.to_numpy()
    err_ideal = R.pred.to_numpy() - R.meas.to_numpy()
    role = R.role.to_numpy()
    dev, lock = role == "dev", role == "locked"
    binary, tern = R.n_comp.to_numpy() == 2, R.n_comp.to_numpy() == 3
    dev_rows = np.where(dev)[0]
    from sklearn.model_selection import GroupKFold
    oof = {k: np.zeros(n) for k in MODELS}
    for f, (tr, te) in enumerate(GroupKFold(n_splits=5).split(dev_rows, groups=sysarr[dev_rows])):
        tr_sys = set(sysarr[dev_rows[tr]])
        Ptr = train_set(R, P, tr_sys)
        g, w = targets(R, Ptr)
        model, base = fit_lgb(Ptr, g, w), fit_baselines(Ptr, g, w)
        te_rows = dev_rows[te]
        Pte = P[P.row.isin(te_rows)]
        for k, gv in predict_all(model, base, Pte).items():
            oof[k][te_rows] = corr_by_row(Pte, gv, n)[te_rows]
        print(f"fold {f}: train systems {len(tr_sys)}, train pair rows {len(Ptr)}, test systems {len(set(sysarr[te_rows]))}", flush=True)

    def err_from(corr):
        return R.mm.to_numpy() / (R.v_ideal.to_numpy() + corr) * 1000.0 - R.meas.to_numpy()

    Ptr = train_set(R, P, set(sysarr[dev]))
    g, w = targets(R, Ptr)
    model, base = fit_lgb(Ptr, g, w), fit_baselines(Ptr, g, w)
    cps = pd.DataFrame({"s": sysarr[Ptr.row.to_numpy()], "c": Ptr.class_pair.to_numpy()}).drop_duplicates()
    cp_counts = {dp.CLASS_PAIRS[int(c)]: int(k) for c, k in cps.groupby("c").size().items() if int(c) >= 0}
    print("training systems per class pair:", cp_counts, flush=True)

    E = {"ideal": err_ideal, **{k: err_from(oof[k]) for k in MODELS}}
    readings = dict(seed=SEED, n_boot=N_BOOT, n_rows=n, n_systems=int(R.system_id.nunique()), dev_systems=int(R.system_id[dev].nunique()),
                    locked_systems=int(R.system_id[lock].nunique()), train_pair_rows=int(len(Ptr)), class_pair_training_systems=cp_counts)
    xx_row = np.array([x[0] * x[1] if len(x) == 2 else np.nan for x in R.x])
    ep = dev & binary & (xx_row < XMIN)
    ve = (R.v_meas - R.v_ideal).to_numpy()
    readings["R-MLD-5"] = dict(endpoint_rows_dev=int(ep.sum()), dev_binary_rows=int((dev & binary).sum()),
                               endpoint_row_share=float(ep.sum() / (dev & binary).sum()),
                               endpoint_systems=int(R.system_id[ep].nunique()),
                               share_of_dev_binary_ideal_abs_error_in_endpoint_rows=float(np.abs(err_ideal[ep]).sum() / np.abs(err_ideal[dev & binary]).sum()),
                               mean_abs_err_ideal_endpoint_kg_m3=float(np.abs(err_ideal[ep]).mean()) if ep.any() else None,
                               mean_abs_ve_endpoint_cm3_mol=float(np.abs(ve[ep]).mean()) if ep.any() else None,
                               real_fuel_and_thermoml2020="not run, time")
    r1 = {m: gate(R, E, dev, "ideal", m) for m in MODELS}
    r1["pass"] = bool(r1["ml"]["mae_b"] < r1["ml"]["mae_a"] and r1["ml"]["lo"] > 0)
    readings["R-MLD-1"] = r1
    print("R-MLD-1", json.dumps(r1), flush=True)
    print("R-MLD-5 (dev, pre-locked)", json.dumps(readings["R-MLD-5"]), flush=True)
    (OUT / "readings_dev.json").write_text(json.dumps(readings, indent=1, default=float))
    if args.dev_only:
        return

    Pl = P[P.row.isin(np.where(lock)[0])]
    pl = predict_all(model, base, Pl)
    E_all = {"ideal": err_ideal}
    for k in MODELS:
        E_all[k] = err_from(np.where(lock, corr_by_row(Pl, pl[k], n), oof[k]))
    rl = {m: gate(R, E_all, lock, "ideal", m) for m in MODELS}
    per_i, per_m = macro(R, E_all["ideal"], lock), macro(R, E_all["ml"], lock)
    worse = float(((per_m - per_i) > 0.2).mean())
    rl["share_systems_worse_gt_0.2"] = worse
    rl["pass"] = bool(rl["ml"]["mae_b"] < rl["ml"]["mae_a"] and rl["ml"]["lo"] > 0 and worse <= 0.25)
    readings["R-MLD-2"] = rl
    print("R-MLD-2", json.dumps(rl), flush=True)

    def part(mm):
        return dict(rows=int(mm.sum()), systems=int(R.system_id[mm].nunique()),
                    **{f"mae_{k}": float(np.abs(E_all[k][mm]).mean()) for k in ("ideal",) + MODELS})
    r3 = {"all": part(tern), "dev_oof": part(tern & dev), "locked": part(tern & lock)}
    r3["all"]["macro"] = gate(R, E_all, tern)
    r3["pass"] = bool(r3["all"]["mae_ml"] <= r3["all"]["mae_ideal"])
    readings["R-MLD-3"] = r3

    m = {}
    for k in ("ml", "ideal"):
        m[k] = dict(locked_margin_kg_m3=cf.two_sided_margin(E_all[k][lock], sysarr[lock], 0.10),
                    dev_oof_margin_kg_m3=cf.two_sided_margin(E_all[k][dev], sysarr[dev], 0.10))
        m[k]["locked_coverage_by_dev_margin"] = list(cf.coverage(E_all[k][lock], sysarr[lock], m[k]["dev_oof_margin_kg_m3"], two_sided=True))
        m[k]["locked_coverage_by_locked_margin"] = list(cf.coverage(E_all[k][lock], sysarr[lock], m[k]["locked_margin_kg_m3"], two_sided=True))
    readings["R-MLD-4"] = m
    adopt = bool(r1["pass"] and rl["pass"])
    readings["adoption"] = dict(ships_variant_ml_pair_excess_volume=adopt, promotion_offered=bool(adopt and r3["pass"]))

    rows_out = R[["system_id", "T_K", "t_bin", "n_comp", "role"]].copy()
    for k in ("ideal",) + MODELS:
        rows_out[f"err_{k}"] = E_all[k]
    rows_out.to_csv(OUT / "rows_scored.csv", index=False, float_format="%.6g")
    model.booster_.save_model(str(MODEL_DIR / "model.txt"))
    meta = dict(model="LightGBM pair g(d_i,d_j,T); V_ML = V_ideal + sum_{i<j} x_i x_j g (cm3/mol)", features=dp.FEATURES,
                class_pairs=dp.CLASS_PAIRS, lgb_params=dp.LGB_PARAMS, winsor_cm3_mol=WINS, xmin=XMIN, seed=SEED,
                class_pair_training_systems=cp_counts,
                input_sha256={"results/v4/L/density_lab_residuals.csv": hashlib.sha256(COMMITTED.read_bytes()).hexdigest(),
                              "data/processed/v2/thermoml_clean_hc.csv": hashlib.sha256(THERMOML.read_bytes()).hexdigest()},
                cv_scores=dict(dev_oof_mae_ideal=r1["ml"]["mae_a"], dev_oof_mae_ml=r1["ml"]["mae_b"],
                               locked_mae_ideal=rl["ml"]["mae_a"], locked_mae_ml=rl["ml"]["mae_b"]),
                margin_kg_m3=m, adoption=readings["adoption"], training_scope="lab binaries (ThermoML), 167 dev systems")
    (MODEL_DIR / "meta.json").write_text(json.dumps(meta, indent=1, default=float))
    (OUT / "readings.json").write_text(json.dumps(readings, indent=1, default=float))
    print(json.dumps({k: readings[k] for k in ("R-MLD-3", "R-MLD-4", "adoption")}, indent=1, default=float))


if __name__ == "__main__":
    main()
