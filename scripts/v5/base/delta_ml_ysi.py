"""v5 delta_ml_ysi: learned group-contribution YSI gap-filler, nested grouped CV and the declared readings.

Declared readings: results/v5/BASE/delta_ml_ysi/README.md (committed before this script existed).
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_ml_ysi.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "v5" / "BASE" / "delta_ml_ysi"
MODEL_DIR = ROOT / "models" / "v5" / "ml_ysi"
BASE = ROOT / "results" / "v5" / "BASE"
V5 = ROOT / "data" / "real_fuels" / "v5"
SEED = 42
NBOOT = 2000

from fuelblend.ml import ysi_model as YM                                           # noqa: E402

warnings.simplefilter("ignore")


def _load_script(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------------ models and baselines
def lgbm_fit_predict(Xtr, ytr, Xte):
    import lightgbm as lgb
    m = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.05, num_leaves=7, min_child_samples=5, subsample=0.8, subsample_freq=1,
                          colsample_bytree=0.8, reg_lambda=1.0, random_state=SEED, deterministic=True, force_row_wise=True, n_jobs=1, verbose=-1)
    m.fit(Xtr, ytr, sample_weight=YM.weights(ytr))
    return m.predict(Xte)


def class_mean_predict(cls_tr, ytr, cls_te):
    g = float(np.mean(ytr))
    cm = {c: float(np.mean(ytr[cls_tr == c])) for c in np.unique(cls_tr)}
    return np.array([cm.get(c, g) for c in cls_te])


def class_linear_predict(cls_tr, nc_tr, ytr, cls_te, nc_te):
    g = float(np.mean(ytr))
    cm = {c: float(np.mean(ytr[cls_tr == c])) for c in np.unique(cls_tr)}
    lin = {}
    for c in ("paraffin", "iso", "cyclo"):
        m = cls_tr == c
        if m.sum() >= 3 and len(np.unique(nc_tr[m])) >= 2:
            b, a = np.polyfit(nc_tr[m], ytr[m], 1)
            lin[c] = (a, b)
    out = []
    for c, n in zip(cls_te, nc_te):
        if c in lin:
            out.append(max(0.0, lin[c][0] + lin[c][1] * n))
        else:
            out.append(cm.get(c, g))
    return np.array(out)


def boot_diff(errA, errB, groups, rng):
    """95 % interval of mean(errA - errB) with groups resampled."""
    ug, inv = np.unique(groups, return_inverse=True)
    d = np.asarray(errA, float) - np.asarray(errB, float)
    s = np.bincount(inv, d, len(ug))
    n = np.bincount(inv, minlength=len(ug)).astype(float)
    idx = rng.integers(0, len(ug), (NBOOT, len(ug)))
    est = s[idx].sum(1) / n[idx].sum(1)
    return float(d.mean()), float(np.percentile(est, 2.5)), float(np.percentile(est, 97.5))


def conf_q(scores, alpha=0.10):
    s = np.sort(np.asarray(scores, float))
    k = math.ceil((len(s) + 1) * (1 - alpha))
    return float(s[min(k, len(s)) - 1])


def rel(y, p):
    return np.abs(np.asarray(p) - np.asarray(y)) / np.asarray(y)


# ------------------------------------------------------------------------------------------------ readings 1 / 1x
def run_cv(df, X, y, groups, cls, nc):
    pred = {k: np.full(len(y), np.nan) for k in ("ridge", "lgbm", "const", "classmean", "classlin")}
    alphas = []
    for tr, te in GroupKFold(5).split(X, y, groups):
        m = YM.fit_ridge(X[tr], y[tr], groups[tr])
        alphas.append(m.alpha)
        pred["ridge"][te] = m.predict(X[te])
        pred["lgbm"][te] = lgbm_fit_predict(X[tr], y[tr], X[te])
        pred["const"][te] = y[tr].mean()
        pred["classmean"][te] = class_mean_predict(cls[tr], y[tr], cls[te])
        pred["classlin"][te] = class_linear_predict(cls[tr], nc[tr], y[tr], cls[te], nc[te])
    return pred, alphas


def reading_1(df, y, groups, cls, pred, alphas, rng):
    ae = {k: np.abs(v - y) for k, v in pred.items()}
    re = {k: rel(y, v) for k, v in pred.items()}
    out = dict(n=len(y), n_groups=int(len(np.unique(groups))), alpha_per_fold=alphas,
               mae={k: float(v.mean()) for k, v in ae.items()}, median_rel_err={k: float(np.median(v)) for k, v in re.items()},
               mean_rel_err={k: float(v.mean()) for k, v in re.items()})
    gains = {}
    for k in ("const", "classmean", "classlin", "lgbm"):
        m, lo, hi = boot_diff(ae[k], ae["ridge"], groups, rng)
        gains[k] = dict(mae_baseline_minus_ridge=m, lo95=lo, hi95=hi)
    out["gain_vs"] = gains
    out["per_class_median_rel_err_ridge"] = {c: float(np.median(re["ridge"][cls == c])) for c in np.unique(cls)}
    out["per_class_n"] = {c: int((cls == c).sum()) for c in np.unique(cls)}
    p1 = out["median_rel_err"]["ridge"] <= 0.20 and gains["classmean"]["lo95"] > 0
    out["pass"] = bool(p1)
    out["rule"] = "median relative error <= 0.20 AND bootstrap 95% lower bound of (MAE_classmean - MAE_ridge) > 0"
    return out


def reading_1x(df, X, y, groups, cls, nc, smi, rng):
    tr = nc <= 10
    nonar_ring = np.array([YM.is_nonaromatic_ring(s) for s in smi])
    sat = np.array([YM.is_saturated_acyclic_or_cyclic(s) for s in smi])
    cyc_max = nc[tr & nonar_ring].max()
    te = ((sat & (nc >= 11) & (nc <= 12)) | (nonar_ring & (nc > cyc_max))) & ~tr
    te |= ((sat & (nc >= 11) & (nc <= 12)))
    te &= ~tr
    m = YM.fit_ridge(X[tr], y[tr], groups[tr])
    p = m.predict(X[te])
    pl = class_linear_predict(cls[tr], nc[tr], y[tr], cls[te], nc[te])
    pm = class_mean_predict(cls[tr], y[tr], cls[te])
    pg = lgbm_fit_predict(X[tr], y[tr], X[te])
    yt = y[te]
    re = rel(yt, p)
    sc = np.abs(yt - p) / (p + 50.0)
    out = dict(n_train=int(tr.sum()), n_test=int(te.sum()), alpha=m.alpha, cyclic_train_max_nC=int(cyc_max),
               test_compounds=[dict(smiles=s, nC=int(n), cls=c, ysi=float(a), ridge=float(b), classlin=float(c2), classmean=float(c3), lgbm=float(d))
                               for s, n, c, a, b, c2, c3, d in zip(smi[te], nc[te], cls[te], yt, p, pl, pm, pg)],
               median_rel_err=dict(ridge=float(np.median(re)), classlin=float(np.median(rel(yt, pl))), classmean=float(np.median(rel(yt, pm))),
                                   lgbm=float(np.median(rel(yt, pg)))),
               mae=dict(ridge=float(np.abs(yt - p).mean()), classlin=float(np.abs(yt - pl).mean()), classmean=float(np.abs(yt - pm).mean()),
                        lgbm=float(np.abs(yt - pg).mean())))
    g = np.arange(te.sum())
    d, lo, hi = boot_diff(np.abs(yt - pl), np.abs(yt - p), g, rng)
    out["classlin_minus_ridge_mae"] = dict(mean=d, lo95=lo, hi95=hi)
    out["linear_rule_better_excluding_zero"] = bool(hi < 0)
    out["pass"] = bool(out["median_rel_err"]["ridge"] <= 0.15)
    out["rule"] = "median relative error <= 0.15"
    out["q_extrap"] = conf_q(sc) if len(sc) else float("nan")
    out["signed_err_ridge_mean"] = float((p - yt).mean())
    return out


# ------------------------------------------------------------------------------------------------ final model
def final_model(df, X, y, groups, cls, nc, smi, q_oof, q_extrap):
    m = YM.fit_ridge(X, y, groups)
    class_max = {c: int(nc[cls == c].max()) for c in np.unique(cls)}
    scaff = sorted({s for s in (YM.scaffold_generic(x) for x in smi) if s})
    meta = dict(ridge=m.to_dict(), features=YM.FEATURES, class_max_nC=class_max, train_scaffolds=scaff, q_oof=q_oof, q_extrap=q_extrap,
                measured={s: float(v) for s, v in zip(smi, y)}, strict_carbon_margin=YM.STRICT_CARBON_MARGIN,
                training_input=str(YM.YSI_CSV.relative_to(ROOT)).replace("\\", "/"), training_input_sha256=sha(YM.YSI_CSV), n_train=int(len(y)))
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    (MODEL_DIR / "model.json").write_text(json.dumps(meta, indent=1, sort_keys=True), encoding="utf-8")
    YM._load.cache_clear()
    return m, meta


# ------------------------------------------------------------------------------------------------ fuels (R-MLY-2, 3, 5)
def fuel_components():
    cb = _load_script("cells_t2_bulk_mod", "scripts/v5/base/cells_t2_bulk.py")
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.eval.composition_route import kind_bins, map_bins
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    t = pd.read_csv(BASE / "t2_targets.csv")
    fu = pd.read_csv(BASE / "t2_fuels.csv").set_index("entity")
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    ents = sorted(t[t.property == "smoke_point_mm"].entity.unique())
    out = {}
    for ent in ents:
        key = fu.loc[ent, "source"]
        fid = t[t.entity == ent].fuel_id.iloc[0]
        fid = int(fid) if key == "ajftd" else fid
        c = comps[key]
        mf = map_bins(kind_bins(key, c[c.fuel_id == fid]), reps)
        out[ent] = [(Chem.MolToSmiles(Chem.MolFromSmiles(reg.smiles[int(i)]), isomericSmiles=False), reg.ik[int(i)], float(w)) for i, w in zip(mf.idx, mf.w)]
    return cb, out


def ysi_route(smi, w, ysi_vals, scale_ml=None, is_ml=None):
    """Smoke point from mass fractions w and per-component YSI (NaN components dropped and renormalised)."""
    from fuelblend.baselines.sooting import smoke_point_from_ysi
    w = np.asarray(w, float)
    y = np.asarray(ysi_vals, float).copy()
    if scale_ml is not None:
        y = np.where(is_ml, y * scale_ml, y)
    ok = np.isfinite(y)
    mw = np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smi])
    ww = w[ok] / w[ok].sum()
    x = (ww / mw[ok]) / (ww / mw[ok]).sum()
    return float(smoke_point_from_ysi(float(np.dot(x, y[ok])), float(np.dot(x, mw[ok]))))


def reading_fuels(cb, fuel_comp, model, meta, rng):
    from fuelblend.api.predict_route import predict_fuel
    from fuelblend.eval.group_stats import logo_constant
    rows, per = [], {}
    for ent, comp in fuel_comp.items():
        smi = [c[0] for c in comp]
        w = np.array([c[2] for c in comp])
        pr = YM.predict_ysi(smi)
        src = pr.source.to_numpy()
        ysi = pr.ysi.to_numpy(float)
        meas = src == "measured"
        strict = src == "ml_predicted"
        ext = src == "ml_predicted_extrapolated"
        cov_own, cov_strict, cov_full = w[meas].sum(), w[meas | strict].sum(), w[meas | strict | ext].sum()
        is_ml = strict | ext
        sp = ysi_route(smi, w, ysi) if cov_full >= 0.90 else np.nan
        d_ml = np.nan
        if np.isfinite(sp):
            a, b = ysi_route(smi, w, ysi, 0.7, is_ml), ysi_route(smi, w, ysi, 1.3, is_ml)
            d_ml = max(abs(a - sp), abs(b - sp))
        per[ent] = dict(entity=ent, n_comp=len(comp), cov_own=cov_own, cov_strict=cov_strict, cov_full=cov_full, mass_ml_strict=w[strict].sum(),
                        mass_ml_extrap=w[ext].sum(), smoke_ml_filled=sp, sens_mm=d_ml, evaluable_full=bool(cov_full >= 0.90),
                        evaluable_strict=bool(cov_strict >= 0.90), evaluable_own=bool(cov_own >= 0.90))
        try:
            r = predict_fuel([(c[1], c[2]) for c in comp], composition_kind="reconstructed", fuel_type="jet")
            per[ent]["smoke_bulk"] = r["smoke_point_mm"]["value"]
        except Exception as e:                                                     # noqa: BLE001
            per[ent]["smoke_bulk"] = np.nan
            per[ent]["bulk_error"] = f"{type(e).__name__}: {e}"
    pf = pd.DataFrame(per.values())
    r2 = dict(n_fuels=len(pf), evaluable_own=int(pf.evaluable_own.sum()), evaluable_full=int(pf.evaluable_full.sum()),
              evaluable_strict=int(pf.evaluable_strict.sum()), pass_full_ge_32=bool(pf.evaluable_full.sum() >= 32),
              mean_cov_own=float(pf.cov_own.mean()), mean_cov_strict=float(pf.cov_strict.mean()), mean_ml_mass_extrapolated=float(pf.mass_ml_extrap.mean()),
              n_fuels_extrap_mass_gt_0p10=int((pf.mass_ml_extrap > 0.10).sum()))
    # R-MLY-3: BASE scoring subset
    d = pd.read_csv(BASE / "rows_t2_bulk.csv")
    d = d[(d.property == "smoke_point_mm") & (d.model == "ysi_route")].merge(pf[["entity", "smoke_ml_filled", "smoke_bulk", "sens_mm"]], on="entity")
    hd = cb.one_per_group(d[~d.low_coverage]).copy()
    hd = hd[np.isfinite(hd.smoke_ml_filled)]
    const = logo_constant(hd.measured, hd.group, pd.Series(["-"] * len(hd), index=hd.index))
    e_ml, e_b, e_c = hd.smoke_ml_filled - hd.measured, hd.smoke_bulk - hd.measured, np.asarray(const) - hd.measured.to_numpy()
    r3 = dict(n_groups=int(len(hd)), mae_ml_ysi=float(e_ml.abs().mean()), mae_bulk=float(np.nanmean(np.abs(e_b))), mae_logo_constant=float(np.abs(e_c).mean()),
              bias_ml_ysi=float(e_ml.mean()), bias_bulk=float(np.nanmean(e_b)),
              ml_minus_const_mae_boot=boot_diff(np.abs(e_ml.to_numpy()), np.abs(e_c), hd.group.to_numpy(), rng),
              ml_minus_bulk_mae_boot=boot_diff(np.abs(e_ml.to_numpy()), np.abs(e_b.to_numpy()), hd.group.to_numpy(), rng),
              measured_mean=float(hd.measured.mean()), pred_ml_range=[float(hd.smoke_ml_filled.min()), float(hd.smoke_ml_filled.max())],
              note="descriptive; these fuels were already read in BASE")
    s = pf.sens_mm.dropna()
    r5 = dict(n=int(len(s)), frac_le_1mm=float((s <= 1.0).mean()), max_mm=float(s.max()), median_mm=float(s.median()),
              claim_fills_immaterial=bool((s <= 1.0).mean() >= 0.90))
    return pf, r2, r3, r5


# ------------------------------------------------------------------------------------------------ R-MLY-4
def reading_mixtures(model, meta, rng):
    ds = _load_script("delta_nf_mod", "scripts/v5/base/delta_newdata_freeze_smoke.py")
    from fuelblend.baselines.sooting import mixture_mw, smoke_point_from_ysi
    from fuelblend.data import loaders_smoke_lit as LS
    ctx = ds.get_ctx()
    smi_map = ds.smiles_map()
    d = LS.load_smoke_lit()
    d = d[d.property == "smoke_point_mm"].copy()
    d["scope"] = [ds.smoke_scope(r) for r in d.itertuples()]
    keep = d[~d.qualifier.isin(["reference", "reported_mean"]) & (d.scope == "smoke_mixtures")].copy()
    keep["key"] = [(r.source_id, r.system_id, r.fractions, r.fraction_basis, r.scope) for r in keep.itertuples()]
    ref = pd.read_csv(BASE / "delta_newdata_freeze_smoke" / "rows_smoke.csv")
    ref = ref[ref.scope == "smoke_mixtures"].copy()
    rows = []
    for k, g in keep.groupby(keep.key.astype(str), sort=False):
        r = g.iloc[0]
        iks, fr, basis = json.loads(r.comp_inchikeys), json.loads(r.fractions), r.fraction_basis
        smi = [smi_map[i] for i in iks]
        if basis == "volume":
            idx = [ctx.reg.ik2i.get(i) for i in iks]
            if any(i is None for i in idx):
                continue
            _, rho = ctx.prov.lookup(np.array(idx, np.int64)[None, :], ds.T_REF_K)
            fr = list(np.array(fr) * np.asarray(rho)[0])
            basis = "mass"
        fr = np.array(fr, float)
        mw = np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smi])
        x = fr / fr.sum() if basis == "mole" else (fr / mw) / (fr / mw).sum()
        pr = YM.predict_ysi(smi)
        y = pr.ysi.to_numpy(float)
        sp = float(smoke_point_from_ysi(float(np.dot(x, y)), float(np.dot(x, mw)))) if np.all(np.isfinite(y)) else np.nan
        rows.append(dict(source=r.source_id, system_id=r.system_id, sample=" | ".join(sorted(set(g.sample_id))), measured=float(g.value.mean()),
                         ml_filled=sp, n_ml=int((pr.source != "measured").sum()), n_ml_extrap=int((pr.source == "ml_predicted_extrapolated").sum()),
                         x_ml=float(x[(pr.source != "measured").to_numpy()].sum()), cond_lo=np.nan, cond_hi=np.nan))
    mx = pd.DataFrame(rows).merge(ref[["source", "system_id", "sample", "pred_ysi", "pred_bulk", "group"]], on=["source", "system_id", "sample"], how="left")
    mx["has_ysi_own"] = np.isfinite(mx.pred_ysi.astype(float))
    ev = mx[np.isfinite(mx.ml_filled)]
    own = ev[ev.has_ysi_own]
    out = dict(n_points=int(len(mx)), n_evaluable_own=int(mx.has_ysi_own.sum()), n_evaluable_ml=int(len(ev)),
               mae_ysi_without_ml_on_own=float((own.pred_ysi - own.measured).abs().mean()) if len(own) else None,
               bias_ysi_without_ml_on_own=float((own.pred_ysi - own.measured).mean()) if len(own) else None,
               mae_ysi_ml_filled_all=float((ev.ml_filled - ev.measured).abs().mean()), bias_ysi_ml_filled_all=float((ev.ml_filled - ev.measured).mean()),
               mae_ysi_ml_filled_same_subset=float((own.ml_filled - own.measured).abs().mean()) if len(own) else None,
               mae_bulk_designated_all=float((ev.pred_bulk - ev.measured).abs().mean()), mae_bulk_same_subset=float((own.pred_bulk - own.measured).abs().mean()) if len(own) else None,
               n_with_extrapolated_fill=int((ev.n_ml_extrap > 0).sum()))
    return out, mx


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    df = YM.load_training()
    smi = df.smiles.to_numpy()
    y = df.ysi.to_numpy(float)
    X = np.array([YM.features(s) for s in smi])
    groups = np.array([YM.group_key(s) for s in smi])
    cls = np.array([YM.mol_class(s) for s in smi])
    nc = X[:, 0].astype(int)
    print(f"n={len(y)} compounds, {len(np.unique(groups))} groups, classes {pd.Series(cls).value_counts().to_dict()}")

    pred, alphas = run_cv(df, X, y, groups, cls, nc)
    r1 = reading_1(df, y, groups, cls, pred, alphas, rng)
    pd.DataFrame(dict(smiles=smi, ysi=y, group=groups, cls=cls, nC=nc, **{f"pred_{k}": v for k, v in pred.items()})).to_csv(OUT / "oof_predictions.csv", index=False)
    q_oof = conf_q(np.abs(y - pred["ridge"]) / (pred["ridge"] + 50.0))
    cover = int((np.abs(y - pred["ridge"]) <= q_oof * (pred["ridge"] + 50.0)).sum())
    r1["conformal_q_oof"] = q_oof
    r1["conformal_oof_cover_counts"] = f"{cover}/{len(y)}"

    r1x = reading_1x(df, X, y, groups, cls, nc, smi, rng)
    model, meta = final_model(df, X, y, groups, cls, nc, smi, q_oof, r1x["q_extrap"] if np.isfinite(r1x["q_extrap"]) else 0.0)
    (MODEL_DIR / "meta.json").write_text(json.dumps(dict(
        model="ridge group-contribution on 16 RDKit counts, StandardScaler, weight 1/(YSI+50)^2", alpha=model.alpha, features=YM.FEATURES,
        coef=dict(zip(YM.FEATURES, model.coef_.tolist())), training_input_sha256=meta["training_input_sha256"], n_train=meta["n_train"],
        cv=dict(median_rel_err_ridge=r1["median_rel_err"]["ridge"], mae_ridge=r1["mae"]["ridge"], mae_classmean=r1["mae"]["classmean"], pass_R_MLY_1=r1["pass"]),
        extrapolation=dict(median_rel_err_ridge=r1x["median_rel_err"]["ridge"], n_test=r1x["n_test"], pass_R_MLY_1x=r1x["pass"]),
        seed=SEED, note="ridge fit minimises sum sw*(y-f)^2 + (alpha/n)*|beta|^2 on standardised features, sw = w/sum(w)"), indent=1), encoding="utf-8")

    cb, fuel_comp = fuel_components()
    pf, r2, r3, r5 = reading_fuels(cb, fuel_comp, model, meta, rng)
    pf.to_csv(OUT / "rows_fuels.csv", index=False)
    r4, mx = reading_mixtures(model, meta, rng)
    mx.to_csv(OUT / "rows_smoke_mixtures_ml.csv", index=False)
    pd.DataFrame(r1x["test_compounds"]).to_csv(OUT / "rows_extrapolation_test.csv", index=False)
    g = pd.DataFrame([dict(smiles=s, **dict(zip(["ysi", "lo", "hi", "source"], [r.ysi, r.lo, r.hi, r.source])))
                      for s in sorted({c[0] for v in fuel_comp.values() for c in v}) for r in [YM.predict_ysi([s]).iloc[0]]])
    g.to_csv(OUT / "gap_molecules_predictions.csv", index=False)

    adoption = dict(R_MLY_1=r1["pass"], R_MLY_1x=r1x["pass"], offer_variant_ysi_route_ml_filled=bool(r1["pass"] and r1x["pass"]),
                    designated_smoke_model="bulk regression (unchanged, user decision)")
    readings = {"R-MLY-1": r1, "R-MLY-1x": {k: v for k, v in r1x.items() if k != "test_compounds"}, "R-MLY-2": r2, "R-MLY-3": r3, "R-MLY-4": r4, "R-MLY-5": r5,
                "adoption": adoption, "inputs": dict(ysi_hc_sha256=meta["training_input_sha256"], seed=SEED, n_boot=NBOOT)}
    (OUT / "readings.json").write_text(json.dumps(readings, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: v for k, v in readings.items() if k not in ("inputs",)}, indent=1, default=float)[:9000])


if __name__ == "__main__":
    main()
