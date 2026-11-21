"""Track A step 3: baselines (a), (c1), (c2), (d), (e) and the monotone LightGBM residual, LOGO on dev.
Holdout is NOT touched here (step4_holdout.py scores it once)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import RidgeCV

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

from fuelblend.eval import composition_route as CR  # noqa: E402

CLIP = (-100.0, -20.0)
LGB_COLS = ["n_heavy_pct", "n_pct", "meanC", "arom_pct"]
LGB_MONO = [1, 1, 1, 0]


def features(ent: dict) -> dict:
    b = ent["bins"]
    tot = b.mass.sum()
    n = b[b.cls == "n_paraffin"]
    f = dict(n_pct=100 * n.mass.sum() / tot, n_heavy_pct=100 * n[n.n_carbon >= 14].mass.sum() / tot,
             arom_pct=100 * b[b.cls == "aromatic"].mass.sum() / tot, meanC=float((b.n_carbon * b.mass).sum() / tot))
    for c in range(9, 21):
        f[f"nC{c}"] = 100 * n[n.n_carbon == c].mass.sum() / tot
    return f


def c1_pred(ent: dict, tdes: dict) -> float:
    m = ent["mapped"]
    tm = np.array([(tdes.get(k) or {}).get("tm_K", np.nan) for k in m.ik], float)
    w = m.frac.to_numpy(float)
    ok = np.isfinite(tm)
    return float((w[ok] * tm[ok]).sum() / w[ok].sum() - 273.15) if ok.any() else np.nan


def fp_const(tr, te):
    return np.full(len(te), tr.freezing_C.mean())


def fp_ridge(tr, te):
    cols = ["n_pct", "n_heavy_pct", "arom_pct", "meanC"]
    mu, sd = tr[cols].mean(), tr[cols].std().replace(0, 1)
    m = RidgeCV(alphas=[0.1, 1, 10, 100, 1000]).fit((tr[cols] - mu) / sd, tr.freezing_C)
    return m.predict((te[cols] - mu) / sd)


def fp_control(tr, te):
    cands = [f"nC{c}" for c in range(9, 21)] + ["arom_pct"]
    y = tr.freezing_C.to_numpy()
    best, best_mae, best_inc = None, np.inf, True
    for f in cands:
        x = tr[f].to_numpy()
        if np.ptp(x) < 1e-9:
            continue
        inc = bool(np.corrcoef(x, y)[0, 1] >= 0)
        err = []
        for i in range(len(tr)):
            k = np.arange(len(tr)) != i
            ir = IsotonicRegression(increasing=inc, out_of_bounds="clip").fit(x[k], y[k])
            err.append(abs(ir.predict([x[i]])[0] - y[i]))
        if np.mean(err) < best_mae:
            best, best_mae, best_inc = f, float(np.mean(err)), inc
    ir = IsotonicRegression(increasing=best_inc, out_of_bounds="clip").fit(tr[best].to_numpy(), y)
    return ir.predict(te[best].to_numpy())


def make_lgbm(base_col: str):
    import lightgbm as lgb

    def fp(tr, te):
        btr = tr[base_col].clip(*CLIP).to_numpy()
        bte = te[base_col].clip(*CLIP).to_numpy()
        m = lgb.LGBMRegressor(n_estimators=60, learning_rate=0.05, num_leaves=4, min_child_samples=8, subsample=1.0,
                              monotone_constraints=LGB_MONO, verbose=-1, random_state=0)
        m.fit(tr[LGB_COLS], tr.freezing_C.to_numpy() - btr)
        return bte + m.predict(te[LGB_COLS])
    return fp


def logo(df: pd.DataFrame, fit_predict) -> np.ndarray:
    out = np.full(len(df), np.nan)
    for i in range(len(df)):
        tr = df.drop(df.index[i])
        out[i] = fit_predict(tr, df.iloc[[i]])[0]
    return out


def stats(y, p, const):
    e, ec = p - y, const - y
    rng = np.random.default_rng(0)
    d = np.abs(ec) - np.abs(e)
    boot = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
    lo = float(np.quantile(boot, .025))
    return dict(n=len(y), mae=float(np.abs(e).mean()), bias=float(e.mean()), rmse=float(np.sqrt((e ** 2).mean())), frac_under=float((p < y).mean()),
                mae_const=float(np.abs(ec).mean()), diff_vs_const=float(d.mean()), ci_lo=lo, ci_hi=float(np.quantile(boot, .975)),
                margin90=float(np.quantile(y - p, 0.9)), verdict="skill" if lo > 0 else "no_skill_shown")


def main() -> None:
    ds = pd.read_csv(C.OUT / "baseline_b_preds.csv", dtype={"fuel_id": str})
    ents = C.entities(ds)
    s5 = C.load("s5_v6c", "scripts/v5/s5_freezing_nalkane.py")
    g1 = C.load("g1_v6c", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    from fuelblend.api.predict_route import _tmap_freezing
    from fuelblend.baselines import flash_unifac as fu
    from fuelblend.data.tile_isomers import get_registry_and_reps, tile_mean_tmap
    tdes = _tmap_freezing()
    reg, reps = get_registry_and_reps()
    tm_e = tile_mean_tmap(reg, reps, tdes)
    ds = pd.concat([ds.reset_index(drop=True), pd.DataFrame([features(ents[e]) for e in ds.entity])], axis=1)
    ds["pred_c1"] = [c1_pred(ents[e], tdes) for e in ds.entity]
    ds["pred_e"] = [C.unifac_designated(ents[e], s5, g1, CR, fu, tm_e) for e in ds.entity]
    ds.to_csv(C.OUT / "features_preds.csv", index=False, float_format="%.8g")

    dev = ds[ds.role == "dev"].reset_index(drop=True)
    oof = {}
    for name, fn in [("a_constant", fp_const), ("c2_ridge", fp_ridge), ("d_control_curve", fp_control)]:
        oof[name] = logo(dev, fn)
    for name, col in [("b_unifac_svl", "pred_b"), ("c1_weighted_tm", "pred_c1"), ("e_isomer_tiles", "pred_e")]:
        oof[name] = dev[col].to_numpy(float)
    y = dev.freezing_C.to_numpy()
    base_name = min(("b_unifac_svl", "e_isomer_tiles"), key=lambda k: float(np.nanmean(np.abs(oof[k] - y))))
    base_col = {"b_unifac_svl": "pred_b", "e_isomer_tiles": "pred_e"}[base_name]
    oof["f_lgbm_residual_on_" + base_name] = logo(dev, make_lgbm(base_col))
    const = oof["a_constant"]
    res = []
    for name, p in oof.items():
        for tag, mask in (("all", np.ones(len(y), bool)), ("no_fuel64", (dev.fuel_id != "64").to_numpy())):
            ok = mask & np.isfinite(p)
            res.append(dict(model=name, subset=tag, n_nonfinite=int((mask & ~np.isfinite(p)).sum()), **stats(y[ok], p[ok], const[ok])))
    r = pd.DataFrame(res)
    r.to_csv(C.OUT / "dev_logo_summary.csv", index=False, float_format="%.4f")
    pd.DataFrame(oof).assign(entity=dev.entity, y=y).to_csv(C.OUT / "dev_oof_predictions.csv", index=False, float_format="%.4f")
    json.dump(dict(base=base_name, base_col=base_col), open(C.OUT / "chosen_base.json", "w"))
    pd.set_option("display.width", 200)
    print(r[["model", "subset", "n", "mae", "bias", "mae_const", "ci_lo", "ci_hi", "margin90", "verdict"]].round(2).to_string())


if __name__ == "__main__":
    main()
