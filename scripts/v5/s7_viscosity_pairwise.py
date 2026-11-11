"""v5 S7: pairwise model bootstrap and split-conformal recalibration of the viscosity GP interval.

Reads existing v4 outputs only (no refit).  Readings R-S7-1..4 are declared in results/v5/S7/README.md.
GP predictive sd is recovered from the stored 90% bounds: sd = (ln hi - ln lo) / (2 * 1.645) (symmetric in ln space).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "results" / "v4" / "V"
OUT = ROOT / "results" / "v5" / "S7"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
Z90 = 1.645
V4_SEED = 20261110
MODELS = {"bilinear": "Bilinear pairwise", "ridge": "Structured ridge", "gp": "Structured GP"}
PAIRS = [("bilinear", "ridge"), ("bilinear", "gp"), ("ridge", "gp")]
TCOND = [("T = 20 C", 20.0), ("T = -20 C", -20.0), ("T = -40 C", -40.0)]
CONDS = [("overall (all T)", None)] + TCOND
MIN_DEV_FUELS = 9


def rmse(e):
    return float(np.sqrt(np.mean(np.square(e))))


def conformal_rank(n, level=0.9):
    return math.ceil(round((n + 1) * level, 9))


def conformal_factor(scores, level=0.9):
    """ceil((n+1)*level)-th smallest score (inf if that rank exceeds n)."""
    s = np.sort(np.asarray(scores, float))
    k = conformal_rank(len(s), level)
    return float("inf") if k > len(s) else float(s[k - 1])


def dev_factor(df, ratio_col, cond_T=None):
    """Factor from DEV fuels only (locked rows are filtered out here, whatever the caller passes).
    Score per fuel = max over that fuel's rows (in the condition) of the nonconformity ratio."""
    d = df[~df.locked]
    if cond_T is not None:
        d = d[d.T_C == cond_T]
    return conformal_factor(d.groupby("fuel_id")[ratio_col].max().to_numpy()), d.fuel_id.nunique()


def paired_boot(err_a, err_b, fid, n=5000, seed=0, metric=rmse):
    uc, inv = np.unique(fid, return_inverse=True)
    idx = [np.where(inv == c)[0] for c in range(len(uc))]
    rng = np.random.default_rng(seed)
    d = np.empty(n)
    for b in range(n):
        s = np.concatenate([idx[c] for c in rng.integers(0, len(uc), len(uc))])
        d[b] = metric(err_a[s]) - metric(err_b[s])
    return metric(err_a) - metric(err_b), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    D = pd.read_csv(V4 / "viscosity_predictions_v4.csv")
    lb = pd.read_csv(V4 / "viscosity_leaderboard_v4.csv")
    locked_ids = set(json.load(open(SPLITS))["locked"])
    assert set(D.fuel_id[D.locked]) == {f for f in locked_ids if f.startswith("CRC647_")} and D.fuel_id[~D.locked].nunique() == 45
    y = np.log(D.y.to_numpy())
    err = {k: np.log(D[f"pred_{v}"].to_numpy()) - y for k, v in MODELS.items()}
    err["recal"] = np.log(D["pred_Recalibrated Refutas"].to_numpy()) - y
    lo, hi = np.log(D.gp90_lo.to_numpy()), np.log(D.gp90_hi.to_numpy())
    mu = np.log(D["pred_Structured GP"].to_numpy())
    D["gp_sd"] = (hi - lo) / (2 * Z90)
    assert np.abs((lo + hi) / 2 - mu).max() < 1e-9 and (D.gp_sd > 0).all()
    D["ratio_gp"] = np.abs(mu - y) / D.gp_sd
    D["abs_res_ridge"] = np.abs(err["ridge"])
    D["abs_res_bilinear"] = np.abs(err["bilinear"])
    fid = D.fuel_id.to_numpy()
    splits = [("dev (LOO by fuel)", ~D.locked.to_numpy()), ("locked (untouched)", D.locked.to_numpy())]

    def cmask(cT):
        return np.ones(len(D), bool) if cT is None else (D.T_C.to_numpy() == cT)

    # ---- reproduction check of the logged v4 numbers
    rep = []
    for sname, sm in splits:
        for cname, cT in CONDS:
            mk = sm & cmask(cT)
            row = lb[(lb.model == "Structured GP") & (lb.split == sname) & (lb.condition == cname)].iloc[0]
            inside = (y[mk] >= lo[mk]) & (y[mk] <= hi[mk])
            c_re, w_re = round(float(inside.mean()), 3), round(float(np.mean(hi[mk] - lo[mk])), 3)
            rep.append(dict(check="gp_cov90", split=sname, condition=cname, model="Structured GP",
                            logged=row.gp_cov90, recomputed=c_re, match=bool(abs(row.gp_cov90 - c_re) < 1e-9)))
            rep.append(dict(check="gp_mean_width90_log", split=sname, condition=cname, model="Structured GP",
                            logged=row.gp_mean_width90_log, recomputed=w_re,
                            match=bool(abs(row.gp_mean_width90_log - w_re) < 1e-9)))
            for k, v in MODELS.items():
                lr = lb[(lb.model == v) & (lb.split == sname) & (lb.condition == cname)].iloc[0]
                d, l, h = paired_boot(err[k][mk], err["recal"][mk], fid[mk], n=2000, seed=V4_SEED)
                s_re = f"[{l:.4f}, {h:.4f}]"
                rep.append(dict(check="d_log_rmse_vs_recal", split=sname, condition=cname, model=v,
                                logged=f"{lr.d_log_rmse_vs_recal:.4f} {lr.d_log_rmse_vs_recal_ci95}",
                                recomputed=f"{d:.4f} {s_re}",
                                match=bool(abs(round(d, 4) - lr.d_log_rmse_vs_recal) < 1e-9
                                           and s_re == lr.d_log_rmse_vs_recal_ci95)))
    rep = pd.DataFrame(rep)
    rep.to_csv(OUT / "s7_reproduction_check.csv", index=False)
    print("reproduction:", int(rep.match.sum()), "/", len(rep), "match")

    # ---- R-S7-1 pairwise paired bootstrap (5000, seed 0)
    pw = []
    for sname, sm in splits:
        for cname, cT in CONDS:
            mk = sm & cmask(cT)
            for a, b in PAIRS:
                d, l, h = paired_boot(err[a][mk], err[b][mk], fid[mk], n=5000, seed=0)
                pw.append(dict(split=sname, condition=cname, n_fuels=len(np.unique(fid[mk])), n_obs=int(mk.sum()),
                               model_a=MODELS[a], model_b=MODELS[b], rmse_a=round(rmse(err[a][mk]), 4),
                               rmse_b=round(rmse(err[b][mk]), 4), d_log_rmse_a_minus_b=round(d, 4),
                               ci95_lo=round(l, 4), ci95_hi=round(h, 4), differs=bool(l > 0 or h < 0)))
    pw = pd.DataFrame(pw)
    pw.to_csv(OUT / "s7_pairwise_bootstrap.csv", index=False)

    # ---- R-S7-2 split-conformal factors (dev only)
    fac = []
    ratio_cols = {"GP": ("ratio_gp", "gp"), "ridge": ("abs_res_ridge", "ridge"), "bilinear": ("abs_res_bilinear", "bilinear")}
    factors = {}
    for iname, (col, _) in ratio_cols.items():
        for cname, cT in TCOND:
            s, nf = dev_factor(D, col, cT)
            assert nf >= MIN_DEV_FUELS
            factors[(iname, cT)] = s
            fac.append(dict(interval=iname, scope=f"per condition: {cname}", n_dev_fuels=nf,
                            rank_k=conformal_rank(nf), factor_s=round(s, 4)))
        s, nf = dev_factor(D, col, None)
        fac.append(dict(interval=iname, scope="global over all T (sensitivity, not the declared recipe)",
                        n_dev_fuels=nf, rank_k=conformal_rank(nf), factor_s=round(s, 4)))
    fac = pd.DataFrame(fac)
    fac.to_csv(OUT / "s7_recalibration_factors.csv", index=False)

    # ---- coverage / widths before and after
    cw = []

    def add(iname, kind, centre, half, sname, sm):
        for cname, cT in CONDS:
            mk = sm & cmask(cT)
            ins = np.abs(centre[mk] - y[mk]) <= half[mk]
            cw.append(dict(interval=iname, kind=kind, split=sname, condition=cname,
                           n_fuels=len(np.unique(fid[mk])), n_obs=int(mk.sum()),
                           coverage=round(float(ins.mean()), 3), mean_width_log=round(float(np.mean(2 * half[mk])), 3)))

    sd = D.gp_sd.to_numpy()
    for sname, sm in splits:
        add("GP", "before (+-1.645 sd, v4)", mu, Z90 * sd, sname, sm)
        for iname, (col, key) in ratio_cols.items():
            s_row = np.array([factors[(iname, t)] for t in D.T_C])
            centre = mu if key == "gp" else y + err[key]
            half = s_row * (sd if key == "gp" else 1.0)
            add(iname, "after (split-conformal, per condition, dev-derived s)", centre, half, sname, sm)
    cw = pd.DataFrame(cw)
    cw.to_csv(OUT / "s7_coverage_widths.csv", index=False)

    # ---- R-S7-3 verdict (coverage = overall, all T)
    g = cw[(cw.interval == "GP") & cw.kind.str.startswith("after") & (cw.condition == "overall (all T)")]
    dev_c = float(g[g.split.str.startswith("dev")].coverage.iloc[0])
    lock_c = float(g[g.split.str.startswith("locked")].coverage.iloc[0])
    pd.DataFrame([dict(dev_cov_overall=dev_c, locked_cov_overall=lock_c,
                       adopted=bool(dev_c >= 0.88 and lock_c >= 0.80))]).to_csv(OUT / "s7_verdict_rs73.csv", index=False)

    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(pw.drop(columns=["n_fuels", "n_obs"]).to_string(index=False))
    print(fac.to_string(index=False))
    print(cw.to_string(index=False))
    print("dev", dev_c, "locked", lock_c, "adopted", dev_c >= 0.88 and lock_c >= 0.80)
    print(rep[~rep.match].to_string(index=False))


if __name__ == "__main__":
    main()
