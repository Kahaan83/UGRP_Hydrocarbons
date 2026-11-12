"""v5 BASE Task 7: T3 AJFTD bulk-descriptor cells (aggregated from S4/S4b, never re-run) and the CRC 647 pour-point baseline.

Locked S4 rows are aggregated from the committed locked_scores.csv (label already_read); `cmd_locked` is never called.
Outputs: results/v5/BASE/cells_t3_bulk.csv, rows_pour_point.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row, no_eval  # noqa: E402
from fuelblend.eval.group_stats import logo_constant, verdict  # noqa: E402

R = ROOT / "results"
OUT = R / "v5" / "BASE"
TIER = "T3_ajftd"
PROP = {"flash": "flash_point_C", "freezing": "freezing_point_C", "heat": "heat_of_combustion_net_MJ_kg", "visc_m20": "viscosity_kin_cSt"}
UNVERIFIED = {"api_2b7_t10", "butler_t10"}                 # flash correlations whose coefficients are unverified (audit B4)
BASELINE_VARIANT = "designated"


def s4_row(r, label: str, ev: str) -> dict:
    """One S4 score row -> scoreboard row.  S4 stores mae - mae_const (positive = worse); the scoreboard stores trivial - model."""
    prop = PROP[r.target]
    var = "unverified_coefficients" if r.model in UNVERIFIED else BASELINE_VARIANT
    d, lo, hi = -r.mae_diff_vs_const, -r.ci_hi, -r.ci_lo
    trivial = "logo_constant"
    row = make_row(prop, TIER, f"ajftd_bulk_{r.scheme}_{r.subset}", r.model, var, "aggregated_v5", label, ev, trivial_model=trivial,
                   margin90=r.margin90, coverage90=r.coverage90,
                   fixed=dict(n_groups=int(r.n_groups), n_rows=int(r.n_fuels), bias=r.bias_pred_minus_meas, mae=r.mae, mae_trivial=r.mae_const,
                              diff_trivial_minus_model=d, ci_lo=lo, ci_hi=hi, verdict="not_applicable" if r.model == "const" else verdict(int(r.n_groups), lo, hi)))
    return row


def s4_rows() -> list[dict]:
    rows = []
    d = pd.read_csv(R / "v5/S4/dev_scores.csv")
    for r in d[(d.scheme == "logo") & (d.subset == "all")].itertuples():
        rows.append(s4_row(r, "dev", "results/v5/S4/dev_scores.csv"))
    lk = pd.read_csv(R / "v5/S4/locked_scores.csv")
    for r in lk[lk.subset == "all"].assign(scheme="locked").itertuples():
        rows.append(s4_row(r, "already_read", "results/v5/S4/locked_scores.csv"))
    return rows


def flash_method_stratified() -> list[dict]:
    p = pd.read_csv(R / "v5/S4/dev_predictions.csv", low_memory=False)
    p = p[(p.scheme == "logo") & (p.target == "flash")]
    base = p[p.model == "const"].set_index("fuel_id")
    rows = []
    for m in ("ridge", "hgb"):
        q = p[p.model == m].set_index("fuel_id")
        c = base.loc[q.index]
        const_m = pd.Series(logo_constant(q.y.to_numpy(), q["cell"].to_numpy(), q.flash_family.to_numpy()), index=q.index)
        rows.append(make_row("flash_point_C", TIER, "ajftd_bulk_logo_all", m, "designated", "recomputed_v5", "dev", "results/v5/S4/dev_predictions.csv",
                             err=q.pred - q.y, groups=q["cell"], unsafe=lambda e: e > 5.0, err_trivial=const_m - q.y,
                             trivial_model="logo_constant_within_method"))
    return rows


def heat_measured_rows() -> list[dict]:
    m = pd.read_csv(R / "v5/S4b/s4b_metrics.csv")
    b = pd.read_csv(R / "v5/S4b/s4b_bootstrap.csv")
    rows = []
    for r in m.itertuples():
        lab = "already_read" if "locked" in r.reading else "dev"
        bb = b[(b.reading == r.reading) & (b.model == r.model)]
        fx = dict(n_groups=int(r.n_groups), n_rows=int(r.n_fuels), mae=r.mae_model, mae_trivial=r.mae_C_within)
        if len(bb):
            x = bb.iloc[0]
            fx.update(diff_trivial_minus_model=x.diff_Cwithin_minus_model, ci_lo=x.ci_lo, ci_hi=x.ci_hi, verdict=verdict(int(r.n_groups), x.ci_lo, x.ci_hi))
        elif r.model != "C_all":
            fx.update(diff_trivial_minus_model=r.diff_Cwithin_minus_model, ci_lo=r.ci_lo, ci_hi=r.ci_hi, verdict=verdict(int(r.n_groups), r.ci_lo, r.ci_hi))
        rows.append(make_row("heat_of_combustion_net_MJ_kg", TIER, r.reading, r.model, "designated", "aggregated_v5", lab,
                             "results/v5/S4b/s4b_metrics.csv; results/v5/S4b/s4b_bootstrap.csv", trivial_model="logo_constant_within_method", fixed=fx,
                             extra=dict(note=str(r.note) if isinstance(r.note, str) else "")))
    return rows


# ------------------------------------------------------------------------------------------------------------- pour point (CRC 647, n = 57)
def pour_rows() -> tuple[list[dict], pd.DataFrame]:
    j = pd.read_csv(ROOT / "data/real_fuels/v4/jet_properties_v4.csv")
    w = j.pivot_table(index="fuel_id", columns="property", values="value", aggfunc="median")
    dens = j[(j.property == "density_kg_m3") & (j.T_C.round(1) == 15.6)].groupby("fuel_id").value.median()
    w["density15"] = dens
    cols = ["distillation_T10_C", "distillation_T50_C", "distillation_T90_C", "aromatics_vol_pct", "density15"]
    d = w.dropna(subset=["pour_point_C", "freezing_point_C"] + cols).copy()
    n = len(d)
    y = d.pour_point_C.to_numpy()
    X = d[cols].to_numpy(float)
    idx = np.arange(n)
    grp = d.index.to_numpy()
    pred_const, pred_ridge, pred_off = np.zeros(n), np.zeros(n), np.zeros(n)
    gap = (d.freezing_point_C - d.pour_point_C).to_numpy()

    def ridge_fit(Xt, yt, a):
        mu, sd = Xt.mean(0), Xt.std(0) + 1e-9
        Z = (Xt - mu) / sd
        beta = np.linalg.solve(Z.T @ Z + a * np.eye(Z.shape[1]), Z.T @ (yt - yt.mean()))
        return lambda Xn: yt.mean() + ((Xn - mu) / sd) @ beta

    for i in idx:
        tr = idx != i
        pred_const[i] = y[tr].mean()
        best_a, best = None, np.inf
        for a in (0.1, 1.0, 10.0):                                    # inner leave-one-out
            e = 0.0
            for k in idx[tr]:
                tr2 = tr & (idx != k)
                e += (ridge_fit(X[tr2], y[tr2], a)(X[[k]])[0] - y[k]) ** 2
            if e < best:
                best, best_a = e, a
        pred_ridge[i] = ridge_fit(X[tr], y[tr], best_a)(X[[i]])[0]
        pred_off[i] = d.freezing_point_C.iloc[i] - gap[tr].mean()
    out = d.reset_index()[["fuel_id", "pour_point_C", "freezing_point_C"]].assign(pred_const=pred_const, pred_ridge=pred_ridge, pred_freezing_offset=pred_off)
    ev = "results/v5/BASE/rows_pour_point.csv"
    et = pred_const - y
    rows = []
    for m, p, tm, var, nt in (("const", pred_const, "", "designated", ""), ("ridge", pred_ridge, "logo_constant", "designated", ""),
                              ("freezing_offset", pred_off, "logo_constant", "uses_measured_freezing", "uses another measured property")):
        r = make_row("pour_point_C", "CRC647_r3", "crc647_measured_bulk_n57", m, var, "recomputed_v5", "dev", ev, err=p - y, groups=grp,
                     err_trivial=et if m != "const" else None, trivial_model=tm,
                     extra=dict(note="CRC 647 only; features: D86 T10/T50/T90, aromatics vol %, density 15.6 C; leave-one-fuel-out, inner LOO alpha in [0.1,1,10]"
                                + ("; " + nt if nt else "")))
        if m == "freezing_offset":
            r["verdict"] = "not_applicable"
        rows.append(r)
    return rows, out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pr, pout = pour_rows()
    pout.to_csv(OUT / "rows_pour_point.csv", index=False, float_format="%.6g")
    rows = s4_rows() + flash_method_stratified() + heat_measured_rows() + pr + [
        no_eval("density_kg_m3", TIER, "-", "density is an input to the bulk models, not a target"),
        no_eval("smoke_point_mm", TIER, "-", "see Task 8 carried v4 smoke-point cell"),
    ]
    df = frame(rows)
    df.to_csv(OUT / "cells_t3_bulk.csv", index=False, float_format="%.8g")
    print(df[df.status != "no_evaluation"].groupby(["property", "label"]).size().to_string())


if __name__ == "__main__":
    main()
