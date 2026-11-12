"""v5 BASE Task 8: carried cells - CRC 647 r3 (dev 45 / fold_b 12, recomputed from the v4 per-fuel files), FACE diesel (S6),
designs (D4), octane (v4 O), smoke (v4 SP), plus the explicit no_evaluation cells of the CRC 647, T3 and diesel columns.

Nothing is refitted except leave-one-group-out constants (the trivial model).  fold_b = the former CRC 647 locked 12, read five
times; its trivial constant is the dev mean.  Output: results/v5/BASE/cells_carried.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row, no_eval  # noqa: E402
from fuelblend.eval.group_stats import logo_constant, verdict  # noqa: E402

R = ROOT / "results"
OUT = R / "v5" / "BASE"
CRC = "CRC647_r3"
S1_MARGIN = {"T_lt_263K": 0.4154, "T_263_283K": 0.6084}   # results/v5/S1/uncertainty_v5.yaml (unsafe side, cold viscosity)
OVER5 = lambda e: e > 5.0                                  # noqa: E731
MARG = yaml.safe_load(open(ROOT / "configs/uncertainty_v4.yaml"))
MARG = MARG.get("margins", MARG)


def margin(key: str) -> float:
    m = MARG.get(key)
    return float(m["value"]) if m else np.nan


def trivial(y: np.ndarray, groups: np.ndarray, role: np.ndarray, strata=None) -> np.ndarray:
    """dev rows: leave-one-group-out constant over dev; fold_b rows: the dev mean (within strata when given)."""
    y, role = np.asarray(y, float), np.asarray(role)
    st = np.zeros(len(y), dtype=object) if strata is None else np.asarray(strata).astype(object)
    out = np.full(len(y), np.nan)
    dev = role == "dev"
    out[dev] = logo_constant(y[dev], np.asarray(groups)[dev], st[dev])
    for s in np.unique(st[~dev]):
        m, tr = (~dev) & (st == s), dev & (st == s)
        out[m] = y[tr].mean() if tr.any() else y[dev].mean()
    return out


def coverage(err, groups, m: float) -> float:
    """Fraction of groups whose max |err| is inside the margin (the v4 per-fuel conformal score)."""
    if not np.isfinite(m):
        return np.nan
    s = pd.Series(np.abs(np.asarray(err, float))).groupby(np.asarray(groups)).max()
    return float((s <= m).mean())


def crc_rows(prop, df, model, variant, err, meas, ev, *, T=None, strata=None, unsafe=None, mkey=None, ln=False, note=""):
    """dev and fold_b rows for one model.  df needs fuel_id and role ('dev'|'fold_b')."""
    triv_pred = trivial(np.log(meas) if ln else meas, df.fuel_id.to_numpy(), df.role.to_numpy(), strata)
    et = triv_pred - (np.log(meas) if ln else meas)
    out = []
    for lab in ("dev", "fold_b"):
        k = (df.role == lab).to_numpy()
        if not k.any():
            continue
        m = margin(mkey) if mkey else np.nan
        r = make_row(prop, CRC, "crc647_r3", model, variant, "carried_v4", lab, ev, err=np.asarray(err)[k], groups=df.fuel_id.to_numpy()[k],
                     unsafe=unsafe, err_trivial=et[k], trivial_model="logo_constant" if lab == "dev" else "dev_mean_constant",
                     margin90=m, coverage90=coverage(np.asarray(err)[k], df.fuel_id.to_numpy()[k], m),
                     extra=dict(T_K=np.nan if T is None else T, note=note))
        out.append(r)
    return out


def role_of(split_role: pd.Series) -> pd.Series:
    return np.where(split_role == "locked", "fold_b", "dev")


# ------------------------------------------------------------------------------------------------ CRC 647 r3
def crc_flash() -> list[dict]:
    f = pd.read_csv(R / "v4/R3/F/flash_real_fuels.csv")
    f = f[f.family == "crc647"].assign(role=lambda x: role_of(x.split_role))
    ev = "results/v4/R3/F/flash_real_fuels.csv"
    out = crc_rows("flash_point_C", f, "ideal_le_chatelier", "designated", f.ideal_C - f.measured_C, f.measured_C.to_numpy(), ev,
                   strata=f.method.to_numpy(), unsafe=OVER5)
    out += crc_rows("flash_point_C", f, "unifac_le_chatelier", "unifac", f.unifac_C - f.measured_C, f.measured_C.to_numpy(), ev,
                    strata=f.method.to_numpy(), unsafe=OVER5)
    for r in out:
        r["trivial_model"] = r["trivial_model"].replace("logo_constant", "logo_constant_within_method")
    # LFL overlay (L1): only aggregates exist (no per-fuel file) -> fixed numbers, no verdict
    l1 = pd.read_csv(R / "v5/L1/comparison.csv")
    for lab, sl in (("dev", "D93 dev"), ("fold_b", "D93 locked")):
        for mdl, var in (("ideal", "lfl_overlay"), ("unifac", "unifac_lfl_overlay")):
            q = l1[(l1.scope == "real_fuels_CRC647") & (l1.slice == sl) & (l1.model == mdl)].iloc[0]
            out.append(make_row("flash_point_C", CRC, "crc647_r3", f"{mdl}_le_chatelier", var, "aggregated_v5", lab, "results/v5/L1/comparison.csv",
                                fixed=dict(n_groups=int(q.N), n_rows=int(q.N), bias=q.bias_C_after, mae=q.MAE_C_after, unsafe_rate=q.over_gt5C_rate_after),
                                extra=dict(T_K=np.nan, note="L1 aggregate only; no per-fuel file, so no paired verdict")))
    out.append(no_eval("flash_point_C", CRC, "crc647_r3", "mesitylene fix not run on CRC 647 (Part C Task 10)", model="ideal_le_chatelier",
                       variant="mesitylene_fix", evidence="results/v5/BASE/cells_t2_flash.csv", label="dev"))
    return out


def crc_viscosity() -> list[dict]:
    v = pd.read_csv(R / "v4/R3/V/viscosity_predictions_v4.csv")
    v = v.assign(role=np.where(v.locked, "fold_b", "dev"), T_K=v.T_K.round(2))
    ev = "results/v4/R3/V/viscosity_predictions_v4.csv"
    models = [c[5:] for c in v.columns if c.startswith("pred_")]
    out = []
    for m in models:
        var = "designated" if m == "Recalibrated Refutas" else "carried"
        mk_all = "viscosity_kin_cSt.crc647.recalibrated_refutas.all_T" if m == "Recalibrated Refutas" else None
        e_all = np.log(v[f"pred_{m}"]) - np.log(v.y)
        # overall (all three temperatures): trivial = per-temperature constant; unsafe from the per-row temperature
        uns = np.where(v.T_K < 263.15, e_all < -S1_MARGIN["T_lt_263K"], np.where(v.T_K < 283.15, e_all < -S1_MARGIN["T_263_283K"], False))
        rows = crc_rows("viscosity_kin_cSt", v, m, var, e_all, v.y.to_numpy(), ev, strata=v.T_K.to_numpy(), mkey=mk_all, ln=True,
                        note="ln(pred/meas); all of -40/-20/+20 C; unsafe = under-prediction beyond the S1 margin")
        for r in rows:
            k = (v.role == r["label"]).to_numpy()
            r["unsafe_rate"] = float(np.mean(uns[k]))
        out += rows
        for T, g in v.groupby("T_K"):
            mk = mk_all
            if m == "Recalibrated Refutas" and abs(T - 233.15) < 0.1:
                mk = "viscosity_kin_cSt.crc647.recalibrated_refutas.minus40C"
            un = None if T >= 283.15 else (lambda e, mm=S1_MARGIN["T_lt_263K" if T < 263.15 else "T_263_283K"]: e < -mm)
            out += crc_rows("viscosity_kin_cSt", g, m, var, np.log(g[f"pred_{m}"]) - np.log(g.y), g.y.to_numpy(), ev, T=T, unsafe=un,
                            mkey=mk, ln=True, note="ln(pred/meas)")
    return out


def crc_bulk() -> list[dict]:
    out = []
    d = pd.read_csv(R / "v4/R3/B/density_v4.csv").assign(role=lambda x: role_of(x.split_role))
    ev = "results/v4/R3/B/density_v4.csv"
    out += crc_rows("density_kg_m3", d, "ideal_volume", "designated", d.err_kg_m3, d.measured_kg_m3.to_numpy(), ev, strata=d.T_C.round(2).to_numpy(),
                    mkey="density_kg_m3.crc647_ideal_volume.all_T", note="all temperatures; trivial per temperature")
    for T, g in d.groupby(d.T_C.round(2)):
        out += crc_rows("density_kg_m3", g, "ideal_volume", "designated", g.err_kg_m3, g.measured_kg_m3.to_numpy(), ev, T=round(T + 273.15, 2),
                        mkey="density_kg_m3.crc647_ideal_volume.T15_56C" if abs(T - 15.56) < 0.01 else None)
    h = pd.read_csv(R / "v4/R3/B/hoc_v4.csv").assign(role=lambda x: role_of(x.split_role))
    ev = "results/v4/R3/B/hoc_v4.csv"
    out += crc_rows("heat_of_combustion_net_MJ_kg", h, "ga_comp_H", "designated", h.err_comp_H, h.measured_ncv_MJ_kg.to_numpy(), ev,
                    mkey="heat_of_combustion_net_MJ_kg.crc647.group_additivity_composition_H")
    out += crc_rows("heat_of_combustion_net_MJ_kg", h, "ga_const_H", "const_H", h.err_v2_const_H, h.measured_ncv_MJ_kg.to_numpy(), ev,
                    mkey="heat_of_combustion_net_MJ_kg.crc647.group_additivity_v2_constH13.8")
    z = pd.read_csv(R / "v4/R3/B/freezing_v4.csv").assign(role=lambda x: role_of(x.split_role))
    ev = "results/v4/R3/B/freezing_v4.csv"
    for col, mdl, var in (("err_ideal_sle_mole", "ideal_svl", "designated"), ("err_wat_ideal", "wat_ideal_library", "carried"),
                          ("err_wat_minus_devfit_offset", "wat_minus_devfit_offset", "carried"), ("err_coutinho", "coutinho_library", "carried")):
        out += crc_rows("freezing_point_C", z, mdl, var, z[col], z.measured_C.to_numpy(), ev)
    for p, why in (("smoke_point_mm", "CRC 647 r3 has no smoke-point evaluation in v4"),
                   ("distillation_T10_C", "distillation is an input of the CRC 647 reconstruction, not a target")):
        out.append(no_eval(p, CRC, "crc647_r3", why, label="dev"))
    out.append(no_eval("dcn", CRC, "crc647_r3", "no cetane data for CRC 647 fuels; v4 cetane model uncertified", label="dev"))
    return out


# ------------------------------------------------------------------------------------------------ diesel (S6)
def diesel() -> list[dict]:
    s = pd.read_csv(R / "v5/S6/summary_R-S6-2.csv")
    s = s[s.primary]
    ev = "results/v5/S6/summary_R-S6-2.csv"
    out = []
    for r in s.itertuples():
        if r.n == 0:
            row = no_eval(r.property, "diesel", "face_diesel", f"S6: no verified pure-component {r.property.split('_')[0]} data ({r.lab}); not evaluable",
                          model=r.model, evidence=ev)
            row["verdict"] = "n_too_small"
        else:
            row = make_row(r.property, "diesel", "face_diesel", r.model, "designated", "aggregated_v5", "descriptive", ev,
                           fixed=dict(n_groups=int(r.n), n_rows=int(r.n), bias=r.bias, mae=r.mae,
                                      unsafe_rate=r.over_gt5C_rate if r.property.startswith("flash") else np.nan),
                           extra=dict(note=f"lab {r.lab}; n = {r.n} FACE fuels"))
            row["verdict"] = "n_too_small"
        out.append(row)
    sv = pd.read_csv(R / "v5/S6/surrogate_viscosity_R-S6-3.csv")
    ok = sv[sv.status == "evaluated"]
    out.append(make_row("viscosity_kin_cSt", "diesel", "face_surrogates", "refutas", "designated", "aggregated_v5", "descriptive",
                        "results/v5/S6/surrogate_viscosity_R-S6-3.csv", err=ok.ln_err, groups=ok.fuel_id,
                        n_excluded_nonfinite=int((sv.status != "evaluated").sum()), extra=dict(note="ln(pred/meas); surrogate blends, not real diesel")))
    for p, why in (("freezing_point_C", "diesel has cloud point only, no freezing point"), ("heat_of_combustion_net_MJ_kg", "no heat measurements on FACE diesel"),
                   ("smoke_point_mm", "no smoke-point measurements on FACE diesel"), ("distillation_T10_C", "no distillation model scored on diesel"),
                   ("pour_point_C", "no pour-point model scored on diesel")):
        out.append(no_eval(p, "diesel", "-", why))
    return out


# ------------------------------------------------------------------------------------------------ designs, octane, smoke
def designs() -> list[dict]:
    d = pd.read_csv(R / "v5/D4/designs_v4.csv")
    return [make_row("designs_certified", "designs", "d4", "v4_designer", "designated", "aggregated_v5", "descriptive", "results/v5/D4",
                     fixed=dict(n_rows=len(d), bias=int(d.passed.sum())), extra=dict(note="bias = number of designs certified; " +
                                                                                     "; ".join(f"{a}: {b}" for a, b in zip(d.design_id, d.category))))]


def octane() -> list[dict]:
    o = pd.read_csv(R / "v4/O/octane_eval_v4.csv")
    ev = "results/v4/O/octane_eval_v4.csv"
    met, pdiff = o[o.record_type == "metric"], o[(o.record_type == "paired_diff") & (o.ref_model == "constant_mean")]
    out = []
    for r in met.itertuples():
        q = pdiff[(pdiff.view == r.view) & (pdiff.prop == r.prop) & (pdiff.model == r.model)]
        fx = dict(n_groups=int(r.n_groups), n_rows=int(r.n), bias=r.bias, mae=r.mae, rmse=r.rmse)
        if r.model != "constant_mean" and len(q):
            q = q.iloc[0]
            lo, hi = -q.ci_hi, -q.ci_lo
            fx.update(diff_trivial_minus_model=-q.mae, ci_lo=lo, ci_hi=hi, verdict=verdict(int(q.n_groups), lo, hi))
        out.append(make_row(r.prop, "T1_lab", f"octane_{r.view}", r.model, "designated" if r.model == "base_unchanged" else "carried", "carried_v4",
                            "already_read" if r.view == "locked_thin_test" else "dev", ev, trivial_model="constant_mean", fixed=fx))
    return out


def smoke() -> list[dict]:
    s = pd.read_csv(R / "v4/SP/smoke_point_eval.csv")
    ev = "results/v4/SP/smoke_point_eval.csv"
    out = []
    for r in s[s.slice.str.startswith("ajftd")].itertuples():
        fx = dict(n_groups=np.nan, n_rows=int(r.n), bias=r.bias_mm, mae=r.mae_mm, rmse=r.rmse_mm)
        if np.isfinite(r.mae_diff_vs_const_ci_lo):
            lo, hi = -r.mae_diff_vs_const_ci_hi, -r.mae_diff_vs_const_ci_lo
            fx.update(ci_lo=lo, ci_hi=hi, verdict=verdict(int(r.n), lo, hi))
        out.append(make_row("smoke_point_mm", "T3_ajftd", r.slice, r.model, "designated" if r.model == "arom+diarom_fixed" else "carried", "carried_v4",
                            "dev", ev, trivial_model="const_mean", fixed=fx, extra=dict(note="v4 SP; group = fuel unless slice says otherwise; n_groups not stored")))
    out.append(no_eval("smoke_point_mm", "T3_ajftd", "ysi_route_real_fuels", "v4 SP: YSI route had 0 evaluable real fuels", model="ysi_mix+smoke_point_from_ysi", evidence=ev))
    for p, why in (("density_kg_m3", "density is an input of the T3 bulk models, not a target"),
                   ("distillation_T10_C", "distillation is an input of the T3 bulk models, not a target"),
                   ("pour_point_C", "no pour-point measurements on AJFTD fuels"), ("dcn", "no cetane data on AJFTD fuels")):
        out.append(no_eval(p, "T3_ajftd", "-", why))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = crc_flash() + crc_viscosity() + crc_bulk() + diesel() + designs() + octane() + smoke()
    df = frame(rows)
    df.to_csv(OUT / "cells_carried.csv", index=False, float_format="%.8g")
    print(df.groupby(["tier", "property", "status"]).size().to_string())


if __name__ == "__main__":
    main()
