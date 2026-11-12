"""v5 BASE Task 14: freezing follow-ups (F1x variant, wax-appearance T1 source, D86 lead).

Outputs (results/v5/BASE/delta_freezing/): delta.csv, rows_wax_appearance.csv, rows_d86_ridge.csv, readings.json
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402
from fuelblend.eval.group_stats import logo_constant  # noqa: E402
from fuelblend.mixture.t2_residual import fit_predict_residual  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_freezing"
V5 = ROOT / "data" / "real_fuels" / "v5"
R_GAS = 8.314462618
TIERS = {"two_compilations": ("v4_verified", "two_compilations"), "single_compilation": ("v4_verified", "two_compilations", "single_compilation")}


def f1x_rows() -> tuple[list[dict], dict]:
    s = pd.read_csv(ROOT / "results/v5/S5/s5_predictions.csv")
    s = s[(s.role == "headline") & s.F1x.notna() & s.freezing_C.notna()].copy()
    s["Cl"] = logo_constant(s.freezing_C, s.group)                       # pooled LOGO constant (S5 'C' column is per-set)
    rows, mae = [], {}
    for name, g in [("ajftd_gcxgc", s[s.source == "AJFTD"]), ("av3322", s[s.source == "AV-33-22"]), ("pooled", s)]:
        triv = g.C if name != "pooled" else g.Cl
        r = make_row("freezing_point_C", "T2_measured", name, "F1x_coutinho_complete", "F1x", "recomputed_v5", "descriptive",
                     "results/v5/S5/s5_predictions.csv", err=g.F1x - g.freezing_C, groups=g.group, err_trivial=triv - g.freezing_C,
                     trivial_model="logo_constant")
        raw = make_row("freezing_point_C", "T2_measured", name, "F1", "", "recomputed_v5", "descriptive", "", err=g.F1 - g.freezing_C, groups=g.group)
        r.update(baseline_mae=raw["mae"], baseline_bias=raw["bias"], baseline_unsafe_rate=np.nan)
        rows.append(r)
        mae[name] = r
    return rows, mae


def ideal_sle_K(x, tm, dh) -> float:
    inv = 1.0 / tm - R_GAS / dh * np.log(x)
    return float(np.max(1.0 / inv))


def wax_rows() -> tuple[list[dict], pd.DataFrame]:
    w = pd.read_csv(ROOT / "data/processed/v2/wax_dauphin_coutinho.csv")
    c = pd.read_csv(ROOT / "results/v5/SLE2/constants_extended.csv")
    rows, recs = [], []
    for tier, keep in TIERS.items():
        k = c[c.tier.isin(keep) & c.tm_K.notna() & c.dh_J_mol.notna()].set_index("inchikey")
        for r in w.itertuples():
            iks, x = json.loads(r.inchikeys), np.array(json.loads(r.mole_fractions), float)
            ok = all(i in k.index for i in iks)
            pred = ideal_sle_K(x, k.loc[iks, "tm_K"].to_numpy(float), k.loc[iks, "dh_J_mol"].to_numpy(float)) if ok else np.nan
            recs.append(dict(tier=tier, system_id=r.system_id, system_name=r.system_name, n_comp=r.n_comp, measured_K=r.wat_measured_K,
                             pred_K=pred, method="wax_appearance", status="ok" if ok else "missing_constants"))
        d = pd.DataFrame([q for q in recs if q["tier"] == tier])
        s = d[d.status == "ok"]
        if s.empty:
            rows.append(dict(make_row("freezing_point_C", "T1_lab", "wax_appearance_dauphin_coutinho", "ideal_svl", f"sle2_{tier}", "no_evaluation",
                                      "exploratory", ""), reason=f"no system has constants for all components ({len(d)} systems)"))
            continue
        triv = logo_constant(s.measured_K, s.system_id) if s.system_id.nunique() > 1 else np.full(len(s), np.nan)
        r = make_row("freezing_point_C", "T1_lab", "wax_appearance_dauphin_coutinho", "ideal_svl", f"sle2_{tier}", "recomputed_v5", "exploratory",
                     "data/processed/v2/wax_dauphin_coutinho.csv; results/v5/SLE2/constants_extended.csv", err=s.pred_K - s.measured_K,
                     groups=s.system_id, err_trivial=triv - s.measured_K, trivial_model="logo_constant",
                     extra=dict(method_label="wax_appearance", n_systems_missing_constants=int((d.status != "ok").sum())))
        rows.append(r)
    return rows, pd.DataFrame(recs)


def d86_rows() -> tuple[list[dict], pd.DataFrame]:
    p = pd.read_csv(V5 / "crc_av3322_properties.csv")
    p = p[p.qualifier.isna() & p.value.notna()]
    piv = p[p.property.isin(["freezing_point_C", "distillation_T90_C", "distillation_T100_C"])].groupby(["fuel_id", "property"]).value.median().unstack()
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    npar = comp[comp.group == "n-alkane"].groupby("fuel_id").mass_pct.sum() / comp.groupby("fuel_id").mass_pct.sum()
    piv["n_paraffin_share"] = npar
    fu = pd.read_csv(ROOT / "results/v5/BASE/t2_fuels.csv")              # S3 double-count groups (the composition file's column is one constant)
    grp = fu[fu.source == "av3322"].set_index("fuel_id").group
    piv["group"] = grp.reindex(piv.index).fillna(pd.Series(piv.index, index=piv.index))
    d = piv.dropna(subset=["freezing_point_C", "distillation_T90_C", "distillation_T100_C", "n_paraffin_share"]).reset_index()
    X = d[["distillation_T90_C", "distillation_T100_C", "n_paraffin_share"]].to_numpy(float)
    y, g = d.freezing_point_C.to_numpy(float), d.group.astype(str).to_numpy()
    pred = np.full(len(d), np.nan)
    for gg in np.unique(g):
        te = g == gg
        pred[te] = fit_predict_residual(X[~te], y[~te], g[~te], X[te])
    d["pred_ridge"], d["pred_const"] = pred, logo_constant(y, g)
    r = make_row("freezing_point_C", "T2_measured", "av3322", "ridge_d86_T90_FBP_nparaffin", "d86_lead", "recomputed_v5", "exploratory",
                 "data/real_fuels/v5/crc_av3322_properties.csv; crc_av3322_composition.csv", err=d.pred_ridge - y, groups=g,
                 err_trivial=d.pred_const - y, trivial_model="logo_constant")
    return [r], d


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    a, am = f1x_rows()
    b, wax = wax_rows()
    c, dd = d86_rows()
    frame(a + b + c).to_csv(OUT / "delta.csv", index=False)
    wax.to_csv(OUT / "rows_wax_appearance.csv", index=False, float_format="%.4f")
    dd.to_csv(OUT / "rows_d86_ridge.csv", index=False, float_format="%.4f")
    res = dict(f1x_pooled_mae=am["pooled"]["mae"], f1x_pooled_n=am["pooled"]["n_rows"], reproduces_3p07=bool(abs(am["pooled"]["mae"] - 3.07) <= 0.02),
               f1x_verdict={k: v["verdict"] for k, v in am.items()}, ajftd_groups=am["ajftd_gcxgc"]["n_groups"],
               reading_pass=bool(am["av3322"]["verdict"] == "beats_trivial" and am["pooled"]["verdict"] == "beats_trivial" and am["ajftd_gcxgc"]["n_groups"] < 20),
               wax={r["variant"]: dict(n=r["n_groups"], mae=r["mae"], bias=r["bias"], verdict=r["verdict"]) for r in b},
               d86=dict(n=c[0]["n_groups"], mae=c[0]["mae"], mae_const=c[0]["mae_trivial"], verdict=c[0]["verdict"], ci=[c[0]["ci_lo"], c[0]["ci_hi"]]))
    (OUT / "readings.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
