#!/usr/bin/env python3
"""eval_octane_blending_benchmarks.py - Octane blending benchmark on KAUST Table 4 (and Table 5 FACE).

Models (RON, MON):
  - Constant: training-fold mean of the measured value (no-skill reference)
  - Linear volume and linear molar blending, with pure RON/MON taken from measured single-component
    rows (see kaust_registry.pure_props), not typed literals
  - Linear molar + GBM residual, evaluated with GroupKFold over base-fuel families

Groups: the component set with ethanol removed, so a hydrocarbon base fuel and its ethanol blends
(TMB+0 ... TMB+60, PRF + ethanol, TPRF + ethanol) always fall in the same fold.
Single-component anchor rows (TMB+0, CP+0, Hex+0) define the pure values, so they are excluded from
the error metrics. Report text contains computed values only.
Does not touch data/unified/.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold

from _common import OUT_DIR, STAGING_DIR, mol_info
from kaust_registry import COMPONENTS, densities, load_tables, pure_props

NAMES = list(COMPONENTS)
N_SPLITS = 5


def build_table(t3, t4):
    pp = pure_props(t3, t4)
    dens = densities()
    mw = {k: mol_info(v["smiles"])[1] for k, v in COMPONENTS.items()}
    recs = []
    for r in t4.itertuples():
        comp = json.loads(r.composition_vol_pct)
        tot = sum(comp.values())
        v = {k: x / tot for k, x in comp.items()}
        n = {k: v[k] * dens[k][0] / mw[k] for k in v}
        nt = sum(n.values())
        x = {k: n[k] / nt for k in n}
        rec = {"blend_name": r.blend_name, "n_comp": len(comp),
               "group": "+".join(sorted(c for c in comp if c != "ethanol")) or "ethanol",
               "ethanol_vol_frac": v.get("ethanol", 0.0),
               "ron_meas": float(r.ron), "mon_meas": float(r.mon)}
        for t in ("ron", "mon"):
            rec[f"{t}_lin_vol"] = sum(v[k] * pp[k][t] for k in v)
            rec[f"{t}_lin_mol"] = sum(x[k] * pp[k][t] for k in x)
        for k in NAMES:
            rec[f"v_{k}"] = v.get(k, 0.0)
        recs.append(rec)
    return pd.DataFrame(recs), pp, dens


def metrics(err):
    err = np.asarray(err, dtype=float)
    return {"mae": float(np.mean(np.abs(err))), "rmse": float(np.sqrt(np.mean(err ** 2))),
            "bias": float(np.mean(err)), "max_err": float(np.max(np.abs(err)))}


def cv_predictions(df, target):
    X = df[[f"v_{k}" for k in NAMES]].values
    y = df[f"{target}_meas"].values
    base = df[f"{target}_lin_mol"].values
    groups = df["group"].values
    n_splits = min(N_SPLITS, len(set(groups)))
    const = np.zeros(len(df))
    gbm = np.zeros(len(df))
    for tr, te in GroupKFold(n_splits=n_splits).split(X, y, groups):
        const[te] = y[tr].mean()
        m = GradientBoostingRegressor(n_estimators=45, max_depth=3, learning_rate=0.08, random_state=42)
        m.fit(X[tr], y[tr] - base[tr])
        gbm[te] = base[te] + m.predict(X[te])
    return const, gbm, n_splits


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t3, t4, t5 = load_tables()
    full, pp, dens = build_table(t3, t4)
    anchors = full[full["n_comp"] == 1]
    df = full[full["n_comp"] > 1].reset_index(drop=True)

    models = {}
    n_splits = None
    for t in ("ron", "mon"):
        const, gbm, n_splits = cv_predictions(df, t)
        df[f"{t}_pred_const"] = const
        df[f"{t}_pred_gbm"] = gbm
        models[t] = {
            "Constant (train mean)": df[f"{t}_pred_const"],
            "Linear volume": df[f"{t}_lin_vol"],
            "Linear molar": df[f"{t}_lin_mol"],
            "Linear molar + GBM residual (GroupKFold)": df[f"{t}_pred_gbm"],
        }
    rows = []
    for t, d in models.items():
        for name, pred in d.items():
            rows.append({"target": t.upper(), "model": name, "n": len(df),
                         **{k: round(v, 3) for k, v in metrics(pred - df[f"{t}_meas"]).items()}})
    summary = pd.DataFrame(rows)

    fam_rows = []
    for g, sub in df.groupby("group"):
        r = {"group": g, "n": len(sub)}
        for t in ("ron", "mon"):
            r[f"{t}_mae_lin_mol"] = round(float(np.mean(np.abs(sub[f"{t}_lin_mol"] - sub[f"{t}_meas"]))), 3)
            r[f"{t}_mae_lin_vol"] = round(float(np.mean(np.abs(sub[f"{t}_lin_vol"] - sub[f"{t}_meas"]))), 3)
            r[f"{t}_mae_gbm"] = round(float(np.mean(np.abs(sub[f"{t}_pred_gbm"] - sub[f"{t}_meas"]))), 3)
        fam_rows.append(r)
    fam = pd.DataFrame(fam_rows)

    # Density sensitivity: linear-molar MAE if every density is shifted by +/-5 percent
    sens = []
    for scale in (0.95, 1.05):
        mw = {k: mol_info(v["smiles"])[1] for k, v in COMPONENTS.items()}
        errs = {"ron": [], "mon": []}
        for r in t4[t4["blend_name"].isin(df["blend_name"])].itertuples():
            comp = json.loads(r.composition_vol_pct)
            tot = sum(comp.values())
            n = {k: (x / tot) * dens[k][0] * (scale if k == "diisobutylene" else 1.0) / mw[k] for k, x in comp.items()}
            nt = sum(n.values())
            for t in errs:
                pred = sum(n[k] / nt * pp[k][t] for k in n)
                errs[t].append(pred - float(getattr(r, t)))
        sens.append({"dib_density_scale": scale,
                     "ron_mae_lin_mol": round(float(np.mean(np.abs(errs["ron"]))), 3),
                     "mon_mae_lin_mol": round(float(np.mean(np.abs(errs["mon"]))), 3)})
    sens = pd.DataFrame(sens)

    # Table 5: FACE + ethanol, linear volume between the measured FACE+0 row and pure ethanol
    eth = pp["ethanol"]
    f_rows = []
    for base, sub in t5.groupby("base_fuel"):
        z = sub[sub["ethanol_vol_pct"] == 0].iloc[0]
        for r in sub[sub["ethanol_vol_pct"] > 0].itertuples():
            fe = r.ethanol_vol_pct / 100.0
            f_rows.append({"base_fuel": base, "ethanol_vol_pct": r.ethanol_vol_pct,
                           "ron_meas": r.ron, "mon_meas": r.mon,
                           "ron_lin_vol": (1 - fe) * z.ron + fe * eth["ron"],
                           "mon_lin_vol": (1 - fe) * z.mon + fe * eth["mon"]})
    face = pd.DataFrame(f_rows)
    face["ron_err"] = face["ron_lin_vol"] - face["ron_meas"]
    face["mon_err"] = face["mon_lin_vol"] - face["mon_meas"]
    face_m = {t: metrics(face[f"{t}_err"]) for t in ("ron", "mon")}

    df.to_csv(OUT_DIR / "octane_blending_predictions.csv", index=False)
    summary.to_csv(OUT_DIR / "octane_blending_summary.csv", index=False)
    fam.to_csv(OUT_DIR / "octane_blending_by_family.csv", index=False)
    face.to_csv(OUT_DIR / "octane_face_ethanol_linear_volume.csv", index=False)

    # ---------------- report: computed values only ----------------
    def S(t, m):
        return summary[(summary.target == t) & (summary.model == m)].iloc[0]

    lines = ["# KAUST octane blending benchmark (computed)", "",
             "All numbers below are produced by `scripts/external_ingest/eval_octane_blending_benchmarks.py`.",
             f"Data: KAUST Table 4, {len(full)} rows, of which {len(anchors)} single-component anchor rows "
             f"({', '.join(anchors['blend_name'])}) define the pure values and are excluded from the metrics; "
             f"{len(df)} blends scored. Validation: GroupKFold, {n_splits} splits, groups = component set with "
             f"ethanol removed ({df['group'].nunique()} groups, so ethanol blends share a fold with their base fuel).",
             "", "## Pure values used (from measured rows)", "", "| Component | RON | MON | Source | Density kg/m3 | Density source |",
             "|---|---|---|---|---|---|"]
    for k in NAMES:
        lines.append(f"| {k} | {pp[k]['ron']} | {pp[k]['mon']} | {pp[k]['source']} | {dens[k][0]:.1f} | {dens[k][1]} |")
    lines += ["", "Diisobutylene uses the Table 3 row 72 isomer (2,4,4-trimethyl-1-pentene) as a proxy; no pure "
              "diisobutylene row exists. n-Heptane uses the PRF scale definition (not tabulated in Table 3).",
              "", "## Overall metrics (ON)", "", "| Target | Model | n | MAE | RMSE | Bias | Max |", "|---|---|---|---|---|---|---|"]
    for r in summary.itertuples():
        lines.append(f"| {r.target} | {r.model} | {r.n} | {r.mae} | {r.rmse} | {r.bias:+.3f} | {r.max_err} |")
    lines += ["", "## Skill checks (computed)", ""]
    for t in ("RON", "MON"):
        lm, c, g, lv = (S(t, "Linear molar"), S(t, "Constant (train mean)"),
                        S(t, "Linear molar + GBM residual (GroupKFold)"), S(t, "Linear volume"))
        lines.append(f"- {t}: linear molar MAE {lm.mae} vs constant {c.mae} -> "
                     f"{'beats' if lm.mae < c.mae else 'does not beat'} the constant; "
                     f"GBM residual MAE {g.mae} vs linear molar {lm.mae} -> "
                     f"{'improves' if g.mae < lm.mae else 'does not improve'} on it; "
                     f"linear volume MAE {lv.mae} vs linear molar {lm.mae}.")
    lines += ["", "## Per base-fuel family (MAE, ON)", "",
              "| Group | n | RON lin. mol | RON lin. vol | RON GBM | MON lin. mol | MON lin. vol | MON GBM |",
              "|---|---|---|---|---|---|---|---|"]
    for r in fam.itertuples():
        lines.append(f"| `{r.group}` | {r.n} | {r.ron_mae_lin_mol} | {r.ron_mae_lin_vol} | {r.ron_mae_gbm} | "
                     f"{r.mon_mae_lin_mol} | {r.mon_mae_lin_vol} | {r.mon_mae_gbm} |")
    lines += ["", "## Density sensitivity (linear molar MAE, diisobutylene density scaled)", "",
              "| DIB density scale | RON MAE | MON MAE |", "|---|---|---|"]
    for r in sens.itertuples():
        lines.append(f"| {r.dib_density_scale} | {r.ron_mae_lin_mol} | {r.mon_mae_lin_mol} |")
    lines += ["", "## FACE gasoline + ethanol (Table 5), linear volume between measured FACE+0 and ethanol", "",
              f"n = {len(face)} blends; RON MAE {face_m['ron']['mae']:.3f} (bias {face_m['ron']['bias']:+.3f}, "
              f"max {face_m['ron']['max_err']:.3f}); MON MAE {face_m['mon']['mae']:.3f} "
              f"(bias {face_m['mon']['bias']:+.3f}, max {face_m['mon']['max_err']:.3f}). No model is fitted here.",
              "", "## Scope", "",
              "Model-level benchmark on literature-compiled data. Oxygenate (ethanol) rows are in scope for octane only. "
              "No comparison to literature model errors is made because none is computed in this repository."]
    (OUT_DIR / "octane_blending_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))
    print(fam.to_string(index=False))
    print(sens.to_string(index=False))
    print(face_m)


if __name__ == "__main__":
    main()
