#!/usr/bin/env python
"""Physics baselines + learned residuals for mixture properties.

    python baselines.py            # viscosity baseline + system-grouped residual model

Physics (no fitting):
  viscosity log-linear: ln eta_mix = sum w_i ln eta_i  (Arrhenius-ish)
  viscosity Refutas   : VBN_i = 14.534 ln ln(nu_i+0.8)+10.975, nu_i = eta_i/rho_i
                        eta_pred = nu_mix * rho_mix (Amagat density)
Pure-component eta/rho taken from the unified pure rows, nearest in T (<= 7.5 K).

Residual: HistGradientBoosting on composition aggregates + the physics value,
GroupKFold 5-fold on system_id (no system in both train and test).

Outputs: data/unified/baseline_viscosity.csv, BASELINE_VISCOSITY_REPORT.md
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

UNI = Path("data/unified")
HC_CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]
FEATURES = HC_CLASSES + ["other", "mean_mw", "mean_nc", "max_nc", "n_comp", "T_K", "log_eta_base"]


def parse_lists(s: str):
    return json.loads(s) if isinstance(s, str) and s else []


def nearest(tab: dict, T: float, tol: float = 7.5):
    best = min(tab.items(), key=lambda kv: abs(kv[0] - T)) if tab else None
    return best[1] if best and abs(best[0] - T) <= tol else None


def main():
    d = pd.read_csv(UNI / "measurements_dedup.csv")
    comp = pd.read_csv(UNI / "compounds.csv")
    comp = comp[comp.is_hc.astype(bool)]
    comp_lookup = {r.inchikey: dict(cls=r.cls, mw=r.mw, n_carbon=r.n_carbon)
                   for r in comp.itertuples()}

    pure = d[(d.is_hc) & (d.kind == "pure") & (d.property == "viscosity_dyn_mPas")]
    dens = d[(d.is_hc) & (d.kind == "pure") & (d.property == "density_kg_m3")]
    mix = d[(d.is_hc) & (d.n_comp >= 2) & (d.property == "viscosity_dyn_mPas")]

    eta_tab: dict[str, dict[float, float]] = {}
    for r in pure.itertuples():
        for ik in parse_lists(r.comp_inchikeys)[:1]:
            eta_tab.setdefault(ik, {})[round(float(r.T_K), 1)] = float(r.value)
    rho_tab: dict[str, dict[float, float]] = {}
    for r in dens.itertuples():
        for ik in parse_lists(r.comp_inchikeys)[:1]:
            rho_tab.setdefault(ik, {})[round(float(r.T_K), 1)] = float(r.value)

    # --- baseline pass over mixtures ---------------------------------------
    fail = defaultdict(int)
    recs = []
    for r in mix.itertuples():
        ik = parse_lists(r.comp_inchikeys)
        w = parse_lists(r.w_mass)
        if len(ik) < 2 or len(ik) != len(w) or abs(sum(w) - 1) > 0.2 or np.isnan(r.T_K):
            fail["shape"] += 1
            continue
        etas, rhos, ws, vbns = [], [], [], []
        ok = True
        for i, w_i in zip(ik, w):
            e = nearest(eta_tab.get(i, {}), r.T_K)
            rho = nearest(rho_tab.get(i, {}), r.T_K)
            if not e or not rho or e <= 0:
                fail["no_pure"] += 1
                ok = False
                break
            nu = e / rho * 1000.0   # mPa.s / (kg/m3) *1000 -> cSt (mm2/s)
            if nu <= 0.2:
                fail["nu_low"] += 1
                ok = False
                break
            etas.append(e)
            rhos.append(rho)
            ws.append(w_i)
            vbns.append(w_i * (14.534 * np.log(np.log(nu + 0.8)) + 10.975))
        if not ok or not len(ws):
            continue
        W = sum(ws)
        eta_loglin = float(np.exp(np.dot(np.array(ws) / W, np.log(etas))))
        rho_mix = 1.0 / sum(w / r2 for w, r2 in zip(ws, rhos))
        nu_mix = np.exp(np.exp((sum(vbns) / W - 10.975) / 14.534)) - 0.8
        eta_refutas = float(nu_mix * rho_mix / 1000.0)   # cSt * kg/m3 /1000 -> mPa.s
        recs.append(dict(system_id=r.system_id, kind=r.kind, n_comp=int(r.n_comp),
                         T_K=float(r.T_K), observed=float(r.value),
                         eta_logline=eta_loglin, eta_refutas=eta_refutas,
                         comp_inchikeys=r.comp_inchikeys, w_mass=r.w_mass, sources=r.sources))
    df = pd.DataFrame(recs)
    print("debug fail counts:", dict(fail))
    print(f"viscosity baseline: {len(df)} mixture records with full pure-component "
          f"coverage (of {len(mix)} mixture records)")

    # --- composition features + residual model -----------------------------
    X, y, groups, base = [], [], [], []
    for r in df.itertuples():
        ik = parse_lists(r.comp_inchikeys)
        w = parse_lists(r.w_mass)
        f = {c: 0.0 for c in HC_CLASSES}
        f["other"] = 0.0
        mw_w, nc_w, wt = 0.0, 0.0, 0.0
        max_nc = 0
        for i, w_i in zip(ik, w):
            c = comp_lookup.get(i, {})
            cls = c.get("cls")
            if cls in f:
                f[cls] += w_i
            else:
                f["other"] += w_i
            mw, nc = c.get("mw", np.nan), c.get("n_carbon", np.nan) or 0
            if pd.notna(mw):
                mw_w += w_i * mw
                wt += w_i
            if nc:
                nc_w += w_i * nc
                max_nc = max(max_nc, nc)
        f["mean_mw"] = mw_w / wt if wt else np.nan
        f["mean_nc"] = nc_w / wt if wt else np.nan
        f["max_nc"] = max_nc or np.nan
        f["n_comp"] = len(ik)
        f["T_K"] = r.T_K
        f["log_eta_base"] = np.log(r.eta_logline)
        X.append([f[k] for k in FEATURES])
        y.append(np.log(r.observed) - np.log(r.eta_logline))
        groups.append(r.system_id)
        base.append(r.eta_logline)
    X, y = np.array(X, dtype=float), np.array(y, dtype=float)
    base = np.array(base)

    gk = GroupKFold(n_splits=5)
    pred_res = np.zeros_like(y)
    for tr, te in gk.split(X, y, groups=groups):
        # impute NaNs by column median of the train fold
        med = np.nanmedian(X[tr], axis=0)
        Xi_tr = np.where(np.isnan(X[tr]), med, X[tr])
        Xi_te = np.where(np.isnan(X[te]), med, X[te])
        m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06)
        m.fit(Xi_tr, y[tr])
        pred_res[te] = m.predict(Xi_te)
    obs, pl = df.observed.values, base
    pred = np.exp(np.log(pl) + pred_res)
    df["pred_gbm"] = pred
    df["err_logline_pct"] = 100 * (pl - obs) / obs
    df["err_refutas_pct"] = 100 * (df.eta_refutas - obs) / obs
    df["err_gbm_pct"] = 100 * (pred - obs) / obs

    def mape(col):
        a = df[col].abs()
        return 100 * float((df[col.replace("err", "abs")]).mean()) if False else float(a.mean())

    for split in ("all",):
        s = df
        print(f"[{split}] n={len(s)}  loglin MAPE={s.err_logline_pct.abs().mean():.2f}%  "
              f"refutas MAPE={s.err_refutas_pct.abs().mean():.2f}%  GBM residual MAPE={s.err_gbm_pct.abs().mean():.2f}%")
    by = df.groupby("kind").apply(lambda s: pd.Series({
        "n": len(s), "loglin": s.err_logline_pct.abs().mean(),
        "refutas": s.err_refutas_pct.abs().mean(), "gbm": s.err_gbm_pct.abs().mean()}))
    print(by.to_string())
    df.to_csv(UNI / "baseline_viscosity.csv", index=False)
    lines = [
        "# Viscosity baseline + residual model",
        "",
        f"records: {len(df)} (GroupKFold on system_id, 5 folds)",
        "",
        "| slice | n | loglin MAPE% | refutas MAPE% | GBM residual MAPE% |",
        "|---|---|---|---|---|",
        f"| all | {len(df)} | {df.err_logline_pct.abs().mean():.2f} | {df.err_refutas_pct.abs().mean():.2f} | {df.err_gbm_pct.abs().mean():.2f} |",
    ]
    for k, row in by.iterrows():
        lines.append(f"| {k} | {row['n']:.0f} | {row['loglin']:.2f} | {row['refutas']:.2f} | {row['gbm']:.2f} |")
    (UNI / "BASELINE_VISCOSITY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {UNI/'baseline_viscosity.csv'} and BASELINE_VISCOSITY_REPORT.md")


if __name__ == "__main__":
    main()
