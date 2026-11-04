#!/usr/bin/env python
"""Physics baseline for flash point (methods document, section 4) evaluated on the unified data.

    python baselines_flash.py              # Le Chatelier + Raoult (gamma = 1) + Antoine
    python baselines_flash.py --lfl-temp   # also apply Zabetakis' LFL(T) = LFL25 * (1 - 0.000721 (T-25))
    python baselines_flash.py --window 30  # tolerance (K) outside the Antoine validity range before a record is flagged

Needs: data/unified/measurements_dedup.csv and compounds.csv (run build_unified.py, refit_antoine.py, fill_lfl.py first).

Model: flash point T solves   H(T) = sum_i x_i * Psat_i(T) / (P_total * LFL_i) = 1
       Psat_i = 10**(A - B/(T/K + C)) kPa  (compounds.csv, kPa convention),  P_total = 101.325 kPa,  LFL_i = vol%/100.
Outputs: data/unified/baseline_flash.csv (one row per record) and data/unified/BASELINE_FLASH_REPORT.md.

How to read the report: `err = predicted - measured` in C; **positive err = the model says it is SAFER than it is** (the dangerous
direction).  Errors are split by pure/mixture, by Antoine validity (is the prediction inside the range the constants were fitted on?),
by LFL quality (measured vs 0.55*Cst estimate), by chain length and by whether an aromatic is present.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from build_unified import md

UNI = Path("data/unified")
P_TOTAL = 101.325


def psat(T_K, A, B, C):
    d = T_K + C
    return np.inf if d <= 1.0 else 10.0 ** (A - B / d)


def hazard(T_K, x, A, B, C, lfl_frac, lfl_temp):
    """H(T) - 1 for one blend."""
    tc = T_K - 273.15
    fac = max(1.0 - 0.000721 * (tc - 25.0), 0.05) if lfl_temp else 1.0
    return sum(xi * psat(T_K, a, b, c) / (P_TOTAL * l * fac) for xi, a, b, c, l in zip(x, A, B, C, lfl_frac)) - 1.0


def predict_flash_C(x, A, B, C, lfl_vol_pct, lfl_temp=False, lo=170.0, hi=700.0):
    """Flash point in C, or (nan, reason)."""
    lfl = np.asarray(lfl_vol_pct, float) / 100.0
    f = lambda T: hazard(T, x, A, B, C, lfl, lfl_temp)
    try:
        flo, fhi = f(lo), f(hi)
        if not np.isfinite(flo) or flo > 0:
            return np.nan, "flash point below 170 K bracket"
        if not np.isfinite(fhi) or fhi < 0:
            return np.nan, "no root below 700 K"
        return brentq(f, lo, hi, xtol=1e-4) - 273.15, ""
    except Exception as e:  # noqa: BLE001
        return np.nan, f"solver: {type(e).__name__}"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--lfl-temp", action="store_true")
    ap.add_argument("--window", type=float, default=20.0)
    a = ap.parse_args(argv)
    comp = pd.read_csv(UNI / "compounds.csv").set_index("inchikey")
    dd = pd.read_csv(UNI / "measurements_dedup.csv")
    fl = dd[(dd.property == "flash_point_C") & dd.is_hc.astype(bool) & (dd.kind != "real_fuel")].copy()
    need = ["antoine_A", "antoine_B", "antoine_C", "lfl_vol_pct"]

    rows = []
    for r in fl.itertuples():
        keys, x = json.loads(r.comp_inchikeys), json.loads(r.x_mole)
        cc = comp.loc[keys]
        rec = dict(system_id=r.system_id, kind=r.kind, n_comp=r.n_comp, measured_C=r.value, n_obs=r.n_obs, conflict=r.conflict,
                   max_n_carbon=int(cc.n_carbon.max()), min_n_carbon=int(cc.n_carbon.min()), has_aromatic=bool((cc.cls == "aromatic").any()),
                   classes="+".join(sorted(set(cc.cls))), comp_inchikeys=r.comp_inchikeys, x_mole=r.x_mole)
        if cc[need].isna().any().any():
            miss = [c for c in need if cc[c].isna().any()]
            rows.append({**rec, "predicted_C": np.nan, "status": "missing " + ",".join(miss)}); continue
        pred, why = predict_flash_C(x, cc.antoine_A.values, cc.antoine_B.values, cc.antoine_C.values, cc.lfl_vol_pct.values, a.lfl_temp)
        if np.isnan(pred):
            rows.append({**rec, "predicted_C": np.nan, "status": why}); continue
        Tp = pred + 273.15
        in_rng = bool(((Tp >= cc.antoine_Tlo - a.window) & (Tp <= cc.antoine_Thi + a.window)).all())
        est = cc.lfl_source.fillna("").str.startswith("estimated") if "lfl_source" in cc else pd.Series(True, index=cc.index)
        rows.append({**rec, "predicted_C": pred, "err_C": pred - r.value, "antoine_in_range": in_rng,
                     "lfl_all_measured": bool(not est.any()), "status": "ok"})
    out = pd.DataFrame(rows)
    out.to_csv(UNI / "baseline_flash.csv", index=False)

    ok = out[out.status == "ok"].copy()
    L = ["# Flash point — physics baseline (Le Chatelier + Raoult + Antoine)\n",
         f"Records: {len(out)}   predicted: {len(ok)}   not predictable: {len(out) - len(ok)}   LFL(T) correction: {a.lfl_temp}\n"]
    if len(out) - len(ok):
        L.append("Why records could not be predicted:\n" + md(out[out.status != "ok"].groupby("status").size()))
    if ok.empty:
        L.append("\nNo predictable records - fill LFL (fill_lfl.py) and Antoine (refit_antoine.py) first.")
    else:
        def stats(g):
            e = g.err_C
            return pd.Series(dict(n=len(g), MAE=e.abs().mean(), RMSE=np.sqrt((e ** 2).mean()), bias=e.mean(), median_abs=e.abs().median(),
                                  over_pred_pct=100 * (e > 0).mean(), over_by_5C_pct=100 * (e > 5).mean(), within_5C_pct=100 * (e.abs() <= 5).mean()))
        L.append("\n## Overall\n" + md(stats(ok).round(2).to_frame("value")))
        for title, col in (("by kind (pure / binary / multi)", "kind"), ("Antoine valid at the predicted T?", "antoine_in_range"),
                           ("LFL all measured (True) or any estimated (False)", "lfl_all_measured"), ("aromatic present", "has_aromatic")):
            L.append(f"\n## {title}\n" + md(ok.groupby(col).apply(stats, include_groups=False).round(2)))
        ok["max_C_bin"] = pd.cut(ok.max_n_carbon, [0, 5, 8, 12, 16, 100], labels=["<=C5", "C6-8", "C9-12", "C13-16", ">C16"])
        L.append("\n## by heaviest component (carbon number)\n" + md(ok.groupby("max_C_bin", observed=True).apply(stats, include_groups=False).round(2)))
        L.append("\n## by class composition\n" + md(ok.groupby("classes").apply(stats, include_groups=False).round(2).sort_values("n", ascending=False).head(12)))
        L.append("\n## worst 15 records (largest |error|)\n" + md(ok.reindex(ok.err_C.abs().sort_values(ascending=False).index).head(15)[
            ["kind", "n_comp", "classes", "max_n_carbon", "measured_C", "predicted_C", "err_C", "antoine_in_range", "lfl_all_measured", "conflict"]].round(1), index=False))
        mix = ok[ok.kind != "pure"]
        if len(mix):
            L.append("\n## Does the ideal-mixing assumption hold? (mixtures only)\n"
                     f"mixtures predicted: {len(mix)};  MAE {mix.err_C.abs().mean():.1f} C;  mean err {mix.err_C.mean():+.1f} C.  "
                     "Compare with the pure-compound MAE above: a clearly larger mixture error means non-ideality (gamma != 1) is the thing to learn.")
    (UNI / "BASELINE_FLASH_REPORT.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwrote {UNI/'baseline_flash.csv'} and {UNI/'BASELINE_FLASH_REPORT.md'}")


if __name__ == "__main__":
    main()
