"""v5 delta_ws3_pure_gaps Part A: heavy n-alkane Antoine by carbon-number interpolation, checked against ThermoML vapour pressures.

Declared readings: results/v5/BASE/delta_ws3_pure_gaps/README.md (committed before this script existed).
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_ws3_antoine.py
Writes antoine_knots.csv, antoine_loo.csv, antoine_targets.csv, readings_antoine.json in the delta folder and
data/processed/v5/antoine_interpolated_nalkanes.csv.  Changes no registry.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.thermo.antoine_interp import SOURCE, T_GRID, grid_interp, log10p, param_interp, rms_log10  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_ws3_pure_gaps"
OUT_TABLE = ROOT / "data" / "processed" / "v5" / "antoine_interpolated_nalkanes.csv"
TARGETS = [13, 15, 17, 18, 19]
GAP4 = [12, 14]
THRESH = 0.03
TLO, THI = float(T_GRID.min()), float(T_GRID.max())


def nalkane_smiles(n):
    return "C" * n


def load():
    comp = pd.read_csv(ROOT / "data/unified/compounds.csv")
    ik2smi = dict(zip(comp.inchikey, comp.smiles))
    d = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=["property", "kind", "comp_inchikeys", "T_K", "value", "sources"])
    v = d[(d.property == "vapor_pressure_kPa") & (d.kind == "pure")].copy()
    v["smiles"] = [ik2smi.get(json.loads(k)[0]) for k in v.comp_inchikeys]
    v = v[v.smiles.astype(str).str.fullmatch(r"C+")]
    v["n"] = v.smiles.str.len()
    refit = pd.read_csv(ROOT / "data/processed/antoine_refit.csv")
    refit = refit[refit.smiles.astype(str).str.fullmatch(r"C+")].copy()
    refit["n"] = refit.smiles.str.len()
    reg = pd.read_csv(ROOT / "data/processed/v2/compounds_anchored_vp.csv")
    reg = reg[reg.smiles.astype(str).str.fullmatch(r"C+")].copy()
    reg["n"] = reg.smiles.str.len()
    return v, refit.set_index("n"), reg.set_index("n")


def exploratory_fresh_refit(meas):
    """POST-HOC, NOT A GATE (added after R-W3-2 failed): knots refit freshly on each compound's own measured ThermoML points (>= 6
    points, own-fit rms <= 0.03), then the same leave-one-knot-out.  Shows whether the failure is the stale antoine_refit.csv or the method."""
    from fuelblend.thermo.antoine_interp import refit_abc
    fk = {}
    for n, g in meas.items():
        g = g[g.value > 0]
        if len(g) < 6 or not (8 <= n <= 22):
            continue
        abc = refit_abc(g.T_K, np.log10(g.value))
        if rms_log10(abc, g.T_K, g.value) <= THRESH:
            fk[n] = abc
    out = []
    ks = sorted(fk)
    for n in ks[1:-1]:
        rest = {k: v_ for k, v_ in fk.items() if k != n}
        g = meas[n]
        g = g[(g.T_K >= TLO) & (g.T_K <= THI) & (g.value > 0)]
        abc_g, _ = grid_interp(rest, n)
        e = log10p(abc_g, g.T_K) - np.log10(g.value)
        lo_k, hi_k = max(k for k in rest if k < n), min(k for k in rest if k > n)
        out.append(dict(knots_fresh=json.dumps(ks), n_held_out=n, left_knot=lo_k, right_knot=hi_k, gap=hi_k - lo_k, n_points=len(g),
                        rms_grid=float(np.sqrt(np.mean(e ** 2))) if len(g) else np.nan, bias=float(e.mean()) if len(g) else np.nan))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v, refit, reg = load()
    meas = {n: g for n, g in v.groupby("n")}
    rows = []
    for n in range(8, 23):
        g = meas.get(n)
        r = refit.loc[n] if n in refit.index else None
        abc = None if r is None or not np.isfinite([r.A, r.B, r.C]).all() else (r.A, r.B, r.C)
        npts = 0 if g is None else len(g)
        rms = rms_log10(abc, g.T_K, g.value) if (abc is not None and npts) else np.nan
        knot = bool(abc is not None and npts >= 6 and np.isfinite(rms) and rms <= THRESH)
        rows.append(dict(n=n, n_measured=npts, T_meas_min=None if g is None else g.T_K.min(), T_meas_max=None if g is None else g.T_K.max(),
                         refit_status=None if r is None else r.status, refit_n_points=None if r is None else r.n_points,
                         refit_T_min=None if r is None else r.T_min, refit_T_max=None if r is None else r.T_max,
                         rms_refit_on_own_measured=rms, knot=knot))
    knots_df = pd.DataFrame(rows)
    knots_df.to_csv(OUT / "antoine_knots.csv", index=False)
    knots = {int(r.n): tuple(refit.loc[int(r.n), ["A", "B", "C"]]) for r in knots_df[knots_df.knot].itertuples()}
    kn = sorted(knots)

    # R-W3-2 leave-one-knot-out over interior knots
    loo, allres_grid, allres_par = [], [], []
    for n in kn[1:-1]:
        rest = {k: v_ for k, v_ in knots.items() if k != n}
        g = meas[n]
        g = g[(g.T_K >= TLO) & (g.T_K <= THI) & (g.value > 0)]
        abc_g, _ = grid_interp(rest, n)
        abc_p = param_interp(rest, n)
        lo_k, hi_k = max(k for k in rest if k < n), min(k for k in rest if k > n)
        eg = log10p(abc_g, g.T_K) - np.log10(g.value)
        ep = log10p(abc_p, g.T_K) - np.log10(g.value)
        allres_grid += list(eg)
        allres_par += list(ep)
        loo.append(dict(n_held_out=n, left_knot=lo_k, right_knot=hi_k, gap=hi_k - lo_k, n_points_in_grid=len(g),
                        rms_grid=float(np.sqrt(np.mean(eg ** 2))) if len(g) else np.nan, bias_grid=float(eg.mean()) if len(g) else np.nan,
                        rms_param=float(np.sqrt(np.mean(ep ** 2))) if len(g) else np.nan,
                        A=abc_g[0], B=abc_g[1], C=abc_g[2]))
    loo = pd.DataFrame(loo)
    loo.to_csv(OUT / "antoine_loo.csv", index=False)
    pooled = float(np.sqrt(np.mean(np.square(allres_grid)))) if allres_grid else np.nan
    pooled_p = float(np.sqrt(np.mean(np.square(allres_par)))) if allres_par else np.nan
    gap_ok = {}
    for n in GAP4:
        r = loo[loo.n_held_out == n]
        gap_ok[n] = dict(present=bool(len(r)), rms=float(r.rms_grid.iloc[0]) if len(r) else None,
                         gap=int(r.gap.iloc[0]) if len(r) else None,
                         passed=bool(len(r) and r.rms_grid.iloc[0] <= THRESH and r.gap.iloc[0] == 4))
    gate = bool(np.isfinite(pooled) and pooled <= THRESH and all(x["passed"] for x in gap_ok.values()))

    # targets
    trows = []
    for n in TARGETS + ([11] if 11 not in knots else []):
        abc_g, _ = grid_interp(knots, n)
        abc_p = param_interp(knots, n)
        rr = refit.loc[n] if n in refit.index else None
        rg = reg.loc[n] if n in reg.index else None
        rec = dict(n=n, smiles=nalkane_smiles(n), A=abc_g[0], B=abc_g[1], C=abc_g[2], T_lo=TLO, T_hi=THI)
        for T in (313.15, 373.15):
            rec[f"log10P_interp_{T:.0f}"] = float(log10p(abc_g, T))
            rec[f"log10P_param_{T:.0f}"] = float(log10p(abc_p, T))
            rec[f"log10P_refit_csv_{T:.0f}"] = float(log10p((rr.A, rr.B, rr.C), T)) if rr is not None and np.isfinite(rr.A) else np.nan
            rec[f"log10P_v2_registry_{T:.0f}"] = (float(log10p((rg.antoine_A, rg.antoine_B, rg.antoine_C), T))
                                                  if rg is not None and np.isfinite(rg.antoine_A) else np.nan)
        g = meas.get(n)
        if g is not None and len(g):
            rec["n_measured"] = len(g)
            rec["rms_vs_measured_any_T"] = rms_log10(abc_g, g.T_K, g.value)
            rec["measured_T_range"] = f"{g.T_K.min():.1f}-{g.T_K.max():.1f}"
        trows.append(rec)
    targ = pd.DataFrame(trows)
    targ.to_csv(OUT / "antoine_targets.csv", index=False)

    tab = targ[targ.n.isin(TARGETS)][["n", "smiles", "A", "B", "C", "T_lo", "T_hi"]].copy()
    tab["antoine_source"] = SOURCE
    tab["antoine_unit"] = "log10(P/kPa) = A - B/(T/K + C)"
    tab["knots"] = json.dumps(kn)
    tab["gate_R_W3_2"] = "passed" if gate else "failed"
    tab["note"] = "carbon-number PCHIP of log10 P on the T grid, refit; not wired into any registry" + ("" if gate else "; gate failed: not for use")
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUT_TABLE, index=False)

    explo = exploratory_fresh_refit(meas)
    explo.to_csv(OUT / "antoine_exploratory_fresh_refit.csv", index=False)

    readings = {"R-W3-1": dict(knots=kn, non_knots=knots_df[~knots_df.knot].n.tolist()),
                "exploratory_fresh_refit_NOT_A_GATE": explo.to_dict("records"),
                "R-W3-2": dict(pooled_rms_grid=pooled, pooled_rms_param_comparator=pooled_p, gap4_cases=gap_ok, threshold=THRESH, passed=gate),
                "C11_check": targ[targ.n == 11].to_dict("records")}
    (OUT / "readings_antoine.json").write_text(json.dumps(readings, indent=2, default=float), encoding="utf-8")
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(knots_df.to_string())
    print(loo.to_string())
    print(targ.to_string())
    print(json.dumps(readings, indent=2, default=float))


if __name__ == "__main__":
    main()
