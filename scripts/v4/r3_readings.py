"""v4 vs r2 vs r3 comparison table (reads existing evaluator outputs only; same formulas/seed as r2_readings.py).

Writes results/v4/R3/readings.csv with columns property, slice, quantity, v4, r2, r3.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DIRS = {"v4": ROOT / "results" / "v4", "r2": ROOT / "results" / "v4" / "R2", "r3": ROOT / "results" / "v4" / "R3"}
KS = ("v4", "r2", "r3")
rows = []


def add(prop, sl, q, a, b, c):
    rows.append(dict(property=prop, slice=sl, quantity=q, v4=a, r2=b, r3=c))


def lofo_const_mae(m):
    m = np.asarray(m, float)
    return float(np.mean([abs(np.delete(m, i).mean() - m[i]) for i in range(len(m))]))


# density at 15.56 C
d = {k: pd.read_csv(p / "B" / "density_summary_v4.csv") for k, p in DIRS.items()}
for sl in ("dev", "locked"):
    g = {k: v[v.scope.str.startswith("15.56C D4052") & (v.slice == sl)].iloc[0] for k, v in d.items()}
    for q in ("N", "bias", "bias_lo95", "bias_hi95", "MAE"):
        add("density_15.56C_kg_m3", sl, q, *[g[k][q] for k in KS])

# freezing
fz = {k: pd.read_csv(p / "B" / "freezing_v4.csv") for k, p in DIRS.items()}
cols = {"ideal_sle_mole": "err_ideal_sle_mole", "wat_ideal_mass_ln": "err_wat_ideal",
        "wat_minus_devfit_offset": "err_wat_minus_devfit_offset", "coutinho_library": "err_coutinho"}
for sl in ("dev", "locked"):
    for nm, c in cols.items():
        add("freezing_C", sl, f"{nm} bias", *[fz[k][fz[k].slice == sl][c].mean() for k in KS])
        add("freezing_C", sl, f"{nm} MAE", *[fz[k][fz[k].slice == sl][c].abs().mean() for k in KS])
    dm = fz["v4"][fz["v4"].slice == "dev"].measured_C.to_numpy()
    lk = fz["v4"][fz["v4"].slice == "locked"].measured_C.to_numpy()
    v = lofo_const_mae(dm) if sl == "dev" else float(np.abs(dm.mean() - lk).mean())
    add("freezing_C", sl, "dev-mean constant MAE (dev: leave-one-fuel-out; locked: dev mean)", v, v, v)


def boot(k):
    rng = np.random.default_rng(20261111)
    r = fz[k][fz[k].slice == "dev"]
    m = r.measured_C.to_numpy(); e = r.err_ideal_sle_mole.abs().to_numpy()
    c = np.array([abs(np.delete(m, i).mean() - m[i]) for i in range(len(m))])
    dd = e - c
    bs = [dd[rng.integers(0, len(dd), len(dd))].mean() for _ in range(5000)]
    return f"{np.percentile(bs, 2.5):.3f}/{dd.mean():.3f}/{np.percentile(bs, 97.5):.3f}"


add("freezing_C", "dev", "MAE(ideal SLE) - MAE(dev-mean const), fuel bootstrap lo95/mean/hi95", np.nan, boot("r2"), boot("r3"))

# viscosity
lb = {k: pd.read_csv(p / "V" / "viscosity_leaderboard_v4.csv") for k, p in DIRS.items()}
for sl, nm in (("dev", "dev (LOO by fuel)"), ("locked", "locked (untouched)")):
    for cond in ("overall (all T)", "T = 20 C", "T = -20 C", "T = -40 C"):
        for q in ("bias_log", "log_rmse"):
            vals = [lb[k][(lb[k].model == "Raw Refutas") & (lb[k].split == nm) & (lb[k].condition == cond)][q].iloc[0] for k in KS]
            add("viscosity_ln", sl, f"raw Refutas {q} [{cond}]", *vals)
    for mname in ("Recalibrated Refutas", "Structured ridge"):
        vals = [lb[k][(lb[k].model == mname) & (lb[k].split == nm) & (lb[k].condition == "overall (all T)")].log_rmse.iloc[0] for k in KS]
        add("viscosity_ln", sl, f"{mname} log_rmse [overall]", *vals)

# flash
fl = {k: pd.read_csv(p / "F" / "flash_real_fuels.csv") for k, p in DIRS.items()}
res = {}
for k, f in fl.items():
    f = f[f.fuel_id.str.startswith("CRC647") & f.statistics_allowed & (f.ideal_status == "ok")]
    dv, lk = f[f.split_role != "locked"], f[f.split_role == "locked"]
    e = dv.ideal_err_C.to_numpy(); mm = dv.measured_C.to_numpy()
    off = e.mean()
    res[k] = {("dev", "N"): len(dv), ("locked", "N"): len(lk),
              ("dev", "raw ideal bias"): e.mean(), ("locked", "raw ideal bias"): lk.ideal_err_C.mean(),
              ("dev", "raw ideal MAE"): np.abs(e).mean(), ("locked", "raw ideal MAE"): lk.ideal_err_C.abs().mean(),
              ("dev", "over-prediction rate >5 C (raw)"): float(np.mean(e > 5)), ("locked", "over-prediction rate >5 C (raw)"): float(np.mean(lk.ideal_err_C > 5)),
              ("dev", "dev-fit constant offset MAE (leave-one-fuel-out)"): float(np.mean([abs(e[i] - np.delete(e, i).mean()) for i in range(len(e))])),
              ("locked", "dev-fit constant offset MAE (offset from all dev)"): float((lk.ideal_err_C - off).abs().mean()),
              ("dev", "dev-mean constant MAE (leave-one-fuel-out)"): lofo_const_mae(mm),
              ("locked", "dev-mean constant MAE (dev mean)"): float(np.abs(mm.mean() - lk.measured_C).mean()),
              ("dev", "dev-fit offset value (C)"): off}
for key in res["v4"]:
    add("flash_C", key[0], key[1], *[res[k].get(key) for k in KS])

# flash residual correlations (dev)
pr = pd.read_csv(ROOT / "data" / "real_fuels" / "v4" / "jet_properties_v4.csv")
ar = pr[(pr.property == "aromatics_vol_pct") & pr.fuel_id.str.startswith("CRC647")].groupby("fuel_id").value.mean()
t10 = pr[(pr.property == "distillation_T10_C") & pr.method.str.startswith("ASTM D86")].groupby("fuel_id").value.mean()
diag = {}
for k in KS:
    f = fl[k][fl[k].fuel_id.str.startswith("CRC647") & fl[k].statistics_allowed].set_index("fuel_id")
    dv = f[f.split_role != "locked"]
    diag[k] = (dv.ideal_err_C.corr(ar.reindex(dv.index)), dv.ideal_err_C.corr(t10.reindex(dv.index)),
               dv.ideal_err_C.corr(dv.unexplained_fraction_dropped))
for i, nm in enumerate(("aromatics_vol_pct", "D86_T10", "dropped_lump_mass")):
    add("flash_C", "dev", f"corr(raw ideal err, {nm})", *[diag[k][i] for k in KS])

out = pd.DataFrame(rows)
out.to_csv(DIRS["r3"] / "readings.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 90, "display.max_rows", 300)
print(out.to_string(index=False))

# reproduction check against results/v4/R2/readings.csv
old = pd.read_csv(DIRS["r2"] / "readings.csv")
new = out[["property", "slice", "quantity", "v4", "r2"]]
mg = old.merge(new, on=["property", "slice", "quantity"], suffixes=("_old", "_new"), how="outer", indicator=True)
bad = 0
for _, r in mg.iterrows():
    if r["_merge"] != "both":
        print("ROW MISMATCH", r["property"], r["slice"], r["quantity"], r["_merge"]); bad += 1; continue
    for c in ("v4", "r2"):
        a, b = r[f"{c}_old"], r[f"{c}_new"]
        if isinstance(a, str) or isinstance(b, str):
            same = str(a) == str(b)
        else:
            same = (pd.isna(a) and pd.isna(b)) or float(a) == float(b)
        if not same:
            print("CELL DIFF", r["property"], r["slice"], r["quantity"], c, a, b); bad += 1
print(f"reproduction check vs results/v4/R2/readings.csv: {len(old)} old rows, {bad} differing cells/rows")
