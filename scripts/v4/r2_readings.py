"""v4 review fix 2: v4-vs-r2 comparison table behind results/v4/R2/README.md (reads existing outputs only).

Writes results/v4/R2/readings.csv: (property, slice, quantity, v4, r2).  The pre-declared readings (README) use the
dev rows; locked rows are listed beside them for information only and no rule was chosen from them.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
V4, R2 = ROOT / "results" / "v4", ROOT / "results" / "v4" / "R2"
rows = []


def add(prop, sl, q, a, b):
    rows.append(dict(property=prop, slice=sl, quantity=q, v4=a, r2=b))


def lofo_const_mae(m):
    m = np.asarray(m, float)
    return float(np.mean([abs(np.delete(m, i).mean() - m[i]) for i in range(len(m))]))


# density at 15.56 C
d = {k: pd.read_csv(p / "B" / "density_summary_v4.csv") for k, p in (("v4", V4), ("r2", R2))}
for sl in ("dev", "locked"):
    g = {k: v[v.scope.str.startswith("15.56C D4052") & (v.slice == sl)].iloc[0] for k, v in d.items()}
    for q in ("N", "bias", "bias_lo95", "bias_hi95", "MAE"):
        add("density_15.56C_kg_m3", sl, q, g["v4"][q], g["r2"][q])

# freezing (ideal SLE mole fraction, WAT variants, constant)
fz = {k: pd.read_csv(p / "B" / "freezing_v4.csv") for k, p in (("v4", V4), ("r2", R2))}
cols = {"ideal_sle_mole": "err_ideal_sle_mole", "wat_ideal_mass_ln": "err_wat_ideal",
        "wat_minus_devfit_offset": "err_wat_minus_devfit_offset", "coutinho_library": "err_coutinho"}
for sl in ("dev", "locked"):
    for nm, c in cols.items():
        add("freezing_C", sl, f"{nm} bias", *[fz[k][fz[k].slice == sl][c].mean() for k in ("v4", "r2")])
        add("freezing_C", sl, f"{nm} MAE", *[fz[k][fz[k].slice == sl][c].abs().mean() for k in ("v4", "r2")])
    dm = fz["v4"][fz["v4"].slice == "dev"].measured_C.to_numpy()
    lk = fz["v4"][fz["v4"].slice == "locked"].measured_C.to_numpy()
    add("freezing_C", sl, "dev-mean constant MAE (dev: leave-one-fuel-out; locked: dev mean)",
        *([lofo_const_mae(dm)] * 2 if sl == "dev" else [float(np.abs(dm.mean() - lk).mean())] * 2))
rng = np.random.default_rng(20261111)
r = fz["r2"][fz["r2"].slice == "dev"]
m = r.measured_C.to_numpy(); e = r.err_ideal_sle_mole.abs().to_numpy()
c = np.array([abs(np.delete(m, i).mean() - m[i]) for i in range(len(m))])
dd = e - c
bs = [dd[rng.integers(0, len(dd), len(dd))].mean() for _ in range(5000)]
add("freezing_C", "dev", "MAE(ideal SLE) - MAE(dev-mean const), fuel bootstrap lo95/mean/hi95",
    np.nan, f"{np.percentile(bs, 2.5):.3f}/{dd.mean():.3f}/{np.percentile(bs, 97.5):.3f}")

# viscosity: raw Refutas
lb = {k: pd.read_csv(p / "V" / "viscosity_leaderboard_v4.csv") for k, p in (("v4", V4), ("r2", R2))}
for sl, nm in (("dev", "dev (LOO by fuel)"), ("locked", "locked (untouched)")):
    for cond in ("overall (all T)", "T = 20 C", "T = -20 C", "T = -40 C"):
        for q in ("bias_log", "log_rmse"):
            vals = [lb[k][(lb[k].model == "Raw Refutas") & (lb[k].split == nm) & (lb[k].condition == cond)][q].iloc[0]
                    for k in ("v4", "r2")]
            add("viscosity_ln", sl, f"raw Refutas {q} [{cond}]", *vals)
    for mname in ("Recalibrated Refutas", "Structured ridge"):
        vals = [lb[k][(lb[k].model == mname) & (lb[k].split == nm) & (lb[k].condition == "overall (all T)")].log_rmse.iloc[0]
                for k in ("v4", "r2")]
        add("viscosity_ln", sl, f"{mname} log_rmse [overall]", *vals)

# flash (ideal Le Chatelier)
fl = {k: pd.read_csv(p / "F" / "flash_real_fuels.csv") for k, p in (("v4", V4), ("r2", R2))}
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
    add("flash_C", key[0], key[1], res["v4"][key], res["r2"].get(key))
out = pd.DataFrame(rows)
out.to_csv(R2 / "readings.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 90, "display.max_rows", 200)
print(out.to_string(index=False))
