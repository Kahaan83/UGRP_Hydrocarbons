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


# ---- flash residual diagnostics + README headline block (regenerated between the HEADLINE markers) -----------------
def _diag():
    pr = pd.read_csv(ROOT / "data" / "real_fuels" / "v4" / "jet_properties_v4.csv")
    ar = pr[(pr.property == "aromatics_vol_pct") & pr.fuel_id.str.startswith("CRC647")].groupby("fuel_id").value.mean()
    t10 = pr[(pr.property == "distillation_T10_C") & pr.method.str.startswith("ASTM D86")].groupby("fuel_id").value.mean()
    o = {}
    for k in ("v4", "r2"):
        f = fl[k][fl[k].fuel_id.str.startswith("CRC647") & fl[k].statistics_allowed].set_index("fuel_id")
        d = f[f.split_role != "locked"]
        o[k] = (d.ideal_err_C.corr(ar.reindex(d.index)), d.ideal_err_C.corr(t10.reindex(d.index)),
                d.ideal_err_C.corr(d.unexplained_fraction_dropped))
    return o


diag = _diag()
for k, names in (("v4", None), ("r2", None)):
    for nm, v in zip(("aromatics_vol_pct", "D86_T10", "dropped_lump_mass"), diag[k]):
        rows.append(dict(property="flash_C", slice="dev", quantity=f"corr(raw ideal err, {nm})", v4=diag["v4"][("aromatics_vol_pct", "D86_T10", "dropped_lump_mass").index(nm)],
                         r2=diag["r2"][("aromatics_vol_pct", "D86_T10", "dropped_lump_mass").index(nm)])) if k == "v4" else None
pd.DataFrame(rows).to_csv(R2 / "readings.csv", index=False)

sens = pd.read_csv(R2 / "sat_split_sensitivity.csv")
sens["err"] = np.where(sens.property.str.startswith("visc"), np.log(sens.predicted / sens.measured), sens.predicted - sens.measured)
within = {p: int((sens[(sens.draw >= 0) & (sens.property == p)].groupby("draw").err.mean().abs() <= 3).sum())
          for p in ("freezing", "flash", "density_15.56C")}


def g(prop, sl, q, col):
    return float(out[(out.property == prop) & (out.slice == sl) & (out.quantity == q)][col].iloc[0])


head = f"""**Headline findings (generated by `scripts/v4/r2_readings.py`)**

* **Flash (safety-relevant): the raw ideal Le Chatelier bias flips from {g('flash_C','dev','raw ideal bias','v4'):+.1f} C (v4, under-predicts) to {g('flash_C','dev','raw ideal bias','r2'):+.1f} C (r2, over-predicts, the unsafe side); {g('flash_C','dev','over-prediction rate >5 C (raw)','r2')*100:.0f}% of dev fuels are over-predicted by more than 5 C (locked {g('flash_C','locked','raw ideal bias','r2'):+.1f} C).** The saturates split does not drive this ({within['flash']} of 32 draws within +/-3 C of zero bias); representative choice and lump handling were not varied.
* Density 15.56 C: dev bias {g('density_15.56C_kg_m3','dev','bias','v4'):+.2f} -> {g('density_15.56C_kg_m3','dev','bias','r2'):+.2f} kg/m3, still outside +/-3; only {within['density_15.56C']} of 32 split draws within +/-3.
* Freezing (ideal SLE, mole): dev bias {g('freezing_C','dev','ideal_sle_mole bias','v4'):+.2f} -> {g('freezing_C','dev','ideal_sle_mole bias','r2'):+.2f} C; {within['freezing']} of 32 split draws within +/-3 C; beating the dev-mean constant not established (bootstrap CI includes 0).
* Viscosity (raw Refutas, dev, ln): bias {g('viscosity_ln','dev','raw Refutas bias_log [overall (all T)]','v4'):+.3f} -> {g('viscosity_ln','dev','raw Refutas bias_log [overall (all T)]','r2'):+.3f}, not below 0.10, so the recalibration was not just correcting the reconstruction.
"""
rd = R2 / "README.md"
txt = rd.read_text(encoding="utf-8")
a, b = txt.index("<!-- HEADLINE:START -->"), txt.index("<!-- HEADLINE:END -->")
rd.write_text(txt[:a] + "<!-- HEADLINE:START -->\n" + head + "\n" + txt[b:], encoding="utf-8")
