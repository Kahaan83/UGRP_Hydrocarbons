"""Sensitivity of the r3 scores to the n-paraffin share prior (the only prior r3 still assumes).

For each n_share: rebuild the r3 composition (rebuild_crc_reconstruction_r3.py --n-share ...) under results/v4/R3/sens/,
run the bulk and flash evaluators (and the viscosity leaderboard with --viscosity), then write
results/v4/R3/n_share_sensitivity.csv.  Nothing is tuned; the evaluators are run unchanged.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PY = sys.executable
SENS = ROOT / "results" / "v4" / "R3" / "sens"
GRID = [("0.10", 0.10), ("0.20", 0.20), ("0.288", None), ("0.40", 0.40), ("0.50", 0.50)]   # None = nominal GC prior (0.2881)


def run(cmd):
    print(">>", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], cwd=ROOT, check=True)


def build_and_eval(tag, val, viscosity):
    d = SENS / f"ns_{tag}"
    d.mkdir(parents=True, exist_ok=True)
    comp = d / "compositions.csv"
    cmd = [PY, "scripts/v4/rebuild_crc_reconstruction_r3.py", "--out-csv", comp, "--out-dir", d]
    if val is not None:
        cmd += ["--n-share", val]
    run(cmd)
    run([PY, "scripts/v4/eval_bulk_props_v4.py", "--compositions", comp, "--out", d / "B"])
    run([PY, "scripts/v4/eval_flash_v4.py", "--compositions", comp, "--out", d / "F"])
    if viscosity:
        run([PY, "scripts/v4/build_viscosity_leaderboard_v4.py", "--compositions", comp, "--out", d / "V"])


def summarise(viscosity):
    rows = []
    for tag, val in GRID:
        d = SENS / f"ns_{tag}"
        ns = val if val is not None else 0.2881   # nominal GC prior as printed by the rebuild script
        ds = pd.read_csv(d / "B" / "density_summary_v4.csv")
        fz = pd.read_csv(d / "B" / "freezing_v4.csv")
        fl = pd.read_csv(d / "F" / "flash_real_fuels.csv")
        fl = fl[fl.fuel_id.str.startswith("CRC647") & fl.statistics_allowed & (fl.ideal_status == "ok")]
        lb = pd.read_csv(d / "V" / "viscosity_leaderboard_v4.csv") if viscosity else None
        for sl in ("dev", "locked"):
            g = ds[ds.scope.str.startswith("15.56C D4052") & (ds.slice == sl)].iloc[0]
            z = fz[fz.slice == sl].err_ideal_sle_mole
            f = fl[(fl.split_role != "locked") if sl == "dev" else (fl.split_role == "locked")].ideal_err_C
            r = dict(n_share=tag, n_share_value=ns, slice=sl,
                     density_15p56C_bias_kg_m3=g["bias"], density_15p56C_MAE_kg_m3=g["MAE"], density_N=g["N"],
                     freezing_ideal_sle_mole_bias_C=z.mean(), freezing_ideal_sle_mole_MAE_C=z.abs().mean(),
                     flash_raw_ideal_bias_C=f.mean(), flash_raw_ideal_MAE_C=f.abs().mean(),
                     flash_frac_overpredicted_gt5C=float((f > 5).mean()), flash_N=len(f))
            if viscosity:
                nm = "dev (LOO by fuel)" if sl == "dev" else "locked (untouched)"
                q = lb[(lb.model == "Raw Refutas") & (lb.split == nm) & (lb.condition == "overall (all T)")].iloc[0]
                r.update(raw_refutas_overall_bias_ln=q["bias_log"], raw_refutas_overall_log_rmse=q["log_rmse"])
            rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results" / "v4" / "R3" / "n_share_sensitivity.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--viscosity", action="store_true", help="also run the viscosity leaderboard per n_share (slow)")
    ap.add_argument("--summarise-only", action="store_true")
    a = ap.parse_args()
    if not a.summarise_only:
        for tag, val in GRID:
            build_and_eval(tag, val, a.viscosity)
    summarise(a.viscosity)
