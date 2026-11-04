#!/usr/bin/env python
"""T3: prediction error stratified by mixture complexity (n_comp buckets).

    python scripts/t3_complexity.py

Evaluates the WP3 physics baselines from data/runs/baseline_preds.csv
(conflict=False rows only; clustered 2000-resample bootstrap by system_id) and
writes reports/T3_COMPLEXITY.md. Flash points compare in absolute degC; relative
error is meaningless below zero C, so flash gets MAE/RMSE only.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from fuelblend.eval import metrics as evm  # noqa: E402

PRED = BASE / "data/runs/baseline_preds.csv"
OUT = BASE / "reports/T3_COMPLEXITY.md"

VERIFY = {  # property -> {baseline: error mode}
    "viscosity_dyn_mPas": {"refutas": "ln", "gn_mass": "ln"},
    "density_kg_m3": {"dens_idvol": "lin"},
    "flash_point_C": {"chat_nominal": "lin"},
}


def bucket(n: int) -> str:
    return "pure" if n <= 1 else "binary" if n == 2 else "ternary" if n == 3 else "4comp" if n == 4 else "5+comp"


def rows_for(df: pd.DataFrame, prop: str, base: str) -> pd.DataFrame:
    d = df[(df["property"] == prop) & df[base].notna() & (df["conflict"] == False)]  # noqa: E712
    d = d[d["value"].notna()]
    if prop == "viscosity_dyn_mPas":
        d = d[(d["value"] > 0) & (d[base] > 0)]
    return d.assign(bucket=d["n_comp"].map(bucket))


def err(y: np.ndarray, yhat: np.ndarray, mode: str):
    if mode == "ln":
        e = np.log(np.clip(y, 1e-9, None)) - np.log(np.clip(yhat, 1e-9, None))
        return dict(mae=float(np.mean(np.abs(e))), rmse=float(np.sqrt(np.mean(e ** 2))), mape=float(np.mean(np.abs(e)) * 100), neg=len(e[e < 0]))
    e = y - yhat
    return dict(mae=float(np.mean(np.abs(e))), rmse=float(np.sqrt(np.mean(e ** 2))), mape=float("nan"), neg=len(e[e < 0]))


def main() -> None:
    df = pd.read_csv(PRED, low_memory=False)
    lines = ["# T3 - baseline error by mixture complexity", "",
             "baselines from data/runs/baseline_preds.csv; conflict rows excluded; cluster bootstrap (2000, system_id).",
             "Viscosity numbers are log-RMSE (ln mPa.s units); density and flash in kg/m3 and degC. Flash is a",
             "temperature whose errors are signed, so bias matters more than RMSE below the table.", ""]

    for prop, bases in VERIFY.items():
        lines += [f"## {prop}", ""]
        for base, mode in bases.items():
            d = rows_for(df, prop, base)
            if not len(d):
                lines += [f"- {base}: no applicable rows", ""]
                continue
            lines.append(f"### baseline `{base}` (n_applicable={len(d)})")
            lines += ["", "| bucket | n | MAE | RMSE | rel(err) | mean signed | 95% CI rmse |",
                      "|---|---|---|---|---|---|---|"]
            for b in ("pure", "binary", "ternary", "4comp", "5+comp"):
                dd = d[d["bucket"] == b]
                if not len(dd):
                    lines.append(f"| {b} | 0 | - | - | - | - | - |")
                    continue
                y, yhat = dd["value"].to_numpy(float), dd[base].to_numpy(float)
                m = err(y, yhat, mode)
                ci = evm.bootstrap_ci(y, yhat, dd["system_id"].to_numpy(),
                                      metric_fn=evm.rmse, n_boot=2000, seed=20261030) if len(dd) >= 10 else {"lo": float("nan"), "hi": float("nan")}
                lines.append(f"| {b} | {len(dd)} | {m['mae']:.4f} | {m['rmse']:.4f} | "
                             f"{(m['mape'] if np.isfinite(m['mape']) else (m['mae'] / max(np.abs(y).mean(), 1e-9) * 100)):.1f}% | "
                             f"{np.mean(yhat - y):+.3f} | [{ci['lo']:.3f}, {ci['hi']:.3f}] |")
            lines.append("")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
