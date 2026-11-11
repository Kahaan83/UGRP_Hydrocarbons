"""Flash-point effect of the 2026-11-09 Antoine reference gate (refit_antoine.py, Ruzicka & Majer Cox curve).

Zero-shot Le Chatelier + Raoult flash (same solver as scripts/v3/audit_flash_anchoring.py) on:
  * the 126 Leenhouts hydrocarbon mixtures (results/eval/flash_loso_predictions.csv)
  * pure-compound flash records of the n-alkanes whose Antoine set changed (data/unified/measurements_dedup.csv)
with three registries:
  before   = extras/fuelblend_data/archive/2026-11-09_pre_cox_gate/compounds_root.csv
  after    = data/unified/compounds.csv
  anchored = data/processed/v2/compounds_anchored_vp.csv   (behind CLM-12; not changed by the gate)

Writes results/v4/antoine_cox_gate/{flash_before_after.csv, summary.md}. Nothing in results/v2|v3 is touched.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "v4" / "antoine_cox_gate"
P_TOTAL_KPA = 101.325
REGISTRIES = {
    "before": ROOT / "extras/fuelblend_data/archive/2026-11-09_pre_cox_gate/compounds_root.csv",
    "after": ROOT / "data/unified/compounds.csv",
    "anchored": ROOT / "data/processed/v2/compounds_anchored_vp.csv",
}


def solve_flash_point(x_mole, antoine_list, lfl_list):
    """Copy of scripts/v3/audit_flash_anchoring.solve_flash_point (LFL temperature-corrected), deg C."""
    def haz(T_K):
        fac = max(1.0 - 0.000721 * (T_K - 273.15 - 25.0), 0.05)
        h = 0.0
        for xi, (A, B, C), lfl in zip(x_mole, antoine_list, lfl_list):
            d = T_K + C
            if d <= 1.0:
                return np.nan
            h += xi * 10.0 ** (A - B / d) / (P_TOTAL_KPA * lfl * fac)
        return h - 1.0
    try:
        return float(brentq(haz, 150.0, 600.0) - 273.15)
    except Exception:
        return np.nan


def predict(reg: pd.DataFrame, keys: list[str], x: np.ndarray) -> float:
    if not all(k in reg.index for k in keys):
        return np.nan
    r = reg.loc[keys]
    abc = list(zip(r.antoine_A.astype(float), r.antoine_B.astype(float), r.antoine_C.astype(float)))
    if not np.isfinite(np.asarray(abc)).all():
        return np.nan
    lfl = [v / 100.0 if pd.notna(v) else 0.01 for v in r.lfl_vol_pct]   # same 1 % fallback as the R2 audit
    return solve_flash_point(x, abc, lfl)


def metrics(err: pd.Series) -> dict:
    e = err.dropna()
    return dict(n=len(e), MAE=round(e.abs().mean(), 2), bias=round(e.mean(), 2), RMSE=round(float(np.sqrt((e ** 2).mean())), 2),
                over_5C_pct=round(100 * (e > 5).mean(), 1), max_over=round(e.max(), 2))


def main() -> None:
    sys.stdout.reconfigure(encoding="utf-8")          # Windows console default (cp1252) cannot print "−"
    OUT.mkdir(parents=True, exist_ok=True)
    regs = {n: pd.read_csv(p).set_index("inchikey") for n, p in REGISTRIES.items()}
    changed = [k for k in regs["after"].index.intersection(regs["before"].index)
               if not np.allclose(regs["after"].loc[k, ["antoine_A", "antoine_B", "antoine_C"]].astype(float),
                                  regs["before"].loc[k, ["antoine_A", "antoine_B", "antoine_C"]].astype(float), equal_nan=True)]

    mix = pd.read_csv(ROOT / "results/eval/flash_loso_predictions.csv")
    assert len(mix) == 126
    rows = []
    for r in mix.itertuples():
        keys, x = json.loads(r.comp_inchikeys), np.asarray(json.loads(r.x_mole), float)
        rows.append(dict(set="mixture_126", system_id=r.system_id, comp_inchikeys=r.comp_inchikeys, measured_C=r.value,
                         touches_changed=any(k in changed for k in keys),
                         **{f"pred_{n}": predict(reg, keys, x) for n, reg in regs.items()}))

    dedup = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", low_memory=False)
    pure = dedup[(dedup.property == "flash_point_C") & (dedup.n_comp == 1)]
    for r in pure.itertuples():
        keys = json.loads(r.comp_inchikeys)
        if keys[0] in changed:
            rows.append(dict(set="pure_changed", system_id=r.system_id, comp_inchikeys=r.comp_inchikeys, measured_C=r.value,
                             touches_changed=True, **{f"pred_{n}": predict(reg, keys, np.ones(1)) for n, reg in regs.items()}))

    df = pd.DataFrame(rows)
    for n in regs:
        df[f"err_{n}"] = df[f"pred_{n}"] - df.measured_C
    df.to_csv(OUT / "flash_before_after.csv", index=False)

    lines = ["# Antoine reference gate: flash-point effect (2026-11-09)", "",
             f"Antoine sets changed by the gate: {len(changed)} "
             f"({', '.join(regs['after'].loc[changed, 'smiles'])}).",
             "Error = predicted − measured (°C); positive = over-prediction (the unsafe direction). LFL fallback 1 % as in the R2 audit.", ""]
    for label, sub in (("126 mixtures, all", df[df.set == "mixture_126"]),
                       ("126 mixtures, only those containing a changed compound", df[(df.set == "mixture_126") & df.touches_changed]),
                       ("pure flash records of the changed compounds", df[df.set == "pure_changed"])):
        t = pd.DataFrame({n: metrics(sub[f"err_{n}"]) for n in regs}).T
        md = ["| registry | " + " | ".join(t.columns) + " |", "|---" * (len(t.columns) + 1) + "|"]
        md += [f"| {i} | " + " | ".join(str(v) for v in row) + " |" for i, row in t.iterrows()]
        lines += [f"## {label}", "", *md, ""]
    (OUT / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
