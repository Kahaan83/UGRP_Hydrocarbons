"""v5 BASE: build scoreboard rows (schema in results/v5/BASE/README.md) from error arrays or from pre-aggregated numbers."""
from __future__ import annotations

import numpy as np
import pandas as pd

from fuelblend.eval.group_stats import cell_metrics, paired_group_bootstrap, verdict

SCHEMA = ("property tier set model variant status n_groups n_rows bias mae rmse unsafe_rate trivial_model mae_trivial "
          "diff_trivial_minus_model ci_lo ci_hi verdict margin90 coverage90 n_excluded_censored n_excluded_nonfinite "
          "n_low_coverage evidence_path label").split()


def make_row(prop, tier, set_, model, variant, status, label, evidence, *, err=None, groups=None, unsafe=None,
             err_trivial=None, trivial_model="", margin90=np.nan, coverage90=np.nan, n_excluded_censored=0,
             n_excluded_nonfinite=None, n_low_coverage=0, fixed=None, extra=None) -> dict:
    """err / err_trivial: pred - measured (same length as groups).  Non-finite model errors are excluded and counted.
    fixed: pre-aggregated values (n_groups, n_rows, bias, mae, rmse, unsafe_rate) used when no per-row errors exist."""
    row = dict(property=prop, tier=tier, set=set_, model=model, variant=variant, status=status, n_groups=np.nan, n_rows=0,
               bias=np.nan, mae=np.nan, rmse=np.nan, unsafe_rate=np.nan, trivial_model=trivial_model, mae_trivial=np.nan,
               diff_trivial_minus_model=np.nan, ci_lo=np.nan, ci_hi=np.nan, verdict="not_applicable", margin90=margin90,
               coverage90=coverage90, n_excluded_censored=n_excluded_censored, n_excluded_nonfinite=0 if n_excluded_nonfinite is None else n_excluded_nonfinite,
               n_low_coverage=n_low_coverage, evidence_path=evidence, label=label)
    if err is not None:
        e, g = np.asarray(err, float), np.asarray(groups).astype(str)
        fin = np.isfinite(e)
        if n_excluded_nonfinite is None:
            row["n_excluded_nonfinite"] = int((~fin).sum())
        if err_trivial is not None:
            et = np.asarray(err_trivial, float)
            fin = fin & np.isfinite(et)
            et = et[fin]
        e, g = e[fin], g[fin]
        row.update(cell_metrics(e, g, unsafe))
        if err_trivial is not None and e.size:
            d, lo, hi = paired_group_bootstrap(np.abs(et), np.abs(e), g)
            row.update(mae_trivial=float(np.abs(et).mean()), diff_trivial_minus_model=d, ci_lo=lo, ci_hi=hi,
                       verdict=verdict(row["n_groups"], lo, hi))
    if fixed:
        row.update(fixed)
    if row["status"] != "no_evaluation" and np.isfinite(row["n_groups"]) and row["n_groups"] < 8 and row["verdict"] != "not_applicable":
        row["verdict"] = "n_too_small"
    if extra:
        row.update(extra)
    return row


def no_eval(prop, tier, set_, reason, *, model="-", variant="designated", evidence="", label="descriptive") -> dict:
    r = make_row(prop, tier, set_, model, variant, "no_evaluation", label, evidence)
    r["reason"] = reason
    return r


def frame(rows) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    extra = [c for c in df.columns if c not in SCHEMA]
    return df[SCHEMA + extra]
