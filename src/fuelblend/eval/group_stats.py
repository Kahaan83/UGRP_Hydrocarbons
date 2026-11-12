"""v5 BASE: shared group statistics, plus the method-family classifiers moved from the S4/S4b scripts (bodies unchanged).

Error convention everywhere: predicted - measured.  One fuel per double-count group is the caller's job; these functions treat the
`groups` argument as the resampling / hold-out unit.
"""
from __future__ import annotations

import re
from typing import Callable, Optional

import numpy as np
import pandas as pd

BOOT_B, BOOT_SEED, MIN_GROUPS = 5000, 0, 8
HEAT = "heat_of_combustion_net_MJ_kg"

# ---- moved from scripts/v5/s4b_heat_measured.py (class and evidence per normalised method code) ----
MEAS, CALC, UNC = "MEASURED", "CALCULATED", "UNCLASSIFIED"

# Class and evidence per normalised method code. Evidence is documentation in the repo or AJFTD page metadata
# on disk. No class was decided from agreement with any model.
_AUD = "docs/audit_fable_v5_2026-11-11.md finding B3"
_NOPAGE = ("saved AJFTD pages extras/fuelblend_data/debug/detail_*.html carry only the code in the 'Test Method' "
           "column (no legend or tooltip); extras/fuelblend_data/data/raw/ajftd/*.csv add no description; "
           "docs/audit_fable_v5_2026-11-11.md line 112 lists the code as unverified")
CLASS_EVIDENCE = {
    "D3338": (CALC, f"{_AUD}: 'D3338 is an estimate from density, aromatics and distillation'"),
    "D4809": (MEAS, f"{_AUD} and its next-step 3 (line 96): 'D4809/D240 (bomb calorimeter)'"),
    "D240": (MEAS, f"{_AUD} and its next-step 3 (line 96): 'D4809/D240 (bomb calorimeter)'"),
    "D4529": (CALC, f"{_AUD} / S4b README name D4529 as calculated; no rows with this code exist in "
                    "ajftd_fuels_long.csv heat rows (raw extras scrape lists it under the sulfur-free-basis metric only)"),
    "400C": (UNC, f"{_NOPAGE}; audit B3: 'meaning of the code not verified'"),
    "400D": (UNC, f"{_NOPAGE}; audit B3: 'meaning of the code not verified'"),
    "400E": (UNC, f"{_NOPAGE}; audit B3: 'meaning of the code not verified'"),
    "400IPA": (UNC, f"{_NOPAGE}; no documentation found"),
    "400IPB": (UNC, f"{_NOPAGE}; no documentation found"),
    "231A": (UNC, f"{_NOPAGE}; audit B3 groups '231A/B' without a meaning"),
    "231B": (UNC, f"{_NOPAGE}; audit B3 groups '231A/B' without a meaning"),
    "NA": (UNC, "method field empty in ajftd_fuels_long.csv"),
}


def norm_code(m) -> str:
    """'D 4809' -> 'D4809'; 'D4809 (SwRI[32])' -> 'D4809'; NaN -> 'NA'."""
    if m is None or (isinstance(m, float) and np.isnan(m)) or pd.isna(m):
        return "NA"
    return re.sub(r"\s+", "", str(m)).split("(")[0] or "NA"


def classify_method(m) -> str:
    return CLASS_EVIDENCE.get(norm_code(m), (UNC, ""))[0]

# ---- moved from scripts/v5/eval_bulk_descriptor.py ----
def flash_family(m: str) -> str:
    u = str(m).upper().replace(" ", "")
    if "93" in u or "IP170" in u:
        return "D93_or_IP170"
    if "D56" in u:
        return "D56"
    if u.startswith("220"):
        return "spec_table_220x"
    return "other"


# ---- statistics ----
def paired_group_bootstrap(abs_err_ref, abs_err_model, groups, n_boot: int = BOOT_B, seed: int = BOOT_SEED) -> tuple[float, float, float]:
    """Paired bootstrap over groups of pooled MAE(ref) - MAE(model); positive = model better.  -> (diff, lo, hi), 95 % percentile."""
    d = np.asarray(abs_err_ref, float) - np.asarray(abs_err_model, float)
    g = pd.Series(np.asarray(groups)).astype(str).to_numpy()
    ug, inv = np.unique(g, return_inverse=True)
    s = np.bincount(inv, weights=d)
    n = np.bincount(inv).astype(float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(ug), size=(n_boot, len(ug)))
    stat = s[idx].sum(1) / n[idx].sum(1)
    return float(d.mean()), float(np.percentile(stat, 2.5)), float(np.percentile(stat, 97.5))


def logo_constant(y, groups, strata=None) -> np.ndarray:
    """Leave-one-group-out mean of y for every row (rows of the scored group never enter), within `strata` when given.
    A stratum with no training row falls back to the leave-one-group-out mean over all strata."""
    y = np.asarray(y, float)
    g = np.asarray(groups).astype(str)
    st = np.zeros(len(y), dtype=object) if strata is None else np.asarray(strata).astype(object)
    out = np.full(len(y), np.nan)
    fin = np.isfinite(y)
    for grp in np.unique(g):
        te = g == grp
        tr_all = fin & ~te
        glob = y[tr_all].mean() if tr_all.any() else np.nan
        for s in np.unique(st[te]):
            m = te & (st == s)
            tr = tr_all & (st == s)
            out[m] = y[tr].mean() if tr.any() else glob
    return out


def verdict(n_groups: int, lo: float, hi: float, min_groups: int = MIN_GROUPS) -> str:
    if n_groups < min_groups:
        return "n_too_small"
    if lo > 0:
        return "beats_trivial"
    if hi < 0:
        return "worse_than_trivial"
    return "no_skill_shown"


def cell_metrics(err, groups, unsafe: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> dict:
    """Pooled error summary of finite errors (pred - meas).  unsafe(err) -> boolean array; unsafe_rate is its mean."""
    e = np.asarray(err, float)
    g = np.asarray(groups).astype(str)
    ok = np.isfinite(e)
    e, g = e[ok], g[ok]
    return dict(n_groups=int(len(np.unique(g))), n_rows=int(e.size),
                bias=float(e.mean()) if e.size else np.nan, mae=float(np.abs(e).mean()) if e.size else np.nan,
                rmse=float(np.sqrt((e ** 2).mean())) if e.size else np.nan,
                unsafe_rate=float(np.mean(unsafe(e))) if (unsafe is not None and e.size) else np.nan)
