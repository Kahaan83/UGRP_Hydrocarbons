"""v5 delta_smoke_cert Task 1.4: real-fuel smoke margin (R-SC-2). Readings declared in
results/v5/BASE/delta_smoke_cert/README.md before any fit; do not change them here.

Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_smoke_cert.py
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.baselines.sooting import BULK_FEATURES, fit_bulk, predict_bulk, smoke_point_from_bulk
from fuelblend.data.loaders_ajftd_v5 import load_ajftd_smoke_naph
from fuelblend.design.smoke_margin_v5 import real_fuel_class_ranges

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "v5" / "BASE" / "delta_smoke_cert"
AROM_LO, AROM_HI = 8.0, 23.0
MIN_GROUP, MIN_FUELS, MIN_GROUPS = 40, 200, 5
LAB_MARGIN_MM = 7.8


def conformal_over(res, alpha: float) -> float:
    """Split-conformal (1-alpha) upper quantile: ceil((n+1)(1-alpha))-th smallest residual (inf if it exceeds n)."""
    r = np.sort(np.asarray(res, float))
    n = len(r)
    k = math.ceil((n + 1) * (1 - alpha))
    return math.inf if (n == 0 or k > n) else float(r[k - 1])


FEATS = (BULK_FEATURES[0],)  # aromatics only: diaromatics is absent from the AJFTD smoke table


def _value_df(d: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({"value": d["smoke_point_mm"].to_numpy(),
                         BULK_FEATURES[0]: d["aromatics_vol_pct"].to_numpy()})


def fit_smoke_margin(df: pd.DataFrame, alpha: float = 0.10, group_col: str = "source_group") -> dict:
    """Smoke margin (mm) from leave-one-group-out refits of the production bulk model form (fit_bulk).

    Rows with a smoke qualifier (censored), no smoke value, or aromatics outside 8-23 vol % are dropped
    (censored rows counted). Residual = predicted - measured (positive = over-prediction = unsafe), predicted from
    aromatics only (diaromatics NaN, as in the design call). Margin = max(fuel-level conformal quantile, largest
    per-group conformal quantile over groups with >= 40 fuels); inf unless >= 200 fuels sit in >= 5 such groups."""
    d = df.copy()
    if "smoke_qualifier" in d:
        q = d["smoke_qualifier"].fillna("").astype(str).str.strip()
    else:
        q = pd.Series("", index=d.index)
    n_cens = int((q != "").sum())
    keep = (q == "") & d["smoke_point_mm"].notna() & d["aromatics_vol_pct"].between(AROM_LO, AROM_HI)
    d = d[keep].reset_index(drop=True)
    res = np.full(len(d), np.nan)
    for g in d[group_col].unique():
        te = (d[group_col] == g).to_numpy()
        tr = d[~te]
        if len(tr) < 3:
            continue
        coef, means = fit_bulk(_value_df(tr), features=FEATS)
        res[te] = predict_bulk(_value_df(d[te]), coef, means, features=FEATS) - d.loc[te, "smoke_point_mm"].to_numpy()
    ok = ~np.isnan(res)
    res_ok, d_ok = res[ok], d[ok]
    sizes = d_ok.groupby(group_col).size()
    per_n = {str(g): int(n) for g, n in sizes.items()}
    per_m = {str(g): conformal_over(res_ok[(d_ok[group_col] == g).to_numpy()], alpha) for g in sizes.index}
    elig = [g for g, n in per_n.items() if n >= MIN_GROUP]
    n_in_elig = int(sum(per_n[g] for g in elig))
    fuel_m = conformal_over(res_ok, alpha)
    passed = n_in_elig >= MIN_FUELS and len(elig) >= MIN_GROUPS
    margin = max([fuel_m] + [per_m[g] for g in elig]) if passed else math.inf
    cov = float(np.mean(res_ok <= margin)) if len(res_ok) else float("nan")
    return dict(margin=margin, fuel_level_margin=fuel_m, per_group_margin=per_m, per_group_n=per_n,
                eligible_groups=elig, n_groups=len(elig), n_groups_all=len(per_n), n_rows=int(len(res_ok)),
                n_rows_in_eligible_groups=n_in_elig, direction="over", coverage_oof=cov, alpha=alpha,
                gate_passed=bool(passed), n_censored_excluded=n_cens, oof_residuals=res_ok.tolist())


def _stored_reference(d: pd.DataFrame) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pred = np.array([smoke_point_from_bulk(a) for a in d["aromatics_vol_pct"]])
    return pred - d["smoke_point_mm"].to_numpy()


def production_predictor_margin(d: pd.DataFrame, alpha: float = 0.10, group_col: str = "source_group") -> dict:
    """R-SC-2c: residuals of smoke_point_from_bulk(aromatics, None) on in-domain fuels not in its 100 training fuels."""
    from fuelblend.baselines.sooting import DATA
    train = pd.read_csv(DATA / "smoke_point_ajftd.csv")
    q = d["smoke_qualifier"].fillna("").astype(str)
    ind = d[q.eq("") & d["smoke_point_mm"].notna() & d["aromatics_vol_pct"].between(AROM_LO, AROM_HI)]
    n_match = int(ind["fuel_id"].isin(train["fuel_id"]).sum())
    held = ind[~ind["fuel_id"].isin(train["fuel_id"])]
    r = _stored_reference(held)
    per_n = {str(g): int((held[group_col] == g).sum()) for g in held[group_col].unique()}
    per_m = {g: conformal_over(r[(held[group_col] == g).to_numpy()], alpha) for g in per_n}
    elig = [g for g, n in per_n.items() if n >= MIN_GROUP]
    fuel_m = conformal_over(r, alpha)
    n_el = sum(per_n[g] for g in elig)
    passed = n_el >= MIN_FUELS and len(elig) >= MIN_GROUPS
    margin = max([fuel_m] + [per_m[g] for g in elig]) if passed else math.inf
    # refit-OLS (full in-domain data, aromatics only) vs stored predictor over aromatics grid
    coef, means = fit_bulk(_value_df(ind), features=FEATS)
    grid = np.arange(8.0, 23.0001, 0.5)
    ols = predict_bulk(pd.DataFrame({BULK_FEATURES[0]: grid}), coef, means, features=FEATS)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stored = np.array([smoke_point_from_bulk(float(a)) for a in grid])
    return dict(margin=margin, fuel_level_margin=fuel_m, per_group_margin=per_m, per_group_n=per_n,
                eligible_groups=elig, n_in_domain=int(len(ind)), n_matched_training_fuels=n_match,
                n_scored=int(len(held)), gate_passed=bool(passed),
                max_abs_pred_diff_refit_vs_stored_mm=float(np.max(np.abs(ols - stored))))


def _json_default(o):
    return str(o)


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = load_ajftd_smoke_naph()
    qual = raw["smoke_qualifier"].fillna("").astype(str)
    d = raw[raw["smoke_point_mm"].notna() | qual.ne("")]
    main_fit = fit_smoke_margin(d)
    sens = fit_smoke_margin(d[d["smoke_point_mm"] != 25.0])  # sensitivity only
    ind = d[d["smoke_qualifier"].fillna("").astype(str).eq("") & d["smoke_point_mm"].notna()
            & d["aromatics_vol_pct"].between(AROM_LO, AROM_HI)]
    leak = _stored_reference(ind)
    pg = {str(g): conformal_over(leak[(ind["source_group"] == g).to_numpy()], 0.10)
          for g in ind["source_group"].unique() if (ind["source_group"] == g).sum() >= MIN_GROUP}
    leaky = dict(fuel_level_margin=conformal_over(leak, 0.10), n=int(len(leak)), per_group=pg)
    leaky["certifying"] = max([leaky["fuel_level_margin"]] + list(pg.values()))
    res = np.asarray(main_fit["oof_residuals"])
    prod = production_predictor_margin(d)
    final = max(main_fit["margin"], prod["margin"])
    rows = [dict(source_group=g, n_fuels=n, eligible=g in main_fit["eligible_groups"],
                 margin_mm=main_fit["per_group_margin"][g])
            for g, n in sorted(main_fit["per_group_n"].items(), key=lambda kv: -kv[1])]
    rows.append(dict(source_group="ALL (fuel-level)", n_fuels=main_fit["n_rows"], eligible=True,
                     margin_mm=main_fit["fuel_level_margin"]))
    rows.append(dict(source_group="CERTIFYING (max)", n_fuels=main_fit["n_rows"], eligible=main_fit["gate_passed"],
                     margin_mm=main_fit["margin"]))
    pd.DataFrame(rows).to_csv(OUT / "margin_fit.csv", index=False)
    r = real_fuel_class_ranges()
    # max margin that certifies = smoke_point_from_bulk(aromatics) - 18 mm (D1655 alternative limit), derived not hard-coded
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cert_levels = {a: round(float(smoke_point_from_bulk(float(a))) - 18.0, 2) for a in (8, 10, 12, 15, 20)}
    m = main_fit["margin"]
    read = dict(
        task="1.4", reading="R-SC-2", alpha=0.10, direction="over", gate_passed=main_fit["gate_passed"],
        gate_rule=">=200 in-domain fuels (counted in groups with >=40) and >=5 groups with >=40 fuels",
        certifying_margin_mm=m, fuel_level_margin_mm=main_fit["fuel_level_margin"],
        per_group_margin_mm=main_fit["per_group_margin"], per_group_n=main_fit["per_group_n"],
        eligible_groups=main_fit["eligible_groups"], n_in_domain_fuels=main_fit["n_rows"],
        n_in_eligible_groups=main_fit["n_rows_in_eligible_groups"], n_groups_all=main_fit["n_groups_all"],
        n_censored_excluded=main_fit["n_censored_excluded"], coverage_oof=main_fit["coverage_oof"],
        oof_residual_mean=float(res.mean()), oof_residual_median=float(np.median(res)),
        oof_residual_q90=float(np.quantile(res, 0.9)),
        sensitivity_smoke25_removed=dict(margin=sens["margin"], fuel_level=sens["fuel_level_margin"],
                                         n_rows=sens["n_rows"], per_group=sens["per_group_margin"]),
        leaky_stored_coefficient_reference=leaky,
        lab_margin_mm=LAB_MARGIN_MM, max_margin_that_certifies_mm_by_arom_vol_pct=cert_levels,
        margin_vs_lab_ratio=(m / LAB_MARGIN_MM) if math.isfinite(m) else None,
        certifies_at={str(a): bool(m < v) for a, v in cert_levels.items()},
        class_mix_range_mass_frac={k: list(v) for k, v in r.items()},
        r_sc_2c=prod, certifying_margin_final=final,
        user_decision_required="R-SC-2b: yaml key not added; needs a quoted user answer")
    (OUT / "readings.json").write_text(json.dumps(read, indent=2, default=_json_default))
    print(json.dumps({k: v for k, v in read.items() if k != "per_group_n"}, indent=1, default=_json_default))
    return read


if __name__ == "__main__":
    main()
