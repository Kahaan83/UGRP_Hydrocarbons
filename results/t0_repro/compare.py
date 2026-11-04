"""T0.1 comparison: stored baseline outputs vs rerun outputs (and PROGRESS_REPORT quotes vs rerun).

Reads only results/t0_repro/reference_backup/ (pre-run), results/t0_repro/rerun_outputs/ (post-run) and
hard-coded numbers transcribed from PROGRESS_REPORT.md. Writes results/t0_repro/comparison.csv.
"""
import numpy as np
import pandas as pd

REF = "results/t0_repro/reference_backup"
NEW = "results/t0_repro/rerun_outputs"
ROWS = []


def add(script, metric, old, new, tol, verdict=None, note=""):
    old_f, new_f = float(old), float(new)
    diff = abs(new_f - old_f)
    if verdict is None:
        verdict = "pass" if diff <= tol else "fail"
    ROWS.append(dict(script=script, metric=metric, old_value=old_f, new_value=new_f, abs_diff=diff,
                     tolerance=tol, verdict=verdict if not note else f"{verdict} ({note})"))


def family_tol(metric):
    """Numeric tolerance by unit of the metric (ground rule: viscosity 0.002 log, density 0.15, flash 0.3 C)."""
    m = metric.lower()
    if "n_obs" in m or "n_systems" in m or "n_rows" in m or "n_fuels" in m or "n_comp" in m:
        return 0.0
    if "over5" in m:
        return 0.05          # the report quotes these to 1 dp
    if "p(model" in m or "p(hgb" in m:
        return 0.02
    if "mape" in m:
        return 0.05
    if "dens" in m or "density" in m:
        return 0.15
    if "visc" in m or "kv_" in m or "rmse_ln" in m:
        return 0.002
    return 0.3          # flash MAE/RMSE/delta/CI


# ---------------------------------------------------------------- grouped 5-fold CV
old = pd.read_csv(f"{REF}/eval/grouped_cv_summary.csv").set_index("property")
new = pd.read_csv(f"{NEW}/eval/grouped_cv_summary.csv").set_index("property")
for prop in old.index:
    for col in old.columns:
        if col == "property" or col not in new.columns:
            continue
        if pd.isna(old.loc[prop, col]):
            continue
        add("grouped_cv_summary.csv", f"{prop}.{col}", old.loc[prop, col], new.loc[prop, col],
            family_tol(("dens" if prop == "dens" else "visc" if prop == "visc" else "flash") + " " + col))

# ---------------------------------------------------------------- flash LOSO
old = pd.read_csv(f"{REF}/eval/flash_loso_summary.csv").iloc[0]
new = pd.read_csv(f"{NEW}/eval/flash_loso_summary.csv").iloc[0]
for col in old.index:
    if col == "property" or pd.isna(old[col]):
        continue
    add("flash_loso_summary.csv", f"flash.{col}", old[col], new[col], family_tol(f"flash {col}"))
fl_new = pd.read_csv(f"{NEW}/eval/flash_loso_summary.csv").iloc[0]
gc_new = pd.read_csv(f"{NEW}/eval/grouped_cv_summary.csv").set_index("property")

# per-row predictions of flash LOSO: identical file?
op = pd.read_csv(f"{REF}/eval/flash_loso_predictions.csv")
np_ = pd.read_csv(f"{NEW}/eval/flash_loso_predictions.csv")
d = (np_["hgb_flash"] - op["hgb_flash"]).abs().max()
add("flash_loso_predictions.csv", "max|delta hgb_flash| (per row)", 0.0, d, 1e-6)

# ---------------------------------------------------------------- hindcast
def hindcast_summary(path):
    res = pd.read_csv(path)
    rows = []
    for (route, kind, T), g in res.groupby(["route", "kind", "T_C"]):
        row = dict(route=route, kind=kind, T_C=T, n_fuels=len(g))
        for mdl in ["physics", "extrap_deepsets", "extrap_pairwise"]:
            e = (g[mdl] - g["measured"]) if kind == "density" else np.log(g[mdl] / g["measured"])
            row[f"{mdl}_rmse"] = float(np.sqrt((e ** 2).mean()))
            row[f"{mdl}_bias"] = float(e.mean())
        rows.append(row)
    return pd.DataFrame(rows)


oh = pd.read_csv(f"{REF}/hindcast/hindcast_predictions.csv")
nh = pd.read_csv(f"{NEW}/hindcast/hindcast_predictions.csv")
for mdl in ["physics", "extrap_deepsets", "extrap_pairwise"]:
    add("hindcast_predictions.csv", f"max|delta {mdl}| (per row)", 0.0,
        (nh[mdl] - oh[mdl]).abs().max(), 0.01, note="float32 write precision")
ohs, nhs = hindcast_summary(f"{REF}/hindcast/hindcast_predictions.csv"), \
           hindcast_summary(f"{NEW}/hindcast/hindcast_predictions.csv")
ohs, nhs = ohs.set_index(["route", "kind", "T_C"]), nhs.set_index(["route", "kind", "T_C"])
for key in ohs.index:
    for col in ohs.columns:
        add("hindcast_predictions.csv (derived)", f"{key[1]}@{key[2]:.0f}C.{col}", ohs.loc[key, col],
            nhs.loc[key, col], family_tol(col))
# not a reported metric, but a visible difference: Monte-Carlo spreads
for mdl in ["extrap_deepsets", "extrap_pairwise", "physics"]:
    add("hindcast_predictions.csv", f"mean {mdl}_mc_spread", oh[f"{mdl}_mc_spread"].mean(),
        nh[f"{mdl}_mc_spread"].mean(), 1e-3, verdict="deviation-not-reported",
        note="MC-Dirichlet spread depends on --mc / draw sequence; not quoted in the report")

# ---------------------------------------------------------------- PROGRESS_REPORT vs rerun
REPORT = [
    # (section, metric, reported value, rerun value, tolerance)
    ("Â§3.3 flash LOSO", "base_mae", 6.16, fl_new["base_mae"], 0.3),
    ("Â§3.3 flash LOSO", "base_rmse", 10.24, fl_new["base_rmse"], 0.3),
    ("Â§3.3 flash LOSO", "base_over5_pct", 19.8, fl_new["base_over5_pct"], 0.05),
    ("Â§3.3 flash LOSO", "base_max_over", 42.4, fl_new["base_max_over"], 0.3),
    ("Â§3.3 flash LOSO", "hgb_mae", 2.75, fl_new["hgb_mae"], 0.3),
    ("Â§3.3 flash LOSO", "hgb_rmse", 3.88, fl_new["hgb_rmse"], 0.3),
    ("Â§3.3 flash LOSO", "hgb_over5_pct", 9.5, fl_new["hgb_over5_pct"], 0.05),
    ("Â§3.3 flash LOSO", "hgb_max_over", 11.4, fl_new["hgb_max_over"], 0.3),
    ("Â§3.3 flash LOSO", "delta mae", -3.43, fl_new["hgb-base_mae"], 0.3),
    ("Â§3.3 flash LOSO", "ci_lo", -6.03, fl_new["ci_95_lo"], 0.3),
    ("Â§3.3 flash LOSO", "ci_hi", -1.26, fl_new["ci_95_hi"], 0.3),
    ("Â§3.3 flash LOSO", "P(model<base)", 1.000, fl_new["P(hgb<base)"], 0.02),
    ("Â§3.4 grouped CV visc", "base log-RMSE", 0.0699, gc_new.loc["visc", "base_rmse"], 0.002),
    ("Â§3.4 grouped CV visc", "base MAPE %", 5.43, gc_new.loc["visc", "base_mape_pct"], 0.05),
    ("Â§3.4 grouped CV visc", "model log-RMSE", 0.0438, gc_new.loc["visc", "model_rmse"], 0.002),
    ("Â§3.4 grouped CV visc", "model MAPE %", 2.48, gc_new.loc["visc", "model_mape_pct"], 0.05),
    ("Â§3.4 grouped CV visc", "delta", -0.0263, gc_new.loc["visc", "diff_metric"], 0.002),
    ("Â§3.4 grouped CV visc", "ci_lo", -0.0331, gc_new.loc["visc", "ci_95_lo"], 0.002),
    ("Â§3.4 grouped CV visc", "ci_hi", -0.0194, gc_new.loc["visc", "ci_95_hi"], 0.002),
    ("Â§3.4 grouped CV dens", "base RMSE", 2.3859, gc_new.loc["dens", "base_rmse"], 0.15),
    ("Â§3.4 grouped CV dens", "base MAE", 1.32, gc_new.loc["dens", "base_mae"], 0.15),
    ("Â§3.4 grouped CV dens", "model RMSE", 1.9908, gc_new.loc["dens", "model_rmse"], 0.15),
    ("Â§3.4 grouped CV dens", "model MAE", 0.93, gc_new.loc["dens", "model_mae"], 0.15),
    ("Â§3.4 grouped CV dens", "delta", -0.4073, gc_new.loc["dens", "diff_metric"], 0.15),
    ("Â§3.4 grouped CV dens", "ci_lo", -0.6126, gc_new.loc["dens", "ci_95_lo"], 0.15),
    ("Â§3.4 grouped CV dens", "ci_hi", -0.2366, gc_new.loc["dens", "ci_95_hi"], 0.15),
    ("Â§3.4 grouped CV flash", "base MAE", 6.16, gc_new.loc["flash", "base_mae"], 0.3),
    ("Â§3.4 grouped CV flash", "base over5 %", 19.8, gc_new.loc["flash", "base_over5_pct"], 0.05),
    ("Â§3.4 grouped CV flash", "model MAE", 3.44, gc_new.loc["flash", "model_mae"], 0.3),
    ("Â§3.4 grouped CV flash", "model over5 %", 13.5, gc_new.loc["flash", "model_over5_pct"], 0.05),
    ("Â§3.4 grouped CV flash", "delta", -2.74, gc_new.loc["flash", "diff_metric"], 0.3),
    ("Â§3.4 grouped CV flash", "ci_lo", -5.78, gc_new.loc["flash", "ci_95_lo"], 0.3),
    ("Â§3.4 grouped CV flash", "ci_hi", -0.15, gc_new.loc["flash", "ci_95_hi"], 0.3),
    ("Â§3.4 grouped CV flash", "P(model<base)", 0.980, gc_new.loc["flash", "P(model<base)"], 0.02),
    ("Â§2.4 density baseline", "mixture density base RMSE", 2.3859, gc_new.loc["dens", "base_rmse"], 0.15),
    ("Â§3.2 hindcast", "visc -20C physics rmse_ln", 0.3189, float(nhs.loc[("bulk", "kv", -20.0), "physics_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc -20C deepsets rmse_ln", 0.1633, float(nhs.loc[("bulk", "kv", -20.0), "extrap_deepsets_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc -20C pairwise rmse_ln", 0.1795, float(nhs.loc[("bulk", "kv", -20.0), "extrap_pairwise_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc -40C physics rmse_ln", 0.6878, float(nhs.loc[("bulk", "kv", -40.0), "physics_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc -40C deepsets rmse_ln", 0.4968, float(nhs.loc[("bulk", "kv", -40.0), "extrap_deepsets_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc -40C pairwise rmse_ln", 0.5012, float(nhs.loc[("bulk", "kv", -40.0), "extrap_pairwise_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc 20C physics rmse_ln", 0.1058, float(nhs.loc[("bulk", "kv", 20.0), "physics_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc 20C deepsets rmse_ln", 0.1172, float(nhs.loc[("bulk", "kv", 20.0), "extrap_deepsets_rmse"]), 0.002),
    ("Â§3.2 hindcast", "visc 20C pairwise rmse_ln", 0.1082, float(nhs.loc[("bulk", "kv", 20.0), "extrap_pairwise_rmse"]), 0.002),
    ("Â§3.2 hindcast", "density 15C physics RMSE", 9.73, float(nhs.loc[("bulk", "density", 15.0), "physics_rmse"]), 0.15),
    ("Â§3.2 hindcast", "density 15C physics bias", 0.51, float(nhs.loc[("bulk", "density", 15.0), "physics_bias"]), 0.15),
    ("Â§3.2 hindcast", "density 15C deepsets RMSE", 10.30, float(nhs.loc[("bulk", "density", 15.0), "extrap_deepsets_rmse"]), 0.15),
    ("Â§3.2 hindcast", "density 15C deepsets bias", 3.64, float(nhs.loc[("bulk", "density", 15.0), "extrap_deepsets_bias"]), 0.15),
    ("Â§3.2 hindcast", "density 15C pairwise RMSE", 9.85, float(nhs.loc[("bulk", "density", 15.0), "extrap_pairwise_rmse"]), 0.15),
    ("Â§3.2 hindcast", "density 15C pairwise bias", 1.22, float(nhs.loc[("bulk", "density", 15.0), "extrap_pairwise_bias"]), 0.15),
]
for sec, met, rep, rerun, tol in REPORT:
    add(f"PROGRESS_REPORT {sec}", met, rep, rerun, tol)

# ---------------------------------------------------------------- baselines_eval / BASELINES.md
def mixture_summary(path):
    """Parse the 'Mixture-only summary' markdown table."""
    txt = open(path, encoding="utf-8").read()
    sec = txt.split("### Mixture-only summary", 1)[1].split("| property | baseline |", 1)[1]
    sec = sec.split("\n\n", 1)[0].strip().splitlines()[1:]     # drop header separator
    rows = {}
    for line in sec:
        c = [x.strip() for x in line.strip("|").split("|")]
        if len(c) < 6:
            continue
        rows[(c[0], c[1], c[2])] = (float(c[3]), float(c[4]), float(c[5].rstrip("%")))
    return rows


try:
    o, n_ = mixture_summary(f"{REF}/extra/BASELINES.md"), mixture_summary(f"{NEW}/extra/BASELINES.md")
    for key in o:
        if key not in n_:
            continue
        tol = family_tol(key[0].replace("viscosity_dyn_mPas", "visc").replace("density_kg_m3", "dens"))
        add("reports/BASELINES.md (mixture summary)", f"{key[0]}.{key[1]}.{key[2]}.n", o[key][0], n_[key][0], 0.0)
        add("reports/BASELINES.md (mixture summary)", f"{key[0]}.{key[1]}.{key[2]}.score", o[key][1], n_[key][1], tol)
except Exception as exc:  # noqa: BLE001
    print("baselines_eval summary parse failed:", exc)

df = pd.DataFrame(ROWS)
df.to_csv("results/t0_repro/comparison.csv", index=False)
print(df.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
print("\nfail/flagged:")
print(df[~df.verdict.str.startswith("pass")].to_string(index=False, float_format=lambda v: f"{v:.6g}"))


