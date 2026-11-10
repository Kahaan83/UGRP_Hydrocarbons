"""v4 Task 10: generate configs/uncertainty_v4.yaml and results/v4/U/README.md from out-of-fold residual files.

Nothing in the yaml is typed by hand.  Residual convention: residual = predicted - measured.
Margins are group-level split-conformal (src/fuelblend/uq/conformal.py), set on DEV groups only; the margin's
coverage on LOCKED groups is a counts-only check and never feeds back into the margin.
Octane and smoke-point per-row out-of-fold residuals were not saved by Tasks 14/15, so this script re-runs those two
(deterministic) evaluations' CV loops and saves the residuals to results/v4/U/*.csv before using them.
Run:  .venv/Scripts/python.exe scripts/v4/build_uncertainty_v4.py
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "v4"))
from fuelblend.uq.conformal import (conformal_rank, coverage, group_scores,  # noqa: E402
                                    one_sided_margin, two_sided_margin)

R = ROOT / "results" / "v4"
OUT = R / "U"
ALPHA = 0.10
HEADER = ("# GENERATED - do not edit. Produced by scripts/v4/build_uncertainty_v4.py from out-of-fold residual files.\n"
          "# Regenerate: .venv/Scripts/python.exe scripts/v4/build_uncertainty_v4.py\n"
          "# Residual convention: predicted - measured.  Margins are group-level split-conformal quantiles set on DEV\n"
          "# groups; value 'inf' means no valid margin exists at this alpha (too few independent groups) and is\n"
          "# deliberately NOT replaced by a typed number.  See results/v4/U/README.md.\n"
          "# CONSUMERS MUST read applied_margin (never the raw quantile in 'value'; a negative quantile is not extra safety).\n"
          "# Entries whose value is the string 'inf' mean NO VALID MARGIN: check isinstance(value, str) and treat as unusable;\n"
          "# do not substitute a number.  Two-sided entries: applied_margin == value (log space for viscosity).\n")
ENTRIES: dict[str, dict] = {}
LOCKED_ROWS: list[dict] = []


def rel(p) -> str:
    return str(Path(p).resolve().relative_to(ROOT)).replace("\\", "/")


def _num(x):
    return None if x is None else float(x)


def add(key, resid, groups, *, kind, direction, unit, model, residual_file, scope, group_unit, caveats,
        locked=None, bias=None, extra=None, log_space=False):
    """kind: 'one' | 'two'.  resid/groups are DEV (out-of-fold) arrays; locked=(resid, groups) or None."""
    resid, groups = np.asarray(resid, float), np.asarray(groups)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        val = (one_sided_margin(resid, groups, ALPHA, direction) if kind == "one"
               else two_sided_margin(resid, groups, ALPHA))
    scores = group_scores(resid, groups, direction, two_sided=(kind == "two"))
    n = int(scores.size)
    k = conformal_rank(n, ALPHA)
    e = dict(value=("inf" if math.isinf(val) else float(val)), unit=unit, direction=direction if kind == "one" else "two_sided",
             alpha=ALPHA, n_groups=n, group_unit=group_unit, conformal_rank_k=k, n_rows=int(np.isfinite(resid).sum()),
             residual_file=residual_file, model=model, scope=scope, residual_convention="predicted - measured")
    if log_space:
        e["space"] = "natural-log residual ln(pred/meas); multiplicative factor = exp(value)"
        if not math.isinf(val):
            e["factor"] = float(math.exp(val))
    if math.isinf(val):
        e["reason"] = (f"only {n} independent groups; a {int((1 - ALPHA) * 100)} % conformal margin needs "
                       f"ceil((n+1)(1-alpha)) = {k} <= n. No valid margin; no substitute number is given.")
    elif kind == "two":
        e["applied_margin"] = float(val)
    else:
        e["applied_margin"] = max(0.0, float(val))
        if val < 0:
            e["applied_margin_note"] = ("computed margin is negative: no/few groups were ever over-predicted at this level; "
                                        "applied_margin is floored at 0 (a negative margin must not be credited as extra safety).")
    if bias is not None:
        e["mean_signed_residual_bias"] = float(bias)
    if locked is not None and len(locked[0]):
        lr, lg = locked
        if kind == "one":
            c, nl = coverage(lr, lg, val, direction)
        else:
            c, nl = coverage(lr, lg, val, two_sided=True)
        e["locked_check"] = dict(groups_covered=c, n_locked_groups=nl,
                                 note="counts only; locked groups never used to set the margin"
                                 + ("; margin is inf so coverage is trivially complete and uninformative" if math.isinf(val) else ""))
        ls = group_scores(lr, lg, direction, two_sided=(kind == "two"))
        e["locked_check"]["locked_group_scores"] = [round(float(s), 4) for s in np.sort(ls)]
        if kind == "one" and not math.isinf(val):
            c2, _ = coverage(lr, lg, max(0.0, val), direction)
            e["locked_check"]["groups_covered_at_applied_margin"] = c2
        LOCKED_ROWS.append(dict(key=key, margin=e["value"], covered=c, n=nl,
                                applied=e["locked_check"].get("groups_covered_at_applied_margin")))
    else:
        e["locked_check"] = None
    e["caveats"] = list(caveats)
    if extra:
        e.update(extra)
    ENTRIES[key] = e
    return e


# ------------------------------------------------------------------------------------------ flash (Task 9)
def flash():
    lab_f = "results/v4/F/flash_lab_loso.csv"
    real_f = "results/v4/F/flash_real_fuels.csv"
    lab = pd.read_csv(ROOT / lab_f)
    real = pd.read_csv(ROOT / real_f)
    d93 = real[real.method == "D93"]
    dev, lk = d93[d93.split_role != "locked"], d93[d93.split_role == "locked"]
    for model, col in (("ideal", "ideal_err_C"), ("unifac", "unifac_err_C")):
        e = add(f"flash_point_C.lab_mixtures.{model}", lab[col], lab.system_id, kind="one", direction="over", unit="degC",
                model=f"Le Chatelier + Antoine + LFL, {'gamma=1' if model == 'ideal' else 'original UNIFAC (Hansen 1991)'}; zero fitted parameters",
                residual_file=lab_f, scope="lab hydrocarbon mixtures (Leenhouts; 126 records), all 31 systems; no dev/locked split exists for lab data",
                group_unit="system (distinct InChIKey set)", bias=lab[col].mean(),
                caveats=["Zero-parameter physics, so leave-one-system-out residuals are genuinely out-of-sample; nothing was tuned on them.",
                         "Mostly binary alkane/aromatic systems (123 of 126 records binary): not representative of a 100+-component jet fuel.",
                         "No locked slice for lab data; instead the margin is cross-checked against the CRC647 D93 fuel residuals (transfer_check)."])
        cc, nn = coverage(real[real.method == "D93"][col], real[real.method == "D93"].fuel_id, e["value"] if e["value"] != "inf" else np.inf, "over")
        e["transfer_check"] = dict(description="fraction of CRC647 D93 fuels whose over-prediction <= the lab margin (different population; a check, not calibration)",
                                   groups_covered=cc, n_groups=nn)
    for model, col in (("ideal", "ideal_err_C"), ("unifac", "unifac_err_C")):
        add(f"flash_point_C.real_fuels_D93_reconstructed.{model}", dev[col], dev.fuel_id, kind="one", direction="over", unit="degC",
            model=f"same physics as lab ({model}) applied to CRC647 reconstructed compositions",
            residual_file=real_f, scope="CRC 647 jet fuels, ASTM D93 closed cup, DEV fuels only (split_role != locked)",
            group_unit="fuel", bias=dev[col].mean(), locked=(lk[col], lk.fuel_id),
            caveats=["Compositions are RECONSTRUCTED from measured aromatics + D86 (ruling R10), not measured by GC: residuals mix reconstruction error with physics error.",
                     "Every fuel is UNDER-predicted (bias about -9 C; D93 closed cup reads warmer than the Le Chatelier prediction), so the conformal over-prediction score is negative. The negative margin only says over-prediction was never observed; it does not make the model accurate.",
                     "applied_margin is floored at 0. Over-prediction can still occur for fuels whose composition is poorly reconstructed or outside CRC647 (diesel/gasoline have no raw flash data).",
                     "ASTM D56 (Tag closed cup, Jet A-1 spec method) has one fuel only (POSF-10325); see the D56 entry."])
    p = real[(real.method == "D56")]
    add("flash_point_C.real_fuels_D56.ideal", p.ideal_err_C.iloc[0:0], p.fuel_id.iloc[0:0], kind="one", direction="over", unit="degC",
        model="ideal Le Chatelier", residual_file=real_f, scope="ASTM D56 fuels with raw flash value",
        group_unit="fuel", caveats=["Exactly one D56 fuel exists (POSF-10325, n=1 headline set, statistics not allowed): residual deliberately NOT used. No D56 margin is possible."],
        extra=dict(n_groups_available=int(len(p))))
    ENTRIES["flash_point_C.real_fuels_D56.ideal"]["n_groups"] = int(len(p))
    ENTRIES["flash_point_C.real_fuels_D56.ideal"]["conformal_rank_k"] = conformal_rank(len(p), ALPHA)
    ENTRIES["flash_point_C.real_fuels_D56.ideal"]["reason"] = ("only 1 independent group (POSF-10325, headline n=1, statistics_allowed=False); "
                                                                 "no valid margin; no substitute number is given.")


# --------------------------------------------------------------------------- density / HoC / freezing (Task 12)
def bulk():
    f = "results/v4/B/density_v4.csv"
    d = pd.read_csv(ROOT / f)
    dev, lk = d[d.slice == "dev"], d[d.slice == "locked"]
    caveat_common = ["Ideal volume mixing of reconstructed CRC647 compositions: error = reconstruction + physics.",
                     "The model is systematically biased low (about -11 kg/m3 on dev); the conformal margin is on |residual| of a biased model, so it is wide. A bias-corrected model would have a narrower margin but its offset would have to be fitted on dev and re-validated; that was not done here.",
                     "Density is only the ideal-volume model; it uses ThermoML pure densities with tiered fallbacks (about 20 % of mass on class regressions)."]
    add("density_kg_m3.crc647_ideal_volume.all_T", dev.err_kg_m3, dev.fuel_id, kind="two", direction="over", unit="kg/m3",
        model="ideal volume mixing, rho = 1/sum(w_i/rho_i(T))", residual_file=f,
        scope="CRC 647 fuels (45 dev), every temperature (-0.2..69.5 C), per-fuel max |err| over its rows",
        group_unit="fuel", bias=dev.err_kg_m3.mean(), locked=(lk.err_kg_m3, lk.fuel_id),
        caveats=caveat_common + ["Resolves the v3 '12.94 vs 2.50' inconsistency: 12.94 kg/m3 was the v3 MAE over 153 fuels (results/v3/G10: MAE 12.94, RMSE 16.75); the 2.50 design margin is unsupported by those residuals (a 90 % quantile below the MAE is impossible). The v4 conformal margin is computed from the 45 real CRC647 dev fuels."])
    d15 = d[np.isclose(d.T_C, 15.56)]
    dev15, lk15 = d15[d15.slice == "dev"], d15[d15.slice == "locked"]
    add("density_kg_m3.crc647_ideal_volume.T15_56C", dev15.err_kg_m3, dev15.fuel_id, kind="two", direction="over", unit="kg/m3",
        model="ideal volume mixing", residual_file=f,
        scope="CRC 647 fuels, ASTM D4052 at 15.56 C only (v3-comparable slice; 45 dev fuels, one row each)",
        group_unit="fuel", bias=dev15.err_kg_m3.mean(), locked=(lk15.err_kg_m3, lk15.fuel_id), caveats=caveat_common)
    h = pd.read_csv(ROOT / "results/v4/B/hoc_v4.csv")
    hd, hl = h[h.split_role != "locked"], h[h.split_role == "locked"]
    for model, col in (("group_additivity_composition_H", "err_comp_H"), ("group_additivity_v2_constH13.8", "err_v2_const_H")):
        add(f"heat_of_combustion_net_MJ_kg.crc647.{model}", hd[col], hd.fuel_id, kind="one", direction="over", unit="MJ/kg",
            model=model, residual_file="results/v4/B/hoc_v4.csv",
            scope="CRC 647 net heat of combustion, 45 dev fuels (min-spec: over-prediction is unsafe)",
            group_unit="fuel", bias=hd[col].mean(), locked=(hl[col], hl.fuel_id),
            caveats=["Reconstructed compositions; group additivity with Tier A/B formation enthalpies (v2 table), no fitted parameters.",
                     "Model under-predicts on average (bias negative), so the over-prediction score can be small or negative; applied_margin is floored at 0."])
    # freezing
    fz = pd.read_csv(ROOT / "results/v4/B/freezing_v4.csv")
    fd, fl = fz[fz.split_role != "locked"].reset_index(drop=True), fz[fz.split_role == "locked"].reset_index(drop=True)
    # leave-one-fuel-out offset: residual_i = err_i - mean(err_{j != i}); locked uses the dev mean (no locked data used)
    e = fd.err_wat_ideal.to_numpy()
    loo = e - (e.sum() - e) / (len(e) - 1)
    add("freezing_point_C.crc647.wat_ideal_minus_loo_offset", loo, fd.fuel_id, kind="one", direction="under", unit="degC",
        model="library ideal WAT minus a constant offset re-estimated leaving each fuel out (locked: dev-mean offset)",
        residual_file="results/v4/B/freezing_v4.csv",
        scope="CRC 647 freezing point (D5972), 45 dev fuels, leave-one-fuel-out offset (max-spec: UNDER-prediction, i.e. predicted colder than measured, is unsafe)",
        group_unit="fuel", bias=float(loo.mean()), locked=(fl.err_wat_ideal - e.mean(), fl.fuel_id),
        caveats=["Task 12 found NO skill vs a constant mean (locked MAE 3.61 vs 3.32 C); this margin is therefore essentially the spread of the fuels' freezing points, not a model-specific accuracy.",
                 "The offset is a fitted constant, here re-estimated out-of-fold; the shipped 'dev-fit' column in freezing_v4.csv is in-sample on dev and was NOT used.",
                 "Reconstructed compositions."])
    m = fd.measured_C.to_numpy()
    cm = (m.sum() - m) / (len(m) - 1) - m
    add("freezing_point_C.crc647.constant_mean_loo", cm, fd.fuel_id, kind="one", direction="under", unit="degC",
        model="no-skill comparator: mean of the other dev fuels' measured freezing points",
        residual_file="results/v4/B/freezing_v4.csv", scope="same as above; comparator, not a recommendation",
        group_unit="fuel", bias=float(cm.mean()), locked=(fd.measured_C.mean() - fl.measured_C, fl.fuel_id),
        caveats=["Included only to show that the physics margin is no tighter than a model that ignores composition entirely."])
    add("freezing_point_C.crc647.coutinho_library", fd.err_coutinho, fd.fuel_id, kind="one", direction="under", unit="degC",
        model="library Coutinho eutectic wax model as shipped", residual_file="results/v4/B/freezing_v4.csv",
        scope="CRC 647 freezing point, 45 dev fuels", group_unit="fuel", bias=fd.err_coutinho.mean(),
        locked=(fl.err_coutinho, fl.fuel_id),
        caveats=["NOT a clean holdout: the library Coutinho parameters were calibrated on Dauphin and CRC data (Task 12), so dev residuals are optimistic.",
                 "Included for completeness; prefer the other entries."])


# ------------------------------------------------------------------------------------------- SLE (Task 13)
def sle():
    f = "results/v4/SLE/sle_eval.csv"
    s = pd.read_csv(ROOT / f)
    s = s[s.scored]
    mg = pd.read_csv(ROOT / "results/v4/SLE/sle_margin.csv").set_index("model")
    ind = {m: dict(point_q90_K=_num(mg.loc[m, "point_q90"]), point_q95_K=_num(mg.loc[m, "point_q95"]),
                   note="Task 13 point-level quantile of (observed - predicted) liquidus on rows where both models exist (points, not independent groups); cluster-bootstrap 90th-percentile 95% CI "
                        f"[{mg.loc[m, 'q90_ci_lo']:.2f}, {mg.loc[m, 'q90_ci_hi']:.2f}] K. Indicative only; NOT a conformal margin.")
           for m in mg.index}
    for model, col, sub in (("ideal_svl", "ideal_svl_K", s[s.ideal_svl_K.notna()]),
                            ("coutinho_lib_common8", "coutinho_lib_K", s[s.ideal_svl_K.notna() & s.coutinho_lib_K.notna()]),
                            ("coutinho_lib_all", "coutinho_lib_K", s[s.coutinho_lib_K.notna()])):
        res = sub[col] - sub.T_obs_K                     # predicted - observed; negative = predicted colder = unsafe for a freezing limit
        extra = {}
        if model == "ideal_svl":
            extra["indicative_not_conformal"] = ind["ideal_svl"]
        if model == "coutinho_lib_common8":
            extra["indicative_not_conformal"] = ind["coutinho_lib"]
        add(f"freezing_point_C.sle_alkane_blends.{model}", res, sub.system_id, kind="one", direction="under", unit="K",
            model={"ideal_svl": "ideal Schroeder-van Laar, mole fractions, verified Tm/dHfus only",
                   "coutinho_lib_common8": "library Coutinho (mass fractions), restricted to the 8 systems where ideal SvL is defined",
                   "coutinho_lib_all": "library Coutinho (mass fractions) on every system it covers"}[model],
            residual_file=f, scope="ThermoML liquidus (SLE) of hydrocarbon binaries/ternaries, P<=150 kPa; no dev/locked split (GroupKFold by system)",
            group_unit="system", bias=res.mean(), extra=extra,
            caveats=["Mostly n-alkane + cyclohexane/toluene blends (about 8 clusters for ideal SvL); not real jet fuel.",
                     "Zero-fitted physics (offset-free rows), so these are out-of-sample; no locked set exists, so no coverage check.",
                     "Unsafe side for a freezing-point limit is predicting a COLDER liquidus than observed."]
            + (["Coutinho library parameters were calibrated partly on Dauphin data, so some of these systems may not be clean holdouts."] if model == "coutinho_lib_all" else []))


# ------------------------------------------------------------------------------------- viscosity (Task 11)
def viscosity():
    f = "results/v4/V/viscosity_predictions_v4.csv"
    v = pd.read_csv(ROOT / f)
    for model, col in (("structured_ridge", "pred_Structured ridge"), ("recalibrated_refutas", "pred_Recalibrated Refutas")):
        v[f"r_{model}"] = np.log(v[col]) - np.log(v["y"])
    dev, lk = v[~v.locked], v[v.locked]
    for model in ("structured_ridge", "recalibrated_refutas"):
        c = f"r_{model}"
        for slice_name, sel in (("all_T", np.ones(len(v), bool)), ("minus40C", np.isclose(v.T_C, -40.0))):
            vd, vl = dev[sel[dev.index]], lk[sel[lk.index]]
            add(f"viscosity_kin_cSt.crc647.{model}.{slice_name}", vd[c], vd.fuel_id, kind="two", direction="over", unit="ln-ratio",
                model=model.replace("_", " ") + " (leave-one-fuel-out, hyperparameters chosen by inner CV)",
                residual_file=f, scope=("CRC 647 kinematic viscosity at +20/-20/-40 C, 45 dev fuels, per-fuel max |ln(pred/meas)| over its rows"
                                         if slice_name == "all_T" else "CRC 647 kinematic viscosity at -40 C only (D455), 45 dev fuels"),
                group_unit="fuel", bias=vd[c].mean(), locked=(vl[c], vl.fuel_id), log_space=True,
                caveats=["Dev rows are leave-one-fuel-out; locked rows predicted by models trained on all 45 dev fuels (Task 11 protocol).",
                         "Reconstructed compositions; only 57 fuels from one source (CRC 647, T range -40..20 C). Not valid for diesel/gasoline or outside the training temperature range.",
                         "Task 11 does not claim a ranking among bilinear / structured ridge / GP; ridge is shown because it is the simplest of the three."]
                + (["Hyperparameters at the edge of their grids (ridge alpha 0.01)."] if model == "structured_ridge" else []))


# ------------------------------------------------------------- octane (Task 15, residuals regenerated) / smoke (Task 14)
def octane():
    import eval_octane_v4 as EO
    from fuelblend.models import octane_blend_v4 as O
    corpus = pd.read_csv(EO.CORPUS)
    sp = json.loads(EO.SPLITS.read_text())
    locked = set(sp["locked"])
    rows = []
    for prop in ("ron", "mon"):
        table, _ = O.prepare_blends(corpus, prop)
        table["is_locked"] = table["blend_id"].isin(locked)
        dev = table[~table.is_locked].reset_index(drop=True)
        lock = table[table.is_locked].reset_index(drop=True)
        folds = list(LeaveOneGroupOut().split(dev, groups=dev.group))
        P, idx = EO.run_view(dev, folds, prop)
        for m in P.columns:
            for i, p in zip(idx, P[m].to_numpy()):
                rows.append(dict(prop=prop, view="dev_leave_one_blendstock_out", model=m, blend_id=dev.blend_id[i], group=dev.group[i],
                                 y=dev.y[i], pred=p, resid=p - dev.y[i]))
        P, idx = EO.run_view(pd.concat([dev, lock], ignore_index=True),
                             [(np.arange(len(dev)), np.arange(len(dev), len(dev) + len(lock)))], prop)
        for m in P.columns:
            for j, p in zip(range(len(lock)), P[m].to_numpy()):
                rows.append(dict(prop=prop, view="locked_fit_on_dev", model=m, blend_id=lock.blend_id[j], group=lock.group[j],
                                 y=lock.y[j], pred=p, resid=p - lock.y[j]))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "octane_oof_residuals.csv", index=False, float_format="%.6g")
    fpath = "results/v4/U/octane_oof_residuals.csv"
    for prop in ("ron", "mon"):
        for m in ("constant_mean", "base_unchanged", "linear_volume_pooled_bov", "ghosh_lumps", "gbm_resid", "gp_resid"):
            a = df[(df.prop == prop) & (df.model == m) & (df.view == "dev_leave_one_blendstock_out")].dropna(subset=["resid"])
            b = df[(df.prop == prop) & (df.model == m) & (df.view == "locked_fit_on_dev")].dropna(subset=["resid"])
            add(f"{prop}.{m}", a.resid, a.group, kind="one", direction="over", unit="octane number", model=m, residual_file=fpath,
                scope=f"Co-Optima bRON blends, {prop.upper()}, dev leave-one-blendstock-out (every held-out blendstock is unseen); min-spec: over-prediction unsafe",
                group_unit="blendstock CAS (BRON:<CAS>)", bias=a.resid.mean(), locked=(b.resid, b.group),
                caveats=["Targets are bRON blends of a base gasoline + ONE blendstock (R17/R19), not market gasolines; only 2 distinct base gasolines exist (RON 90.3/87.7).",
                         "No CI-supported winner among pooled-BOV / Ghosh / GBM / GP (Task 15); all are shown, pick none on margin alone.",
                         "Locked set = 4 blendstock groups (12 blends): coverage counts are nearly uninformative.",
                         "This replaces the hand-typed v3 RON 4.50 / MON 4.00 margins (CLM-28)."]
                + (["base_unchanged and constant_mean are comparators, not candidate models."] if m in ("constant_mean", "base_unchanged") else []))


def smoke():
    import eval_smoke_point as SPE
    from fuelblend.baselines import sooting as S
    a = pd.read_csv(S.DATA / "smoke_point_ajftd.csv")
    a = a[~a.censored].copy().reset_index(drop=True)
    a["tgroup"] = a.fuel_label.astype(str).str.extract(r"\)(T\d+)")[0].fillna("NOCODE")
    y = a.value.to_numpy()
    preds = {"const_mean": np.zeros(len(a)), "selected_in_fold": np.zeros(len(a))}
    for tr, te in GroupKFold(5).split(a, groups=a["tgroup"]):
        tr_, te_ = a.iloc[tr], a.iloc[te]
        preds["const_mean"][te] = tr_.value.mean()
        sel = SPE.select_inside(tr_)
        preds["selected_in_fold"][te] = SPE._pred(te_, SPE.FSETS[sel], SPE._fit(tr_, SPE.FSETS[sel]))
    res = pd.DataFrame(dict(fuel_id=a.fuel_id, fuel_label=a.fuel_label, tgroup=a.tgroup, y=y,
                            pred_const_mean=preds["const_mean"], pred_selected_in_fold=preds["selected_in_fold"]))
    res["resid_const_mean"] = res.pred_const_mean - res.y
    res["resid_selected_in_fold"] = res.pred_selected_in_fold - res.y
    res.to_csv(OUT / "smoke_point_oof_residuals.csv", index=False, float_format="%.6g")
    for m in ("selected_in_fold", "const_mean"):
        add(f"smoke_point_mm.ajftd100_bulk.{m}", res[f"resid_{m}"], res.tgroup, kind="one", direction="over", unit="mm",
            model=("bulk model on aromatics+diaromatics (feature set chosen inside each training fold)" if m == "selected_in_fold" else "no-skill comparator: training-fold mean"),
            residual_file="results/v4/U/smoke_point_oof_residuals.csv",
            scope="AJFTD 100-fuel smoke point (D1322), GroupKFold(5) by inferred T-code (undocumented tank/producer code; unmatched fuels pooled as one group NOCODE); min-spec: over-prediction unsafe",
            group_unit="inferred T-code batch", bias=res[f"resid_{m}"].mean(), locked=None,
            caveats=["Group = inferred T-code, not a documented batch; fuels without a code share one group, so n_groups is an upper bound on independent units.",
                     "Observed smoke points span 19.3-25.0 mm only; the 18/25 mm limits are not probed.",
                     "No locked split exists for AJFTD smoke point: no coverage check. Task 14 found modest skill over the constant mean only when grouped by fuel_id; the T-code grouping is the conservative one.",
                     "The YSI route for pure compounds / real-fuel speciation was not evaluated on real fuels and is excluded."]
            + (["Comparator only."] if m == "const_mean" else []))


# ----------------------------------------------------------------------------------------------- output
FLASH_REC = ("RECOMMENDATION (not a measured result): the D93 reconstructed-CRC entry floors at 0 because every CRC fuel is under-predicted by a "
             "closed-cup/reconstruction offset (about -9.4 C). The lab ideal entry (31 systems) is the only margin with a real over-prediction exposure "
             "(one lab system over-predicts by +10.6 C). Optimizer-designed novel blends must NOT use 0 as the flash margin; use the larger of the lab "
             "applied_margin and the D93 applied_margin unless a different choice is documented.")


def FILT(*w):
    return [list(x) for x in w]


RECOMPUTE = {
    "flash_point_C.lab_mixtures.ideal": dict(column="ideal_err_C", group="system_id", where=[]),
    "flash_point_C.lab_mixtures.unifac": dict(column="unifac_err_C", group="system_id", where=[]),
    "flash_point_C.real_fuels_D93_reconstructed.ideal": dict(column="ideal_err_C", group="fuel_id", where=FILT(("method", "eq", "D93"), ("split_role", "ne", "locked"))),
    "flash_point_C.real_fuels_D93_reconstructed.unifac": dict(column="unifac_err_C", group="fuel_id", where=FILT(("method", "eq", "D93"), ("split_role", "ne", "locked"))),
    "density_kg_m3.crc647_ideal_volume.all_T": dict(column="err_kg_m3", group="fuel_id", where=FILT(("slice", "eq", "dev"))),
    "density_kg_m3.crc647_ideal_volume.T15_56C": dict(column="err_kg_m3", group="fuel_id", where=FILT(("slice", "eq", "dev"), ("T_C", "close", 15.56))),
    "heat_of_combustion_net_MJ_kg.crc647.group_additivity_composition_H": dict(column="err_comp_H", group="fuel_id", where=FILT(("split_role", "ne", "locked"))),
    "heat_of_combustion_net_MJ_kg.crc647.group_additivity_v2_constH13.8": dict(column="err_v2_const_H", group="fuel_id", where=FILT(("split_role", "ne", "locked"))),
    "freezing_point_C.crc647.coutinho_library": dict(column="err_coutinho", group="fuel_id", where=FILT(("split_role", "ne", "locked"))),
}
for _p in ("ron", "mon"):
    for _m in ("constant_mean", "base_unchanged", "linear_volume_pooled_bov", "ghosh_lumps", "gbm_resid", "gp_resid"):
        RECOMPUTE[f"{_p}.{_m}"] = dict(column="resid", group="group", where=FILT(("prop", "eq", _p), ("model", "eq", _m), ("view", "eq", "dev_leave_one_blendstock_out")))
for _m in ("selected_in_fold", "const_mean"):
    RECOMPUTE[f"smoke_point_mm.ajftd100_bulk.{_m}"] = dict(column=f"resid_{_m}", group="tgroup", where=[])
DERIVATION = {
    "freezing_point_C.crc647.wat_ideal_minus_loo_offset": "residual_i = err_wat_ideal_i - mean(err_wat_ideal of the other dev fuels); source column err_wat_ideal in residual_file (dev rows, split_role != locked)",
    "freezing_point_C.crc647.constant_mean_loo": "residual_i = mean(measured_C of the other dev fuels) - measured_C_i; source column measured_C in residual_file (dev rows)",
    "freezing_point_C.sle_alkane_blends.ideal_svl": "residual = ideal_svl_K - T_obs_K on scored rows where ideal_svl_K exists; group = system_id",
    "freezing_point_C.sle_alkane_blends.coutinho_lib_common8": "residual = coutinho_lib_K - T_obs_K on scored rows where both ideal_svl_K and coutinho_lib_K exist; group = system_id",
    "freezing_point_C.sle_alkane_blends.coutinho_lib_all": "residual = coutinho_lib_K - T_obs_K on scored rows where coutinho_lib_K exists; group = system_id",
}
for _m, _c in (("structured_ridge", "pred_Structured ridge"), ("recalibrated_refutas", "pred_Recalibrated Refutas")):
    for _s in ("all_T", "minus40C"):
        DERIVATION[f"viscosity_kin_cSt.crc647.{_m}.{_s}"] = (f"residual = ln({_c}) - ln(y) on rows with locked == False"
                                                              + (" and T_C == -40" if _s == "minus40C" else "") + "; group = fuel_id")


def apply_specs():
    for k, e in ENTRIES.items():
        if k in RECOMPUTE:
            e["recompute"] = RECOMPUTE[k]
        if k in DERIVATION:
            e["residual_derivation"] = DERIVATION[k]


def write():
    ENTRIES["flash_point_C.lab_mixtures.ideal"]["recommendation"] = FLASH_REC
    ENTRIES["flash_point_C.real_fuels_D93_reconstructed.ideal"]["recommendation"] = FLASH_REC
    doc = dict(version="4.0", generated_by="scripts/v4/build_uncertainty_v4.py", alpha=ALPHA,
               residual_convention="predicted - measured; one-sided margins use the unsafe side named in 'direction'",
               margin_rule="conformal quantile k=ceil((n+1)(1-alpha)) of per-group max unsafe-side residuals over n groups; 'inf' if k>n",
               margins=ENTRIES,
               excluded=dict(dvpe_kpa=("EXCLUDED: DVPE targets are UNITS-UNCONFIRMED (Task 15 / SD-03; psi inferred, 51 BOB+oxygenate rows). "
                                       "No margin is produced until units are confirmed from the source."),
                             diesel_cetane_cloud_pour=("No v4 residuals exist (SD-02: 0 diesel fuels qualify); the v3 margins are not carried forward."),
                             pour_point_C="No v4 task evaluated pour point; the v3 hand-typed 4.50 is not carried forward."))
    text = HEADER + yaml.safe_dump(doc, sort_keys=False, default_flow_style=False, width=140)
    (ROOT / "configs" / "uncertainty_v4.yaml").write_text(text, encoding="utf-8")


def fmt(v, d=3):
    return v if isinstance(v, str) else f"{v:.{d}f}"


def readme():
    L = ["# v4 Task 10: conformal uncertainty margins computed from residuals", "",
         "Generated by `scripts/v4/build_uncertainty_v4.py` (do not edit by hand); every number below is read back from `configs/uncertainty_v4.yaml` "
         "or recomputed from the cited residual files. Library: `src/fuelblend/uq/conformal.py`; tests: `tests/test_v4_conformal.py`.", "",
         "## Method",
         "- Residual = predicted - measured. The independent unit is the **group** (fuel for CRC/real fuels, system for lab and SLE data, blendstock CAS for octane, "
         "inferred T-code for AJFTD smoke point); a group's score is its maximum unsafe-side residual, so duplicate rows never inflate n.",
         f"- Margin = the k-th smallest of n group scores, k = ceil((n+1)(1-alpha)), alpha = {ALPHA} (ruling R3: `inf` iff k > n, i.e. fewer than 9 groups at alpha 0.1; "
         "exactly 9 groups gives the maximum). `inf` is written as the string `inf` with a `reason`; no typed substitute is ever supplied.",
         "- Margins are set from **dev** (out-of-fold) residuals only. The **locked** coverage table below is a counts-only check; it was not used to choose or adjust anything.",
         "- One-sided margins can be negative (the model never over-predicted); `applied_margin` floors them at 0. Two-sided margins (density in kg/m3; viscosity in ln-ratio) take the per-group max |residual|.",
         "- Where Tasks 14/15 had not saved per-row residuals (octane, smoke point) this script re-runs their deterministic CV loops and saves the rows to "
         "`results/v4/U/octane_oof_residuals.csv` and `results/v4/U/smoke_point_oof_residuals.csv`.", "",
         "## Entries", "",
         "| key | margin (unit) | n groups | k | direction | bias | locked groups covered |", "|---|---|---|---|---|---|---|"]
    for k, e in ENTRIES.items():
        lc = e["locked_check"]
        cov = f"{lc['groups_covered']}/{lc['n_locked_groups']}" if lc else "n/a"
        val = fmt(e["value"]) + (f" (x{e['factor']:.3f})" if "factor" in e else "")
        L.append(f"| `{k}` | {val} {e['unit']} | {e['n_groups']} | {e['conformal_rank_k']} | {e['direction']} | "
                 f"{fmt(e['mean_signed_residual_bias'], 2) if 'mean_signed_residual_bias' in e else 'n/a'} | {cov} |")
    L += ["", "## What each family says", ""]
    def g(k): return ENTRIES[k]
    fl, fr = g("flash_point_C.lab_mixtures.ideal"), g("flash_point_C.real_fuels_D93_reconstructed.ideal")
    L += [f"### Flash point (Task 9)",
          f"- Lab mixtures (31 systems, ideal Le Chatelier): margin {fmt(fl['value'])} C over {fl['n_groups']} systems; the one large over-prediction (+10.6 C, system 39a914ff0579) sits in the top {fl['n_groups'] - fl['conformal_rank_k'] + 1} of the sorted group scores (the rank k={fl['conformal_rank_k']} excludes the top {fl['n_groups'] - fl['conformal_rank_k']} systems), so it does not set this margin. "
          f"Transfer check on CRC D93 fuels: {fl['transfer_check']['groups_covered']}/{fl['transfer_check']['n_groups']} covered.",
          f"- Real fuels, D93, CRC647 reconstructed compositions, {fr['n_groups']} dev fuels: computed margin {fmt(fr['value'])} C (negative: every fuel is under-predicted, bias {fmt(fr['mean_signed_residual_bias'], 2)} C), applied_margin {fmt(fr.get('applied_margin', 0), 2)} C; "
          f"locked {fr['locked_check']['groups_covered']}/{fr['locked_check']['n_locked_groups']} covered. This is a statement about reconstructed compositions and a closed-cup test; it does not show the model is safe on measured compositions or other fuel types.",
          "- D56 (Tag closed cup): one fuel, margin `inf` (no valid margin).",
          "- **Recommendation for consumers (Task 17), not a measured result:** the D93 reconstructed-CRC entry floors at 0 only because every CRC fuel is under-predicted by a closed-cup/reconstruction offset (about -9.4 C); that offset is not a safety property of the model. "
          f"The lab ideal entry ({fmt(fl['value'], 2)} C, {fl['n_groups']} systems) is the only margin with a real over-prediction exposure (one lab system over-predicts by +10.6 C). "
          "Optimizer-designed novel blends must NOT use 0 as the flash margin; use the larger of the lab applied_margin and the D93 applied_margin unless a different choice is documented.", ""]
    dn, d15 = g("density_kg_m3.crc647_ideal_volume.all_T"), g("density_kg_m3.crc647_ideal_volume.T15_56C")
    L += ["### Density (Task 12) and the 12.94 vs 2.50 inconsistency",
          f"- v3 reported MAE 12.94 kg/m3 (`results/v3/G10`: MAE 12.94, RMSE 16.75 on 153 fuels) next to a v3 design margin of 2.50 kg/m3 that the G10 residuals do not support: a 90 % quantile of 2.50 is impossible when the MAE is 12.94. "
          f"The v4 margin is derived from real CRC647 residuals: {fmt(d15['value'], 2)} kg/m3 at 15.56 C (45 dev fuels) and {fmt(dn['value'], 2)} kg/m3 over all temperatures.",
          f"- The ideal-volume model is biased low (dev mean residual {fmt(d15['mean_signed_residual_bias'], 2)} kg/m3 at 15.56 C), so the two-sided margin is wide; the bias alone accounts for {100 * abs(d15['mean_signed_residual_bias']) / d15['value']:.0f} % of the 15.56 C margin. Locked coverage: {d15['locked_check']['groups_covered']}/{d15['locked_check']['n_locked_groups']} (15.56 C), {dn['locked_check']['groups_covered']}/{dn['locked_check']['n_locked_groups']} (all T).", ""]
    L += ["### Heat of combustion (Task 12)",
          f"- Composition-H group additivity: margin {fmt(g('heat_of_combustion_net_MJ_kg.crc647.group_additivity_composition_H')['value'], 3)} MJ/kg (v2-const-H variant {fmt(g('heat_of_combustion_net_MJ_kg.crc647.group_additivity_v2_constH13.8')['value'], 3)}); supersedes the hand-typed v3 0.28 MJ/kg.", ""]
    L += ["### Freezing point (Tasks 12 and 13)"]
    for k in [k for k in ENTRIES if k.startswith("freezing_point_C")]:
        e = ENTRIES[k]
        L.append(f"- `{k}`: {fmt(e['value'])} {e['unit']} over {e['n_groups']} groups" + ("" if not isinstance(e['value'], str) else f" -- {e['reason']}")
                 + (f" Task 13 indicative (point quantiles, NOT conformal): q90 {fmt(e['indicative_not_conformal']['point_q90_K'], 2)} K, q95 {fmt(e['indicative_not_conformal']['point_q95_K'], 2)} K." if 'indicative_not_conformal' in e else ""))
    fo, fc = g("freezing_point_C.crc647.wat_ideal_minus_loo_offset"), g("freezing_point_C.crc647.constant_mean_loo")
    L += [f"- Task 12: the CRC freezing models show no skill over a constant mean. Consistent with that, the conformal margin of the offset physics model ({fmt(fo['value'], 2)} K) is not tighter than the composition-blind constant-mean comparator ({fmt(fc['value'], 2)} K). "
          f"Task 13's '7 K / 4 K' shorthand corresponds to the stored ideal-SvL point quantiles q95 {fmt(ENTRIES['freezing_point_C.sle_alkane_blends.ideal_svl']['indicative_not_conformal']['point_q95_K'], 2)} K / q90 {fmt(ENTRIES['freezing_point_C.sle_alkane_blends.ideal_svl']['indicative_not_conformal']['point_q90_K'], 2)} K (Coutinho on the same rows: q90 {fmt(ENTRIES['freezing_point_C.sle_alkane_blends.coutinho_lib_common8']['indicative_not_conformal']['point_q90_K'], 2)} K / q95 {fmt(ENTRIES['freezing_point_C.sle_alkane_blends.coutinho_lib_common8']['indicative_not_conformal']['point_q95_K'], 2)} K). They rest on about 8 clusters, are point-level (not group) quantiles, and are kept only under `indicative_not_conformal`; the ideal-SvL conformal margin is `inf` (8 systems < 9).", ""]
    L += ["### Viscosity (Task 11), two-sided in ln space"]
    for k in [k for k in ENTRIES if k.startswith("viscosity")]:
        e = ENTRIES[k]
        L.append(f"- `{k}`: ln-margin {fmt(e['value'])} (factor x{e['factor']:.3f}, i.e. about {100 * (e['factor'] - 1):.0f} %) over {e['n_groups']} fuels; locked {e['locked_check']['groups_covered']}/{e['locked_check']['n_locked_groups']}.")
    L += ["- The -40 C slice is reported separately because Task 11's finding is that raw/recalibrated Refutas collapse at -40 C; the structured ridge is linear in 1000/T by construction.", ""]
    L += ["### Octane (Task 15) and smoke point (Task 14)"]
    for k in [k for k in ENTRIES if k.startswith(("ron.", "mon."))]:
        e = ENTRIES[k]
        L.append(f"- `{k}`: {fmt(e['value'])} ON over {e['n_groups']} blendstock groups; locked {e['locked_check']['groups_covered']}/{e['locked_check']['n_locked_groups']}.")
    L += ["- With 18 dev groups the conformal rank is 18, i.e. the margin equals the **largest** group over-prediction: the finite-sample correction leaves no slack.",
          "- DVPE is excluded: units unconfirmed (SD-03).", ""]
    L += [f"- Smoke point groups: the NOCODE group pools {int((pd.read_csv(OUT / 'smoke_point_oof_residuals.csv').tgroup == 'NOCODE').sum())} fuels with no parseable T-code, so n_groups = 6 is an upper bound on independent units; the result is `inf` either way."]
    for k in [k for k in ENTRIES if k.startswith("smoke_point")]:
        e = ENTRIES[k]
        L.append(f"- `{k}`: {fmt(e['value'])} mm over {e['n_groups']} T-code groups" + (f" -- {e['reason']}" if isinstance(e['value'], str) else "") + ".")
    L += ["", "## Coverage on locked groups (counts only; no tuning)", "", "| entry | margin | locked groups covered / total | covered at applied (floored-at-0) margin |", "|---|---|---|---|"]
    for r in LOCKED_ROWS:
        L.append(f"| `{r['key']}` | {fmt(r['margin'])} | {r['covered']}/{r['n']} | {'' if r['applied'] is None else str(r['applied']) + '/' + str(r['n'])} |")
    L += ["", "Entries without a locked set (lab flash, SLE, AJFTD smoke point) have no row here. A 90 % margin is expected to cover about 90 % of locked groups on average; "
          "with 12 locked fuels (or 4 octane groups) the counts cannot confirm or refute that.", "",
          "## Caveats common to all entries", "",
          "- CRC647 compositions are reconstructed (ruling R10), so every CRC-based margin includes reconstruction error and is conditional on that reconstruction.",
          "- Conformal validity assumes the new group is exchangeable with the calibration groups. A new fuel family (diesel, gasoline, bio-blend) or a new test temperature breaks that assumption; the margins say nothing about them.",
          "- v3 margins in `configs/uncertainty_v3.yaml` are untouched; consumers (Task 17) must read `configs/uncertainty_v4.yaml`.", "",
          "## Claims-register rows that would update (not edited here)", "",
          "- CLM-28 (v3 margins calibrated): flash 4.20, RON 4.50, MON 4.00, DVPE 7.20, DCN 5.50, cloud 3.50, pour 4.50, density 2.50 are replaced by the computed entries above; DVPE, DCN, cloud and pour have no v4 margin.",
          "- CLM-25 (density margin 2.50 kg/m3): the density entries above.",
          "- CLM-26 / CLM-12 / CLM-13 (flash): the flash entries (negative-bias D93 note; no D56 margin).",
          "- CLM-16 / CLM-18 (freezing offset, heat of combustion): the freezing and heat-of-combustion entries."]
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    flash(); bulk(); sle(); viscosity(); octane(); smoke()
    apply_specs(); write(); readme()
    for k, e in ENTRIES.items():
        print(f"{k:70s} {e['value'] if isinstance(e['value'], str) else round(e['value'], 3)!s:>8} n={e['n_groups']}")


if __name__ == "__main__":
    main()
