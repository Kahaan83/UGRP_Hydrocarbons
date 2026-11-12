"""v5 C1: cold-viscosity correction fitted on CRC AV-33-22, scored on AJFTD GCxGC fuels (and the reverse).

Readings R-C1-1..6 and the data/corrections are declared in results/v5/C1/README.md (written before this script was first run).
Sign convention: err = ln(pred) - ln(meas).  A correction adds a ln offset fitted as mean ln(meas/pred_M0) on TRAINING fuels, so
err_K = err_K0 + offset.  The bootstrap difference is MAE(K0) - MAE(K): POSITIVE = the correction is better.

Reuses (does not reimplement): the S3 Refutas-from-measured-composition route (predict_comp, Ctx), the AJFTD mapping and targets
from eval_ajftd_gcxgc_v5.py, the S3 corrected groups/locked flags, the S2 AV-33-22 representative map.  Nothing is written to configs/.
Outputs: results/v5/C1/c1_*.csv
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

V5 = ROOT / "data" / "real_fuels" / "v5"
OUT = ROOT / "results" / "v5" / "C1"
S2 = ROOT / "results" / "v5" / "S2"
S3 = ROOT / "results" / "v5" / "S3"

N_BOOT, SEED, MIN_GROUPS = 5000, 0, 8
BIAS_LIMIT = 0.10
MARGIN_COLD, MARGIN_MID = 0.4154, 0.6084           # S1 two-sided |ln(pred/meas)| margins: T < 263 K, 263-283 K
COLD_T = {"-40C": 233.15, "-20C": 253.15}
T_SNAP_TOL = 0.3
WARM_LO, WARM_HI = 288.15 - 0.3, 373.15 + 0.3        # 15-100 C window for K3


# ---------------------------------------------------------------------------------------------- pure functions (tested)
def margin_for(T_K: float) -> float:
    """Declared S1 margin (ln units) for a temperature; NaN where the README declares none (T >= 283 K)."""
    if T_K < 263.0:
        return MARGIN_COLD
    if T_K <= 283.0:
        return MARGIN_MID
    return np.nan


from fuelblend.eval.cold_offset import apply_k1, fit_k1  # noqa: E402,F401  (moved to src, bodies unchanged)


def fit_k2(T, resid):
    """K2: offset = a + b*1000/T, ordinary least squares over all training rows.  None if fewer than 2 distinct temperatures."""
    T, r = np.asarray(T, float), np.asarray(resid, float)
    ok = np.isfinite(r) & np.isfinite(T)
    x = 1000.0 / T[ok]
    if len(np.unique(np.round(T[ok], 2))) < 2:
        return None
    b, a = np.polyfit(x, r[ok], 1)
    return float(a), float(b)


def apply_k2(ab, T) -> np.ndarray:
    T = np.atleast_1d(np.asarray(T, float))
    return np.full(len(T), np.nan) if ab is None else ab[0] + ab[1] * 1000.0 / T


def walther_extrapolate(T_warm, nu_warm, T_target):
    """K3 (reference): ASTM D341 line ln ln(nu+0.7) = A + B ln T through >= 2 warm points, extrapolated to T_target (cSt)."""
    T_warm, nu_warm = np.asarray(T_warm, float), np.asarray(nu_warm, float)
    ok = np.isfinite(T_warm) & np.isfinite(nu_warm) & (nu_warm > 0)
    if ok.sum() < 2 or len(np.unique(np.round(T_warm[ok], 2))) < 2:
        return np.nan
    B, A = np.polyfit(np.log(T_warm[ok]), np.log(np.log(nu_warm[ok] + 0.7)), 1)
    return float(np.exp(np.exp(A + B * np.log(T_target))) - 0.7)


def split_finite(pred) -> tuple[np.ndarray, int]:
    """Mask of finite predictions and the count excluded (non-finite or non-positive viscosity)."""
    p = np.asarray(pred, float)
    ok = np.isfinite(p) & (p > 0)
    return ok, int((~ok).sum())


def gain_bootstrap(gid, e0, eK, n: int = N_BOOT, seed: int = SEED):
    """Paired bootstrap over groups of MAE(K0) - MAE(K); POSITIVE = correction better.  -> (gain, lo, hi) percentile 95 %."""
    gid, a0, aK = np.asarray(gid), np.abs(np.asarray(e0, float)), np.abs(np.asarray(eK, float))
    ug, inv = np.unique(gid, return_inverse=True)
    G = len(ug)
    s0, sK, cnt = (np.bincount(inv, weights=v, minlength=G) for v in (a0, aK, np.ones(len(a0))))
    ix = np.random.default_rng(seed).integers(0, G, size=(n, G))
    d = s0[ix].sum(1) / cnt[ix].sum(1) - sK[ix].sum(1) / cnt[ix].sum(1)
    return float(a0.mean() - aK.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), float((d > 0).mean())


def verdict(n_groups: int, lo: float, bias_K: float) -> tuple[str, bool]:
    """README success rule: n >= 8, bootstrap interval excludes 0 on the positive side, |bias after correction| <= 0.10 ln."""
    if n_groups < MIN_GROUPS:
        return "n too small", False
    if not lo > 0:
        return "not met (interval does not exclude 0 on the positive side)" + ("" if abs(bias_K) <= BIAS_LIMIT else "; |bias| > 0.10"), False
    if abs(bias_K) > BIAS_LIMIT:
        return "not met (|bias| > 0.10 ln)", False
    return "works", True


def pick_one_per_group(df: pd.DataFrame) -> pd.DataFrame:
    """df columns: entity, group, n_visc_rows.  One entity per group: most viscosity rows, tie by entity id."""
    d = df.sort_values(["group", "n_visc_rows", "entity"], ascending=[True, False, True])
    return d.drop_duplicates("group", keep="first")


# ---------------------------------------------------------------------------------------------- data assembly
def load_module(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def snap(T, nominal) -> float:
    d = np.abs(nominal - T)
    return float(nominal[d.argmin()]) if d.min() <= T_SNAP_TOL else round(float(T), 2)


def assemble():
    """Base table: one row per (entity, T): measured kinematic viscosity, M0 prediction (S3 route), group info."""
    s3 = load_module("s3_for_c1", "scripts/v5/s3_composition_vs_physics.py")
    g1 = load_module("g1_for_c1", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    vlb = load_module("vlb_for_c1", "scripts/v4/build_viscosity_leaderboard_v4.py")
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider

    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    reps_m = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    ik_index = {k: i for i, k in enumerate(reg.ik)}
    ctx = s3.Ctx(reg, prov, vlb.refutas_nu, g1.ideal_sle_C, {})
    fz = pd.read_csv(S3 / "corrected_2026-11-12_fuels.csv").set_index("entity")

    comps, meas, meth = {}, {}, {}
    # AJFTD (G1 mapping, S3 eligibility)
    mapped, _ = g1.build_mapping(reps_m)
    tg = g1.targets(exclude_censored=True)
    tg = tg[(tg.property == "viscosity_kin_cSt") & (tg.T_key > 0)]
    nominal = np.array(sorted(set(np.round(pd.read_csv(V5 / "crc_av3322_properties.csv").query("property=='viscosity_kin_cSt'").T_K.dropna().unique(), 2))))
    tg = tg.assign(T_nom=[snap(t, nominal) for t in tg.T_key])
    for fid, g in mapped.groupby("fuel_id"):
        eid = f"AJFTD:{int(fid)}"
        if eid in fz.index and bool(fz.loc[eid, "M_ok"]):
            comps[eid] = (g.reg_index.to_numpy(np.int64), g.fraction.to_numpy(float))
            q = tg[tg.fuel_id == fid].groupby("T_nom").agg(v=("measured", "median"), m=("method", lambda s: "+".join(sorted(set(s)))))
            meas[eid] = q.v.to_dict()
            meth[eid] = q.m.to_dict()
    # AV-33-22 (S3 route, S2 map)
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    amap = pd.read_csv(S2 / "av3322_representative_map.csv")
    props = pd.read_csv(V5 / "crc_av3322_properties.csv")
    for fid, gdf in comp.groupby("fuel_id"):
        m = gdf.merge(amap[["group", "kind", "carbon_number", "rep_inchikey", "mapped"]], on=["group", "carbon_number"], how="left")
        agg = m[m.mapped.astype(bool)].groupby("rep_inchikey").mass_pct.sum()
        eid = f"AV-33-22:{fid}"
        comps[eid] = (np.array([ik_index[k] for k in agg.index], np.int64), agg.to_numpy(float))
        q = props[(props.fuel_id == fid) & (props.property == "viscosity_kin_cSt") & props.value.notna() & props.T_K.notna()].copy()
        q["Tr"] = q.T_K.round(2)
        gq = q.groupby("Tr").agg(v=("value", "median"), m=("method", lambda s: "+".join(sorted(set(s.astype(str))))), n=("value", "size"))
        meas[eid] = gq.v.to_dict()
        meth[eid] = gq.m.to_dict()
        meth[eid + "#n"] = gq.n.to_dict()
    rows = []
    for eid, (idx, w) in comps.items():
        tt = pd.DataFrame(dict(property="viscosity_kin_cSt", T_K=sorted(meas[eid])))
        if tt.empty:
            continue
        p = s3.predict_comp(ctx, idx, w, tt)
        for r in p.itertuples():
            rows.append(dict(entity=eid, source=eid.split(":")[0], T_K=round(float(r.T_K), 2), measured=meas[eid][float(r.T_K)],
                             method=meth[eid][float(r.T_K)], n_values=meth.get(eid + "#n", {}).get(float(r.T_K), 1), pred_K0=r.pred, status=r.status))
    base = pd.DataFrame(rows)
    base["finite"], _ = split_finite(base.pred_K0)
    base.loc[base.status != "ok", "finite"] = False
    base["err_K0"] = np.where(base.finite, np.log(base.pred_K0.where(base.finite)) - np.log(base.measured), np.nan)
    fz2 = fz.reset_index()[["entity", "group", "group_locked", "locked_ajftd", "label"]]
    base = base.merge(fz2, on="entity", how="left")
    base["n_visc_rows"] = base.groupby("entity").T_K.transform("size")
    return base, dict(ctx=ctx, prov=prov, reg=reg, ik_index=ik_index, comp=comp, amap=amap)


# ---------------------------------------------------------------------------------------------- selection and scoring
def select(base: pd.DataFrame, source: str, T: float, locked: bool = False) -> pd.DataFrame:
    """Scored rows at one temperature: finite K0, one fuel per group (non-locked groups), or the locked AJFTD fuels (kept apart)."""
    d = base[(base.source == source) & (base.T_K == round(T, 2)) & base.finite]
    if locked:
        return d[d.locked_ajftd.fillna(False).astype(bool)]
    d = d[~d.group_locked.fillna(False).astype(bool)]
    return pick_one_per_group(d)


def training_rows(base: pd.DataFrame, source: str, exclude_groups=(), one_per_group: bool = False) -> pd.DataFrame:
    d = base[(base.source == source) & base.finite & ~base.group.isin(list(exclude_groups))]
    if source == "AJFTD":
        d = d[~d.group_locked.fillna(False).astype(bool)]
    if one_per_group:
        d = pd.concat([pick_one_per_group(g) for _, g in d.groupby("T_K")])
    return d


def fit_corrections(train: pd.DataFrame) -> dict:
    resid = -train.err_K0.to_numpy()                       # ln(meas/pred)
    return dict(k1=fit_k1(train.T_K, resid), k2=fit_k2(train.T_K, resid), n_rows=len(train), n_fuels=train.entity.nunique(),
                n_groups=train.group.nunique())


def correct(scored: pd.DataFrame, fit: dict) -> pd.DataFrame:
    s = scored.copy()
    s["off_K1"], s["off_K2"] = apply_k1(fit["k1"], s.T_K.to_numpy()), apply_k2(fit["k2"], s.T_K.to_numpy())
    s["err_K1"], s["err_K2"] = s.err_K0 + s.off_K1, s.err_K0 + s.off_K2
    return s


def k3_errors(base: pd.DataFrame, s: pd.DataFrame) -> pd.Series:
    """K3 error ln(pred/meas) at the row's T from the fuel's own 15-100 C measurements (reference only)."""
    out = []
    for r in s.itertuples():
        w = base[(base.entity == r.entity) & (base.T_K >= WARM_LO) & (base.T_K <= WARM_HI)]
        p = walther_extrapolate(w.T_K, w.measured, r.T_K)
        out.append(np.log(p) - np.log(r.measured) if np.isfinite(p) and p > 0 else np.nan)
    return pd.Series(out, index=s.index)


def cell_metrics(s: pd.DataFrame, k: str, scenario: str, label: str, n_excl_nonfinite: int, extra=None) -> dict:
    """Metrics for one scenario x temperature x correction on the rows where BOTH K0 and K are defined."""
    e = s[np.isfinite(s.err_K0) & np.isfinite(s[f"err_{k}"])]
    n_no_off = int(len(s) - len(e))
    rec = dict(scenario=scenario, T_label=label, correction=k, n_groups=int(e.group.nunique()), n_rows=len(e), n_excluded_nonfinite=n_excl_nonfinite,
               n_excluded_no_correction=n_no_off)
    if extra:
        rec.update(extra)
    if e.empty:
        rec.update(verdict="not evaluable (no rows)", works=False)
        return rec
    e0, eK = e.err_K0.to_numpy(), e[f"err_{k}"].to_numpy()
    gain, lo, hi, pgt0 = gain_bootstrap(e.group, e0, eK)
    T = float(e.T_K.iloc[0])
    mg = margin_for(T)
    v, ok = verdict(rec["n_groups"], lo, eK.mean())
    rec.update(T_K=T, bias_K0=e0.mean(), bias_K=eK.mean(), mae_K0=np.abs(e0).mean(), mae_K=np.abs(eK).mean(), gain_mae_K0_minus_K=gain,
               ci_lo=lo, ci_hi=hi, p_gain_positive=pgt0, margin=mg,
               frac_in_margin_K0=float((np.abs(e0) <= mg).mean()) if np.isfinite(mg) else np.nan,
               frac_in_margin_K=float((np.abs(eK) <= mg).mean()) if np.isfinite(mg) else np.nan,
               verdict=v, works=ok, exploratory=True)
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base, env = assemble()
    base.to_csv(OUT / "c1_base_predictions.csv", index=False, float_format="%.6g")
    nonfinite = base[~base.finite].groupby(["source", "T_K"]).size().to_dict()

    def n_nf(source, T):
        return int(nonfinite.get((source, round(T, 2)), 0))

    metrics, preds, fits = [], [], []

    # ---- R-C1-1 transfer: AV-33-22 (all 38 fuels) -> AJFTD headline groups and the locked fuels apart
    fit_av = fit_corrections(training_rows(base, "AV-33-22"))
    fits.append(dict(scenario="fit on AV-33-22 (all 38 fuels)", n_fuels=fit_av["n_fuels"], n_groups=fit_av["n_groups"], n_rows=fit_av["n_rows"],
                     k2_a=fit_av["k2"][0], k2_b=fit_av["k2"][1], **{f"k1_off_{t:.2f}K": v for t, v in fit_av["k1"].items()}))
    for lab, T in COLD_T.items():
        for locked, scen in ((False, "R-C1-1 transfer AV->AJFTD headline"), (True, "R-C1-1 transfer AV->AJFTD S4-locked (apart)")):
            s = correct(select(base, "AJFTD", T, locked), fit_av)
            s["K3"] = k3_errors(base, s)
            s["scenario"] = scen
            preds.append(s)
            for k in ("K1", "K2"):
                metrics.append(cell_metrics(s, k, scen, lab, n_nf("AJFTD", T)))

    # ---- R-C1-2 within source: leave-one-group-out inside AV-33-22 (the whole group of the scored fuel is held out)
    within = []
    av_groups = sorted(base[base.source == "AV-33-22"].group.unique())
    for g in av_groups:
        fit_g = fit_corrections(training_rows(base, "AV-33-22", exclude_groups=[g]))
        for T in sorted(base[(base.source == "AV-33-22")].T_K.unique()):
            s = select(base, "AV-33-22", T)
            within.append(correct(s[s.group == g], fit_g))
    within = pd.concat(within)
    within["scenario"] = "R-C1-2 within AV-33-22 (leave-one-group-out)"
    for lab, T in COLD_T.items():
        s = within[within.T_K == round(T, 2)].copy()
        s["K3"] = k3_errors(base, s)
        preds.append(s)
        for k in ("K1", "K2"):
            metrics.append(cell_metrics(s, k, "R-C1-2 within AV-33-22 (leave-one-group-out)", lab, n_nf("AV-33-22", T)))

    # ---- R-C1-3 reverse: fit on AJFTD headline groups (non-locked, one fuel per group at each temperature) -> AV-33-22 headline
    tr_aj = training_rows(base, "AJFTD", one_per_group=True)
    fit_aj = fit_corrections(tr_aj)
    fits.append(dict(scenario="fit on AJFTD headline groups (one per group per T)", n_fuels=fit_aj["n_fuels"], n_groups=fit_aj["n_groups"], n_rows=fit_aj["n_rows"],
                     k2_a=fit_aj["k2"][0] if fit_aj["k2"] else np.nan, k2_b=fit_aj["k2"][1] if fit_aj["k2"] else np.nan,
                     **{f"k1_off_{t:.2f}K": v for t, v in fit_aj["k1"].items()}))
    for lab, T in COLD_T.items():
        s = correct(select(base, "AV-33-22", T), fit_aj)
        s["K3"] = k3_errors(base, s)
        s["scenario"] = "R-C1-3 reverse AJFTD->AV-33-22"
        preds.append(s)
        for k in ("K1", "K2"):
            metrics.append(cell_metrics(s, k, "R-C1-3 reverse AJFTD->AV-33-22", lab, n_nf("AV-33-22", T),
                                        extra=dict(train_groups_at_T=int((tr_aj.T_K == round(T, 2)).sum()))))

    # ---- K3 reference rows (never a verdict)
    for sc in ("R-C1-1 transfer AV->AJFTD headline", "R-C1-2 within AV-33-22 (leave-one-group-out)"):
        for lab, T in COLD_T.items():
            s = pd.concat(preds)
            s = s[(s.scenario == sc) & (s.T_K == round(T, 2))]
            s = s.assign(err_K3=s.K3)
            m = cell_metrics(s, "K3", sc, lab, 0)
            m.update(verdict="reference only (uses the scored fuel's own 15-100 C viscosities; never a verdict)", works=False)
            metrics.append(m)

    # ---- source-robust summary per correction x temperature
    M = pd.DataFrame(metrics)
    pd.concat(preds, ignore_index=True).to_csv(OUT / "c1_predictions.csv", index=False, float_format="%.6g")
    M.to_csv(OUT / "c1_metrics.csv", index=False, float_format="%.5g")
    M[["scenario", "T_label", "correction", "n_groups", "gain_mae_K0_minus_K", "ci_lo", "ci_hi", "p_gain_positive"]].to_csv(
        OUT / "c1_bootstrap.csv", index=False, float_format="%.5g")
    pd.DataFrame(fits).to_csv(OUT / "c1_fits.csv", index=False, float_format="%.6g")

    rd = []
    for k in ("K1", "K2"):
        for lab in COLD_T:
            def w(scen):
                r = M[(M.scenario == scen) & (M.T_label == lab) & (M.correction == k)]
                return bool(r.works.iloc[0]) if len(r) else False
            tr, wi, rv = (w("R-C1-1 transfer AV->AJFTD headline"), w("R-C1-2 within AV-33-22 (leave-one-group-out)"), w("R-C1-3 reverse AJFTD->AV-33-22"))
            rd.append(dict(reading="R-C1-1", correction=k, T_label=lab, works=tr))
            rd.append(dict(reading="R-C1-2", correction=k, T_label=lab, works=wi, note="within-source success alone does not count"))
            rd.append(dict(reading="R-C1-3", correction=k, T_label=lab, works=rv, source_robust=bool(tr and rv)))
    pd.DataFrame(rd).to_csv(OUT / "c1_readings.csv", index=False)

    # ---- R-C1-4/5 by temperature: K0 bias and margin fractions for each source; leave-one-group-out K1/K2 on AV for all 9 T
    bt = []
    for src in ("AV-33-22", "AJFTD"):
        for T in sorted(base[base.source == src].T_K.unique()):
            s = select(base, src, T)
            if s.empty:
                continue
            e0 = s.err_K0.to_numpy()
            mg = margin_for(T)
            rec = dict(source=src, T_K=T, T_C=round(T - 273.15, 2), n_groups=s.group.nunique(), n_excluded_nonfinite=n_nf(src, T), bias_K0=e0.mean(),
                       mae_K0=np.abs(e0).mean(), margin=mg, frac_in_margin_K0=float((np.abs(e0) <= mg).mean()) if np.isfinite(mg) else np.nan)
            if src == "AV-33-22":
                w_ = within[within.T_K == round(T, 2)]
                for k in ("K1", "K2"):
                    rec[f"loo_bias_{k}"] = w_[f"err_{k}"].mean()
                    rec[f"loo_mae_{k}"] = w_[f"err_{k}"].abs().mean()
            bt.append(rec)
    pd.DataFrame(bt).to_csv(OUT / "c1_by_temperature.csv", index=False, float_format="%.5g")

    # ---- av3322 methods per T (which values were medianed)
    av = base[base.source == "AV-33-22"].groupby("T_K").agg(n_fuels=("entity", "nunique"), fuels_with_multiple_values=("n_values", lambda v: int((v > 1).sum())),
                                                          methods=("method", lambda m: "; ".join(sorted(set("+".join(m).split("+"))))))
    av.to_csv(OUT / "c1_av3322_temperature_methods.csv", float_format="%.5g")

    pure_check(base, env)
    print(M.drop(columns=["exploratory"], errors="ignore").to_string())
    print(pd.DataFrame(bt).to_string())


def pure_check(base, env) -> None:
    """R-C1-6: the 5 most abundant AV-33-22 representatives; provider pure nu at -20/-40 C versus Boehm 2022 pure data and the provider's data range."""
    from fuelblend.mixture import common
    comp, amap, prov, reg, ik_index = env["comp"], env["amap"], env["prov"], env["reg"], env["ik_index"]
    m = comp.merge(amap[["group", "carbon_number", "rep_inchikey", "mapped"]], on=["group", "carbon_number"], how="left")
    m = m[m.mapped.astype(bool)]
    share = m.groupby(["fuel_id", "rep_inchikey"]).mass_pct.sum().groupby("rep_inchikey").sum() / comp.fuel_id.nunique()
    top = share.sort_values(ascending=False).head(5)
    boehm = pd.read_csv(V5 / "boehm2022_viscosity.csv")
    boehm = boehm[(boehm.kind == "pure") & (boehm.property == "viscosity_kin_cSt")].copy()
    boehm["ik"] = boehm.comp_inchikeys.map(lambda s: json.loads(s)[0])
    fe = common._get_measured_fits()[1]
    rows = []
    for ik, mass in top.items():
        i = ik_index[ik]
        for lab, T in COLD_T.items():
            le, rho = prov.lookup(np.array([[i]]), T)
            nu = float(1000.0 * np.exp(le[0, 0]) / rho[0, 0])
            b = boehm[(boehm.ik == ik) & (np.abs(boehm.T_K - T) < 0.3)]
            rec = dict(inchikey=ik, smiles=reg.smiles[i], mean_mass_pct=mass, T_label=lab, provider_nu_cSt=nu, provider_ln_nu=np.log(nu))
            if len(b):
                bv = float(b.value.median())
                rec.update(boehm_nu_cSt=bv, ln_ratio_provider_over_boehm=float(np.log(nu / bv)), boehm_n=len(b), evaluable=True)
            else:
                rec.update(boehm_nu_cSt=np.nan, ln_ratio_provider_over_boehm=np.nan, boehm_n=0, evaluable=False,
                           note="not evaluable: no Boehm 2022 pure kinematic viscosity for this compound at this T in data/real_fuels/v5/boehm2022_viscosity.csv")
            if ik in fe:
                a, bb, tlo, thi = fe[ik]
                rec.update(provider_fit_T_range=f"{tlo:.1f}-{thi:.1f}", provider_source="measured pure fit (model_table.csv, ln eta linear in ln T)",
                           provider_outside_fit_range=bool(T < tlo or T > thi),
                           provider_extrapolating=bool(T < tlo - 25.0 or T > thi + 25.0))     # beyond the 25 K tolerance the provider switches to its ML predictor
            else:
                rec.update(provider_fit_T_range="none", provider_source="ML predictor (no measured pure viscosity fit)", provider_outside_fit_range=True,
                           provider_extrapolating=True)
            rows.append(rec)
    pd.DataFrame(rows).to_csv(OUT / "c1_pure_check.csv", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
