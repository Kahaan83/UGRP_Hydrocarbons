"""v5 delta_newdata_freeze_smoke: score the shipped freezing and smoke models on the newly ingested literature data.

Declared readings: results/v5/BASE/delta_newdata_freeze_smoke/README.md (committed before this script existed).
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_newdata_freeze_smoke.py
Nothing here changes a model; predict_fuel / smoke_point_from_bulk are called as shipped.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "v5" / "BASE" / "delta_newdata_freeze_smoke"

from fuelblend.api.predict_route import get_ctx, predict_fuel                      # noqa: E402
from fuelblend.data import loaders_freezing_lit as LF                               # noqa: E402
from fuelblend.data import loaders_smoke_lit as LS                                  # noqa: E402
from fuelblend.eval.composition_route import kind_bins, map_bins                    # noqa: E402
from fuelblend.uq.conformal import group_scores, one_sided_margin                                 # noqa: E402

ALPHA = 0.10
SLE_CSV = ROOT / "data" / "processed" / "v5" / "thermoml_sle_hc.csv"
T_REF_K = 288.15
RNG = np.random.default_rng(0)


def _nc(smi):
    return sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(smi).GetAtoms())


def _ik(smi):
    return inchi.MolToInchiKey(Chem.MolFromSmiles(smi))


def loso_const(meas, groups):
    """Leave-one-group-out mean of the measured values (per row: mean over rows of the other groups)."""
    meas, groups = np.asarray(meas, float), np.asarray(groups).astype(str)
    out = np.full(len(meas), np.nan)
    for g in np.unique(groups):
        m = groups == g
        if (~m).any():
            out[m] = meas[~m].mean()
    return out


def margin(resid, groups, direction):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(one_sided_margin(resid, groups, alpha=ALPHA, direction=direction))


def summarize(df, scope, model, pred_col, direction, group_col="group", extra=None):
    """df: rows of the scope, columns measured, censored(bool), pred_col, group_col.  Rows with non-finite prediction are counted unscored."""
    d = df.copy()
    n_scope = len(d)
    cens = d.censored.astype(bool)
    ok = (~cens) & np.isfinite(d[pred_col].astype(float)) & np.isfinite(d.measured.astype(float))
    s = d[ok]
    row = dict(scope=scope, model=model, n_rows_scope=n_scope, n_censored=int(cens.sum()), n_scored=int(len(s)),
               n_unscored=int(n_scope - cens.sum() - len(s)), n_systems=int(s[group_col].nunique()))
    if len(s):
        res = (s[pred_col] - s.measured).to_numpy(float)
        cst = loso_const(s.measured, s[group_col])
        okc = np.isfinite(cst)
        row.update(mae=float(np.abs(res).mean()), bias=float(res.mean()),
                   mae_const=float(np.abs(cst[okc] - s.measured.to_numpy()[okc]).mean()) if okc.any() else np.nan,
                   mae_macro=float(pd.Series(np.abs(res)).groupby(s[group_col].astype(str).to_numpy()).mean().mean()),
                   mae_const_macro=float(pd.Series(np.abs(cst - s.measured.to_numpy())).groupby(s[group_col].astype(str).to_numpy()).mean().mean()) if okc.all() else np.nan,
                   margin_alpha10=margin(res, s[group_col], direction), margin_direction=direction,
                   measured_min=float(s.measured.min()), measured_max=float(s.measured.max()),
                   frac_unsafe=float((res < 0).mean() if direction == "under" else (res > 0).mean()),
                   max_group_unsafe=float(group_scores(res, s[group_col], direction).max()))
    else:
        row.update(mae=np.nan, bias=np.nan, mae_const=np.nan, mae_macro=np.nan, mae_const_macro=np.nan,
                   margin_alpha10=float("inf"), margin_direction=direction, measured_min=np.nan, measured_max=np.nan, frac_unsafe=np.nan)
    if extra:
        row.update(extra(s))
    return row


def bootstrap_diff(s, pred_col, group_col="group", B=2000):
    """System bootstrap of (constant MAE - model MAE), row level; positive = model better."""
    res = (s[pred_col] - s.measured).abs().to_numpy(float)
    cst = np.abs(loso_const(s.measured, s[group_col]) - s.measured.to_numpy(float))
    gs = s[group_col].astype(str).to_numpy()
    ug = np.unique(gs)
    idx = {g: np.flatnonzero(gs == g) for g in ug}
    diffs = []
    for _ in range(B):
        pick = RNG.choice(ug, len(ug))
        ii = np.concatenate([idx[g] for g in pick])
        diffs.append(cst[ii].mean() - res[ii].mean())
    return [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))]


# ======================================================================================================================
# freezing
# ======================================================================================================================
def freeze_call(comp, mole=None, kind="known"):
    """comp: list of (smiles, weight).  Returns dict(svl, f1x, status, sle_missing, warnings)."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = predict_fuel(comp, fuel_type="jet", mole_fractions=mole, composition_kind=kind)
    except Exception as e:                                           # noqa: BLE001 - recorded, never silent
        return dict(svl=np.nan, f1x=np.nan, status=f"error:{type(e).__name__}:{e}", sle_missing=np.nan, warn="")
    f = r["freezing_point_C"]
    v = f.get("value")
    f1 = (f.get("variants") or {}).get("f1x_recalibrated")
    return dict(svl=float(v) if v is not None and np.isfinite(v) else np.nan,
                f1x=float(f1["value"]) if f1 else np.nan, status=f.get("status", ""),
                sle_missing=float(f.get("coverage", {}).get("sle_missing_mapped_mass", np.nan)), warn=";".join(f.get("warnings", [])))


def thermoml_rows():
    d = pd.read_csv(SLE_CSV)
    d = d[d.point_kind == "liquidus"].reset_index(drop=True)
    out = []
    for r in d.itertuples():
        smi, x = json.loads(r.smiles), [float(v) for v in json.loads(r.x_liquid)]
        ncs = [_nc(s) for s in smi]
        c = freeze_call(list(zip(smi, x)), mole=x)
        out.append(dict(dataset="thermoml", source=r.doi, system_id=r.system_id, group=r.system_id, sample=" + ".join(json.loads(r.names)),
                        T_K=r.T_K, x_json=json.dumps(x), n_comp=r.n_comp, all_c8_c16=all(8 <= n <= 16 for n in ncs),
                        measured=float(r.T_K) - 273.15, censored=False, qualifiers="", **{f"pred_{k}": v for k, v in c.items()},
                        comp_scale=r.comp_scale))
    return pd.DataFrame(out)


def ascent_rows():
    p = LF.load_properties()
    p = p[(p.source_id.str.startswith("ascent065A")) & (p.property == "freezing_point_C")].copy()
    rawsmi = {}
    for fn in ("ascent065A_freezing.csv",):
        raw = pd.read_csv(LF.RAW / fn, keep_default_na=False)
        for r in raw.itertuples():
            for s in json.loads(r.smiles):
                rawsmi[_ik(s)] = s
    rows = []
    for r in p.itertuples():
        iks, fr = json.loads(r.comp_inchikeys), json.loads(r.fractions_norm)
        c = freeze_call([(rawsmi[k], f) for k, f in zip(iks, fr)], mole=fr)
        rows.append(dict(dataset="ascent", source=r.source_id, system_id=r.system_id, group=r.system_id, sample=r.sample_id, measured=r.value,
                         censored=False, qualifiers=r.qualifier, kind=r.kind, key=(r.system_id, tuple(round(v, 4) for v in fr)),
                         **{f"pred_{k}": v for k, v in c.items()}))
    d = pd.DataFrame(rows)
    # replicate averaging per (family, system, composition); predictions are identical within a key
    agg = []
    for k, g in d.groupby(d.key.astype(str), sort=False):
        a = g.iloc[0].copy()
        a["measured"] = float(g.measured.mean())
        a["n_obs"] = len(g)
        a["qualifiers"] = ";".join(sorted(set(q for q in g.qualifiers if q)))
        a["sample"] = " | ".join(g["sample"])
        agg.append(a)
    d = pd.DataFrame(agg).drop(columns=["key"]).reset_index(drop=True)
    d["scope_kind"] = np.where(d.kind == "pure", "ascent_pure", "ascent_mixtures")
    return d


def vozka_rows():
    ctx = get_ctx()
    p = LF.load_properties()
    p = p[(p.source_id == "vozka2018_si") & (p.property == "freezing_point_C")]
    rows = []
    for r in p.itertuples():
        bins = LF.blend_composition(r.sample_id)
        bins = kind_bins("ajftd", bins[["hc_class", "n_carbon", "mass_pct"]])
        mf = map_bins(bins, ctx.reps)
        sid = r.sample_id
        grp = "JetA|CAME" if sid in ("Jet A", "CAME") or sid.startswith("C-") else ("JetA|TALL" if sid == "TALL" or sid.startswith("T-") else "JetA|MFAT")
        base = dict(dataset="vozka", source=r.source_id, system_id=r.system_id, group=grp, sample=sid, measured=r.value, censored=False,
                    qualifiers=r.qualifier, unmapped_fraction=mf.unmapped_fraction)
        if mf.low_coverage:
            rows.append({**base, "pred_svl": np.nan, "pred_f1x": np.nan, "pred_status": "low_coverage", "pred_sle_missing": np.nan, "pred_warn": ""})
            continue
        comp = {ctx.reg.ik[i]: float(w) for i, w in zip(mf.idx, mf.w)}
        c = freeze_call(list(comp.items()), kind="reconstructed")
        rows.append({**base, **{f"pred_{k}": v for k, v in c.items()}})
    return pd.DataFrame(rows)


def kittel_freezing():
    p = LF.load_properties()
    k = p[(p.source_id == "kittel2023") & (p.property == "freezing_point_C")]
    return pd.DataFrame(dict(dataset="kittel2023", sample=k.sample_id, measured=k.value, qualifiers=k.qualifier,
                             status="not_scorable_no_composition", censored=k.qualifier.isin(["<", ">"])))


def run_freezing():
    th, asc, voz, kit = thermoml_rows(), ascent_rows(), vozka_rows(), kittel_freezing()
    th["scope_flags"] = np.where(th.all_c8_c16 & (th.T_K <= 273.15), "thermoml_jet_cold;thermoml_all", "thermoml_all")
    asc["scope_flags"] = asc.scope_kind
    voz["scope_flags"] = "vozka"
    rows = pd.concat([th, asc.drop(columns=["scope_kind", "kind", "n_obs"], errors="ignore"), voz], ignore_index=True)
    return rows, kit


def freezing_summary(rows, kit):
    S = []
    scopes = {s: rows[rows.scope_flags.str.split(";").map(lambda l, s=s: s in l)] for s in
              ("thermoml_jet_cold", "thermoml_all", "ascent_mixtures", "ascent_pure", "vozka")}
    for sc, d in scopes.items():
        S.append(summarize(d, sc, "svl_designated", "pred_svl", "under"))
        d2 = d[np.isfinite(d.pred_f1x.astype(float))]
        if len(d2):
            S.append(summarize(d2, sc, "f1x_variant", "pred_f1x", "under"))
    return S, scopes


# ======================================================================================================================
# smoke
# ======================================================================================================================
def smiles_map():
    m = {}
    for fn in LS.FILES:
        raw = pd.read_csv(LS.RAW / fn, keep_default_na=False)
        for r in raw.itertuples():
            for s in json.loads(r.smiles):
                m[_ik(s)] = s
    return m


def smoke_scope(r):
    if r.kind == "real_fuel":
        return "excluded_no_recipe"
    if r.kind == "pure":
        return "smoke_pure"
    if r.source_id == "S4_rainville2025":
        return "smoke_mixtures" if r.kind == "surrogate" else "smoke_rainville_calib"
    return "smoke_mixtures"


def smoke_call(iks, fr, basis, smi_map, ctx):
    """predict_fuel on one recipe -> dict(bulk, ysi, arom, in_domain, ysi_cov, status)."""
    try:
        smi = [smi_map[k] for k in iks]
        if len(iks) == 1:
            fr, basis = [1.0], "mass"                                   # pure compound (the loader stores no fractions for it)
        if basis == "mole":
            res = predict_fuel(list(zip(smi, fr)), mole_fractions=fr)
        elif basis == "mass" or basis == "none":
            res = predict_fuel(list(zip(smi, fr)))
        elif basis == "volume":
            idx = [ctx.reg.ik2i.get(k) for k in iks]
            if any(i is None for i in idx):
                return dict(status="volume_to_mass_needs_density:not_in_registry")
            _, rho = ctx.prov.lookup(np.array(idx, np.int64)[None, :], T_REF_K)
            rho = np.asarray(rho)[0]
            if not np.all(np.isfinite(rho)):
                return dict(status="volume_to_mass_needs_density:no_grid_value")
            w = np.array(fr) * rho
            res = predict_fuel(list(zip(smi, w / w.sum())))
        else:
            return dict(status=f"unknown_basis:{basis}")
    except Exception as e:                                           # noqa: BLE001 - recorded, never silent
        return dict(status=f"error:{type(e).__name__}:{e}")
    sp = res["smoke_point_mm"]
    ysi = (sp.get("variants") or {}).get("ysi_route") or {}
    yv = ysi.get("value")
    return dict(status=sp.get("status", ""), bulk=sp.get("value"), ysi=yv if (ysi.get("status") == "ok" and yv is not None) else np.nan,
                ysi_cov=ysi.get("ysi_own_mass", np.nan), arom=res["aromatics_vol_pct"].get("value"),
                in_domain=not any("out_of_domain" in w for w in sp.get("warnings", [])),
                warn=";".join(sp.get("warnings", [])))


def smoke_rows():
    ctx = get_ctx()
    smi_map = smiles_map()
    d = LS.load_smoke_lit()
    d = d[d.property == "smoke_point_mm"].copy()
    d["scope"] = [smoke_scope(r) for r in d.itertuples()]
    d["excluded_qualifier"] = d.qualifier.isin(["reference", "reported_mean"])
    keep = d[~d.excluded_qualifier].copy()
    keep["key"] = [(r.source_id, r.system_id, r.fractions, r.fraction_basis, r.scope) for r in keep.itertuples()]
    out = []
    for k, g in keep.groupby(keep.key.astype(str), sort=False):
        r = g.iloc[0]
        base = dict(source=r.source_id, system_id=r.system_id, group=r.system_id, sample=" | ".join(sorted(set(g.sample_id))), kind=r.kind, scope=r.scope,
                    basis=r.fraction_basis, n_obs=len(g), measured=float(g.value.mean()),
                    qualifiers=";".join(sorted(set(q for q in g.qualifier if q))), n_comp=int(r.n_comp))
        if r.scope == "excluded_no_recipe":
            out.append({**base, "pred_status": "excluded_no_recipe", "censored": False})
            continue
        c = smoke_call(json.loads(r.comp_inchikeys), json.loads(r.fractions), r.fraction_basis, smi_map, ctx)
        out.append({**base, "censored": False, "pred_status": c.get("status", ""), "pred_bulk": c.get("bulk", np.nan), "pred_ysi": c.get("ysi", np.nan),
                    "ysi_own_cov": c.get("ysi_cov", np.nan), "aromatics_vol_pct": c.get("arom", np.nan), "in_domain": c.get("in_domain", np.nan),
                    "pred_warn": c.get("warn", "")})
    mix = pd.DataFrame(out)
    # Kittel: bulk regression on the measured aromatics vol %
    from fuelblend.baselines.sooting import OutOfDomainWarning, smoke_point_from_bulk
    k = LF.load_properties()
    k = k[k.source_id == "kittel2023"]
    sp, ar = k[k.property == "smoke_point_mm"].set_index("sample_id"), k[k.property == "aromatics_vol_pct"].set_index("sample_id")
    ko = []
    for sid in ar.index:
        a = float(ar.loc[sid, "value"])
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            v = smoke_point_from_bulk(a)
        indom = not any(issubclass(w.category, OutOfDomainWarning) for w in rec)
        if sid in sp.index:
            q = sp.loc[sid, "qualifier"]
            ko.append(dict(source="kittel2023", system_id=f"fuel:kittel2023:{sid.replace(' ', '_')}", group=f"fuel:kittel2023:{sid.replace(' ', '_')}", sample=sid,
                           kind="real_fuel", scope="smoke_kittel", basis="none", n_obs=1, measured=float(sp.loc[sid, "value"]), qualifiers=q,
                           censored=q in ("<", ">"), pred_status="ok", pred_bulk=v, pred_ysi=np.nan, ysi_own_cov=np.nan, aromatics_vol_pct=a, in_domain=indom,
                           pred_warn=""))
        else:
            ko.append(dict(source="kittel2023", system_id=f"fuel:kittel2023:{sid.replace(' ', '_')}", group=f"fuel:kittel2023:{sid.replace(' ', '_')}", sample=sid,
                           kind="real_fuel", scope="smoke_kittel", basis="none", n_obs=0, measured=np.nan, qualifiers="no_smoke_value", censored=False,
                           pred_status="no_measurement", pred_bulk=v, pred_ysi=np.nan, aromatics_vol_pct=a, in_domain=indom, pred_warn=""))
    return pd.concat([mix, pd.DataFrame(ko)], ignore_index=True)


def smoke_summary(rows):
    S, scopes = [], {}
    for sc in ("smoke_mixtures", "smoke_pure", "smoke_rainville_calib", "smoke_kittel"):
        d = rows[rows.scope == sc].copy()
        d = d[d.pred_status != "no_measurement"]
        scopes[sc] = d

        def ex(s):
            if "in_domain" not in s or not len(s):
                return {}
            ind = s.in_domain.astype(bool)
            r = (s.pred_bulk - s.measured).abs()
            cst = np.abs(loso_const(s.measured, s["group"]) - s.measured.to_numpy(float))
            return dict(n_in_domain=int(ind.sum()), n_out_domain=int((~ind).sum()),
                        mae_const_in_domain=float(cst[ind.to_numpy()].mean()) if ind.any() else np.nan,
                        mae_const_out_domain=float(cst[~ind.to_numpy()].mean()) if (~ind).any() else np.nan,
                        mae_in_domain=float(r[ind].mean()) if ind.any() else np.nan, mae_out_domain=float(r[~ind].mean()) if (~ind).any() else np.nan)
        S.append(summarize(d, sc, "bulk_designated", "pred_bulk", "over", extra=ex))
        if sc != "smoke_kittel":
            dy = d[np.isfinite(d.pred_ysi.astype(float))]
            r = summarize(dy, sc, "ysi_variant", "pred_ysi", "over")
            r["n_gate_pass"] = int(len(dy))
            r["gate_coverage"] = float(len(dy) / max(1, int(d.pred_bulk.notna().sum())))
            S.append(r)
    return S, scopes


# ======================================================================================================================
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frz, kit = run_freezing()
    fS, fScopes = freezing_summary(frz, kit)
    smk = smoke_rows()
    sS, sScopes = smoke_summary(smk)
    summary = pd.DataFrame(fS + sS)
    frz.drop(columns=["x_json"], errors="ignore").to_csv(OUT / "rows_freezing.csv", index=False)
    smk.to_csv(OUT / "rows_smoke.csv", index=False)
    kit.to_csv(OUT / "rows_kittel_freezing_not_scorable.csv", index=False)
    summary.to_csv(OUT / "summary.csv", index=False)

    def get(sc, m):
        return summary[(summary.scope == sc) & (summary.model == m)].iloc[0].to_dict()

    r1 = get("thermoml_jet_cold", "svl_designated")
    d1 = fScopes["thermoml_jet_cold"]
    d1 = d1[np.isfinite(d1.pred_svl.astype(float))].copy()
    r1["bootstrap95_const_minus_svl"] = bootstrap_diff(d1, "pred_svl") if len(d1) else None
    pass1 = bool(r1["mae"] < r1["mae_const"])
    rn = {}
    readings = {
        "R-NF-1": dict(scope="thermoml_jet_cold", pass_point_estimate=pass1, svl_mae_C=r1["mae"], const_mae_C=r1["mae_const"], svl_bias_C=r1["bias"],
                       macro_svl=r1["mae_macro"], macro_const=r1["mae_const_macro"], n_rows=r1["n_scored"], n_systems=r1["n_systems"],
                       n_unscored=r1["n_unscored"], bootstrap95_const_minus_svl=r1["bootstrap95_const_minus_svl"],
                       interval_excludes_zero=bool(r1["bootstrap95_const_minus_svl"] and r1["bootstrap95_const_minus_svl"][0] > 0)),
        "R-NF-2": {sc: dict(model="svl_designated", direction="under", n_systems=get(sc, "svl_designated")["n_systems"],
                            margin_C=get(sc, "svl_designated")["margin_alpha10"], finite=bool(np.isfinite(get(sc, "svl_designated")["margin_alpha10"])),
                            n_rows=get(sc, "svl_designated")["n_scored"]) for sc in ("thermoml_jet_cold", "thermoml_all", "ascent_mixtures", "vozka")},
    }
    readings["R-NF-2"]["pass"] = readings["R-NF-2"]["thermoml_jet_cold"]["finite"]
    rb, ry = get("smoke_mixtures", "bulk_designated"), get("smoke_mixtures", "ysi_variant")
    sb = sScopes["smoke_mixtures"]
    sb_b = sb[np.isfinite(sb.pred_bulk.astype(float)) & ~sb.censored.astype(bool)]
    sb_y = sb[np.isfinite(sb.pred_ysi.astype(float))]
    boot_b, boot_y = bootstrap_diff(sb_b, "pred_bulk"), bootstrap_diff(sb_y, "pred_ysi")
    readings["R-NF-3"] = dict(post_hoc_bootstrap95_const_minus_model=dict(bulk=boot_b, ysi=boot_y), scope="smoke_mixtures",
                              bulk=dict(pass_point_estimate=bool(rb["mae"] < rb["mae_const"]), mae_mm=rb["mae"], const_mae_mm=rb["mae_const"], bias_mm=rb["bias"],
                                        n_rows=rb["n_scored"], n_systems=rb["n_systems"], n_unscored=rb["n_unscored"],
                                        n_in_domain=rb.get("n_in_domain"), n_out_domain=rb.get("n_out_domain"),
                                        mae_in_domain=rb.get("mae_in_domain"), mae_out_domain=rb.get("mae_out_domain")),
                              ysi=dict(pass_point_estimate=bool(ry["mae"] < ry["mae_const"]) if ry["n_scored"] else None, mae_mm=ry["mae"], const_mae_mm=ry["mae_const"],
                                       bias_mm=ry["bias"], n_rows=ry["n_scored"], n_systems=ry["n_systems"], gate_coverage=ry.get("gate_coverage")))
    readings["R-NF-4"] = dict(scope="smoke_mixtures", direction="over",
                              bulk=dict(n_systems=rb["n_systems"], margin_mm=rb["margin_alpha10"], finite=bool(np.isfinite(rb["margin_alpha10"])),
                                        measured_span_mm=[rb["measured_min"], rb["measured_max"]]),
                              ysi=dict(n_systems=ry["n_systems"], margin_mm=ry["margin_alpha10"], finite=bool(np.isfinite(ry["margin_alpha10"])),
                                       measured_span_mm=[ry["measured_min"], ry["measured_max"]]))
    readings["counts"] = dict(freezing_rows=len(frz), thermoml_liquidus_rows=int((frz.dataset == "thermoml").sum()), ascent_points=int((frz.dataset == "ascent").sum()),
                              vozka_rows=int((frz.dataset == "vozka").sum()), kittel_freezing_not_scorable=len(kit),
                              kittel_freezing_censored=int(kit.censored.sum()), smoke_points=len(smk),
                              smoke_mixtures_points=int((smk.scope == "smoke_mixtures").sum()))
    (OUT / "readings.json").write_text(json.dumps(readings, indent=2, default=float), encoding="utf-8")
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(summary.to_string())
    print(json.dumps(readings, indent=2, default=float))


if __name__ == "__main__":
    main()
