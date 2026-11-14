"""WS6: GCxGC class x carbon-number matrix models for real-fuel properties.  Readings, gates and settings are declared in
results/v5/BASE/delta_ws6_gcxgc/README.md (committed before this script existed).

Phases (all write into the delta folder):
  table   -> fuel_table.csv (one row per fuel x cell), features.csv
  base    -> power.csv  (baseline errors only: constant / physics / physics+offset; no ML model fitted)
  models  -> oof_G.csv, oof_T.csv, summary.csv, margins.csv, readings.json, proposed_margins.yaml
Run:  PYTHONPATH=src python scripts/v5/base/delta_ws6_gcxgc.py all
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.uq.conformal import coverage, one_sided_margin  # noqa: E402

BASE = ROOT / "results" / "v5" / "BASE"
OUT = BASE / "delta_ws6_gcxgc"
V5 = ROOT / "data" / "real_fuels" / "v5"
CLASSES = ["n_paraffin", "iso_paraffin", "monocycloparaffin", "dicycloparaffin", "tricycloparaffin",
           "alkylbenzene", "cycloaromatic", "diaromatic", "alkene"]
CN_BINS = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17]  # <=8, 9..16, >=17
AV_MAP = {"n-alkane": "n_paraffin", "iso-alkane": "iso_paraffin", "monocycloalkane": "monocycloparaffin",
          "dicycloalkane": "dicycloparaffin", "tricycloalkane": "tricycloparaffin", "alkylbenzene": "alkylbenzene",
          "cycloaromatic": "cycloaromatic", "diaromatic": "diaromatic", "alkene": "alkene"}
AJ_MAP = {c: c for c in CLASSES} | {"diene": "alkene"}
FEATS = [f"{c}|C{b}" for c in CLASSES for b in CN_BINS]
K_FOLDS, N_REP, B_BOOT = 5, 20, 5000
# cell -> (property label, direction of unsafe error, kind)
CELLS = {
    "density_288": ("density_kg_m3", "two", "abs"),
    "flash": ("flash_point_C", "over", "abs"),
    "visc_253": ("viscosity_kin_cSt", "under", "ln"),
    "visc_233": ("viscosity_kin_cSt", "under", "ln"),
    "freezing": ("freezing_point_C", "under", "abs"),
    "smoke": ("smoke_point_mm", "over", "abs"),
}


def cn_bin(n: float) -> int:
    return int(min(max(n, 8), 17))


def matrix(df: pd.DataFrame, cls_col: str, cn_col: str, w_col: str, cmap: dict) -> tuple[np.ndarray, float, float]:
    v = np.zeros(len(FEATS))
    tot = 0.0
    for r in df.itertuples():
        c = cmap.get(getattr(r, cls_col))
        w = float(getattr(r, w_col))
        tot += w
        if c is None:
            continue
        v[CLASSES.index(c) * len(CN_BINS) + CN_BINS.index(cn_bin(getattr(r, cn_col)))] += w
    mapped = v.sum()
    return (v / mapped if mapped > 0 else v), tot, mapped / tot if tot else 0.0


# ------------------------------------------------------------------ phase table
def build_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    fu = pd.read_csv(BASE / "t2_fuels.csv")
    fu = fu[fu.chosen_for_group & ~fu.low_coverage].copy()
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv"),
             "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv")}
    feats, meta = {}, []
    for r in fu.itertuples():
        if r.source == "av3322":
            d = comps["av3322"]; d = d[d.fuel_id == r.fuel_id]
            x, tot, frac = matrix(d, "group", "carbon_number", "mass_pct", AV_MAP)
        elif r.source == "ajftd":
            d = comps["ajftd"]; d = d[d.fuel_id == int(r.fuel_id)]
            x, tot, frac = matrix(d, "hc_class", "n_carbon", "mass_pct", AJ_MAP)
        else:
            d = comps["fuellib"]; d = d[d.fuel_id == r.fuel_id]
            x, tot, frac = matrix(d, "hc_class", "n_carbon", "mass_fraction_norm", AJ_MAP)
        ok = (98.0 <= tot <= 102.0) if r.source != "fuellib" else abs(tot - 1.0) < 0.02 or (98.0 <= tot <= 102.0)
        meta.append(dict(entity=r.entity, source=r.source, fuel_id=r.fuel_id, group=r.group, family="av3322" if r.source == "av3322" else "ajftd",
                         comp_sum=tot, mapped_frac=frac, comp_ok=bool(ok and frac > 0.97)))
        feats[r.entity] = x
    meta = pd.DataFrame(meta)
    F = pd.DataFrame.from_dict(feats, orient="index", columns=FEATS)
    F.index.name = "entity"
    ok_ent = set(meta[meta.comp_ok].entity)

    tg = pd.read_csv(BASE / "t2_targets.csv")
    bulk = pd.read_csv(BASE / "rows_t2_bulk.csv")
    fl = pd.read_csv(BASE / "rows_t2_flash.csv")
    s5 = pd.read_csv(BASE / "../S5/s5_predictions.csv") if (BASE / "../S5/s5_predictions.csv").exists() else pd.read_csv(ROOT / "results/v5/S5/s5_predictions.csv")
    rows = []

    def add(cell, ent, y, p, note=""):
        if ent in ok_ent and np.isfinite(y) and np.isfinite(p):
            rows.append(dict(cell=cell, entity=ent, y=float(y), p=float(p), note=note))

    b = bulk[bulk.chosen & ~bulk.low_coverage & (bulk.status == "ok")]
    for cell, prop, T, model in [("density_288", "density_kg_m3", 288.15, "ideal_volume"), ("visc_253", "viscosity_kin_cSt", 253.15, "refutas_raw"),
                                 ("visc_233", "viscosity_kin_cSt", 233.15, "refutas_raw")]:
        g = b[(b.property == prop) & (b.T_K.round(2) == T) & (b.model == model) & (b.variant == "designated")]
        for r in g.itertuples():
            if prop == "viscosity_kin_cSt":
                add(cell, r.entity, np.log(r.measured), np.log(r.pred))
            else:
                add(cell, r.entity, r.measured, r.pred)
    f2 = fl[fl.chosen & ~fl.low_coverage & (fl.status == "ok")]
    for r in f2.itertuples():
        add("flash", r.entity, r.measured, r.pred_designated, note=str(r.method_family))
    cens = set()
    fz_t = tg[tg.property == "freezing_point_C"].set_index("entity")
    for r in s5.itertuples():
        if r.entity in fz_t.index and np.isfinite(r.F1x):
            ft = fz_t.loc[r.entity]
            add("freezing", r.entity, float(ft.measured), r.F1x)
    # smoke: physics = smoke_point_from_bulk(aromatics vol %) ; targets: AV-33-22 / FuelLib from t2_targets, AJFTD from ajftd_smoke_naph
    from fuelblend.baselines.sooting import smoke_point_from_bulk
    insample = set(pd.read_csv(ROOT / "data/processed/v4/sooting/smoke_point_ajftd.csv").fuel_id.astype(int))
    avb = pd.read_csv(V5 / "crc_av3322_bulk.csv").set_index("fuel_id")
    sn = pd.read_csv(OUT / "inputs" / "ajftd_smoke_naph.csv")
    sn = sn[sn.smoke_qualifier.isna()].drop_duplicates("fuel_id").set_index("fuel_id")
    st = tg[tg.property == "smoke_point_mm"].set_index("entity")
    excl_insample, oodom = [], []
    for ent in feats:
        s = fu[fu.entity == ent].iloc[0]
        y = ar = np.nan
        note = ""
        if s.source == "av3322" and ent in st.index and s.fuel_id in avb.index:
            y, ar = float(st.loc[ent].measured), float(avb.loc[s.fuel_id].coa_aromatics_vol_pct)
        elif s.source == "ajftd":
            fid = int(s.fuel_id)
            if fid in sn.index:
                if fid in insample:
                    excl_insample.append(ent); continue
                y, ar = float(sn.loc[fid].smoke_point_mm), float(sn.loc[fid].aromatics_vol_pct)
        if np.isfinite(y) and np.isfinite(ar):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p = smoke_point_from_bulk(ar)
            if not (8.0 <= ar <= 23.0):
                oodom.append(ent); note = "aromatics_out_of_domain"
            add("smoke", ent, y, p, note=note)
    t = pd.DataFrame(rows).merge(meta[["entity", "source", "group", "family"]], on="entity")
    t.to_csv(OUT / "fuel_table.csv", index=False)
    F.loc[sorted(ok_ent)].to_csv(OUT / "features.csv")
    log = dict(n_fuels_chosen=len(fu), n_comp_ok=len(ok_ent), comp_not_ok=sorted(set(meta.entity) - ok_ent),
               smoke_ajftd_insample_excluded=excl_insample, smoke_aromatics_out_of_domain=oodom,
               counts={c: int((t.cell == c).sum()) for c in CELLS}, counts_by_family={c: t[t.cell == c].family.value_counts().to_dict() for c in CELLS})
    (OUT / "table_log.json").write_text(json.dumps(log, indent=1))
    return t, F


# ------------------------------------------------------------------ models
def pls_fit(Xtr, rtr, k):
    from sklearn.cross_decomposition import PLSRegression
    m = PLSRegression(n_components=k, scale=True)
    m.fit(Xtr, rtr)
    return m


def choose_k(Xtr, rtr, kmax=5, seed=0):
    n = len(rtr)
    kmax = max(1, min(kmax, n - 3))
    rng = np.random.default_rng(seed)
    folds = rng.permutation(n) % 5
    best, bk = np.inf, 1
    for k in range(1, kmax + 1):
        err = []
        for f in range(5):
            te = folds == f
            if te.sum() == 0 or (~te).sum() <= k + 1:
                continue
            m = pls_fit(Xtr[~te], rtr[~te], k)
            err.append(np.abs(m.predict(Xtr[te]).ravel() - rtr[te]))
        e = np.concatenate(err).mean() if err else np.inf
        if e < best - 1e-12:
            best, bk = e, k
    return bk


def fit_predict(Xtr, ytr, ptr, Xte, pte, models):
    """Returns {model: prediction of y on test rows}.  Everything fitted on the training rows only."""
    out = {"C0": np.full(len(Xte), ytr.mean()), "C1": pte.copy()}
    r = ytr - ptr
    mu = r.mean()
    out["C2"] = pte + mu
    if not models:
        return out
    rc = r - mu
    k = choose_k(Xtr, rc)
    pls = pls_fit(Xtr, rc, k)
    out["M1"] = pte + mu + pls.predict(Xte).ravel()
    Str, Ste = pls.transform(Xtr), pls.transform(Xte)
    sd = Str.std(0); sd[sd == 0] = 1.0
    Str, Ste = Str / sd, Ste / sd
    s_r = rc.std() if rc.std() > 0 else 1.0
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
    kern = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(np.ones(k), (0.1, 100.0)) + WhiteKernel(0.5, (1e-3, 5.0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gp = GaussianProcessRegressor(kern, n_restarts_optimizer=2, random_state=0).fit(Str, rc / s_r)
        out["M2"] = pte + mu + gp.predict(Ste) * s_r
        from sklearn.neural_network import MLPRegressor
        pr = np.zeros(len(Xte))
        for sd_ in range(10):
            nn = MLPRegressor(hidden_layer_sizes=(8,), activation="tanh", alpha=1.0, solver="lbfgs", max_iter=500, random_state=sd_).fit(Str, rc / s_r)
            pr += nn.predict(Ste) * s_r / 10
    out["M3"] = pte + mu + pr
    k0 = choose_k(Xtr, ytr - ytr.mean())
    m0 = pls_fit(Xtr, ytr - ytr.mean(), k0)
    out["M0"] = ytr.mean() + m0.predict(Xte).ravel()
    return out


def run_G(t: pd.DataFrame, F: pd.DataFrame, models: bool, n_rep=N_REP) -> pd.DataFrame:
    rows = []
    for cell in CELLS:
        d = t[t.cell == cell].reset_index(drop=True)
        X = F.loc[d.entity].to_numpy(float)
        y, p = d.y.to_numpy(float), d.p.to_numpy(float)
        grp = d.group.to_numpy()
        ug = np.unique(grp)
        for rep in range(n_rep):
            rng = np.random.default_rng(1000 + rep)
            fold_of = dict(zip(rng.permutation(ug), np.arange(len(ug)) % K_FOLDS))
            fo = np.array([fold_of[g] for g in grp])
            for f in range(K_FOLDS):
                te = fo == f
                if te.sum() == 0 or (~te).sum() < 8:
                    continue
                pr = fit_predict(X[~te], y[~te], p[~te], X[te], p[te], models)
                for mname, v in pr.items():
                    for i, vi in zip(np.where(te)[0], v):
                        rows.append((cell, d.entity[i], d.group[i], d.family[i], rep, mname, vi - y[i]))
    return pd.DataFrame(rows, columns=["cell", "entity", "group", "family", "rep", "model", "err"])


def run_T(t: pd.DataFrame, F: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell in CELLS:
        d = t[t.cell == cell].reset_index(drop=True)
        X = F.loc[d.entity].to_numpy(float)
        y, p = d.y.to_numpy(float), d.p.to_numpy(float)
        for src in ("av3322", "ajftd"):
            te = (d.family == src).to_numpy()
            tr = ~te
            if te.sum() < 3 or tr.sum() < 8:
                continue
            # training rows whose group appears in the test source are removed
            tr = tr & ~d.group.isin(set(d.group[te])).to_numpy()
            pr = fit_predict(X[tr], y[tr], p[tr], X[te], p[te], True)
            for mname, v in pr.items():
                for i, vi in zip(np.where(te)[0], v):
                    rows.append((cell, d.entity[i], d.group[i], d.family[i], 0, mname, vi - y[i]))
    return pd.DataFrame(rows, columns=["cell", "entity", "group", "family", "rep", "model", "err"])


# ------------------------------------------------------------------ statistics
def per_fuel(o: pd.DataFrame) -> pd.DataFrame:
    o = o.assign(ae=o.err.abs())
    return o.groupby(["cell", "model", "entity", "group", "family"], as_index=False).agg(ae=("ae", "mean"), err=("err", "mean"))


def boot(d: np.ndarray, B=B_BOOT, seed=0):
    rng = np.random.default_rng(seed)
    n = len(d)
    idx = rng.integers(0, n, size=(B, n))
    s = d[idx].mean(1)
    return float(d.mean()), float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5)), float((s <= 0).mean()), s


def power_table(pf: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell in CELLS:
        a = pf[(pf.cell == cell) & (pf.model == "C2")].set_index("entity").ae
        c = pf[(pf.cell == cell) & (pf.model == "C0")].set_index("entity").ae
        dd = (c - a.reindex(c.index)).dropna()
        n = len(dd)
        rows.append(dict(cell=cell, n=n, mae_C0=c.mean(), mae_C2=a.mean(), sd_paired_diff=dd.std(), mdd_80pct_power=2.8 * dd.std() / np.sqrt(n) if n > 1 else np.nan))
    return pd.DataFrame(rows)


def unsafe_flags(cell, err, margin=None):
    direction = CELLS[cell][1]
    if direction == "over":
        return err > (5.0 if cell == "flash" else (margin if margin is not None else 0.0))
    if direction == "under":
        return -err > (margin if margin is not None else 0.0)
    return np.zeros(len(err), bool)


def summarise(pf: pd.DataFrame, secondary_alpha=0.017) -> tuple[pd.DataFrame, dict]:
    rows, readings = [], {}
    for cell, (prop, direction, _) in CELLS.items():
        g = pf[pf.cell == cell]
        models = [m for m in ["C0", "C1", "C2", "M0", "M1", "M2", "M3"] if m in set(g.model)]
        W = {m: g[g.model == m].set_index("entity") for m in models}
        ents = W["C2"].index
        n = len(ents)
        for m in models:
            e, ae = W[m].err.reindex(ents), W[m].ae.reindex(ents)
            row = dict(cell=cell, model=m, n=n, mae=ae.mean(), bias=e.mean())
            if direction in ("over", "under"):
                mg = one_sided_margin(e.to_numpy(), W[m].group.reindex(ents).to_numpy(), 0.10, direction)
                row["margin90_oof"] = mg
                row["unsafe_rate_at_own_margin"] = float(unsafe_flags(cell, e.to_numpy(), mg).mean()) if np.isfinite(mg) else np.nan
                if cell == "flash":
                    row["over_gt5C_rate"] = float((e > 5).mean()); row["worst_over"] = float(e.max())
                elif cell.startswith("visc"):
                    row["under_gt_0p4154_rate"] = float((-e > 0.4154).mean()); row["worst_under"] = float((-e).max())
                elif cell in ("freezing", "smoke"):
                    row["unsafe_rate_at_zero"] = float(unsafe_flags(cell, e.to_numpy(), 0.0).mean())
            if m.startswith("M"):
                for ref in ("C0", "C1", "C2"):
                    d = (W[ref].ae.reindex(ents) - ae).to_numpy()
                    diff, lo, hi, p1, s = boot(d)
                    row[f"d_vs_{ref}"], row[f"lo_{ref}"], row[f"hi_{ref}"], row[f"p1_{ref}"] = diff, lo, hi, p1
                    if m != "M1":
                        lo2, hi2 = np.percentile(boot(d)[4], [100 * secondary_alpha / 2, 100 * (1 - secondary_alpha / 2)])
                        row[f"lo98_{ref}"], row[f"hi98_{ref}"] = lo2, hi2
            rows.append(row)
        readings[cell] = dict(n=int(n))
    return pd.DataFrame(rows), readings


def holm(p: dict[str, float], alpha=0.05) -> dict[str, bool]:
    order = sorted(p, key=p.get)
    ok, passed = True, {}
    for i, k in enumerate(order):
        ok = ok and p[k] <= alpha / (len(order) - i)
        passed[k] = ok
    return passed


def main(phase="all") -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t, F = build_table() if phase in ("table", "all") else (pd.read_csv(OUT / "fuel_table.csv"), pd.read_csv(OUT / "features.csv", index_col="entity"))
    if phase == "table":
        return
    if phase in ("base", "all"):
        oB = run_G(t, F, models=False)
        pfb = per_fuel(oB)
        pw = power_table(pfb)
        pw.to_csv(OUT / "power.csv", index=False)
        print(pw.to_string())
        if phase == "base":
            return
    if phase == "summary":
        oG, oT = pd.read_csv(OUT / "oof_G.csv"), pd.read_csv(OUT / "oof_T.csv")
    else:
        oG = run_G(t, F, models=True)
        oG.to_csv(OUT / "oof_G.csv", index=False)
        oT = run_T(t, F)
        oT.to_csv(OUT / "oof_T.csv", index=False)
    pf = per_fuel(oG)
    sG, rd = summarise(pf)
    sG.to_csv(OUT / "summary_G.csv", index=False)
    pfT = per_fuel(oT)
    trows = []
    for (cell, fam), g in pfT.groupby(["cell", "family"]):
        for m, gg in g.groupby("model"):
            trows.append(dict(cell=cell, test_family=fam, model=m, n=len(gg), mae=gg.ae.mean(), bias=gg.err.mean()))
    pd.DataFrame(trows).to_csv(OUT / "summary_T.csv", index=False)

    # readings
    res = {}
    pvals = {}
    for cell, (prop, direction, _) in CELLS.items():
        s = sG[(sG.cell == cell)].set_index("model")
        m1 = s.loc["M1"]
        ci_ok = all(m1[f"lo_{r}"] > 0 for r in ("C0", "C1", "C2"))
        safe_ok = True
        if direction in ("over", "under"):
            key = "over_gt5C_rate" if cell == "flash" else "unsafe_rate_at_own_margin"
            safe_ok = bool(m1[key] <= s.loc["C2"][key] + 1e-12)
        res[cell] = dict(n=int(m1.n), mae={m: float(s.loc[m].mae) for m in s.index}, M1_vs={r: [float(m1[f"d_vs_{r}"]), float(m1[f"lo_{r}"]), float(m1[f"hi_{r}"])] for r in ("C0", "C1", "C2")},
                         intervals_above_0=bool(ci_ok), unsafe_not_higher=safe_ok, gate_pass_unadjusted=bool(ci_ok and safe_ok),
                         C2_beats_C1=bool(s.loc["C2"].mae < s.loc["C1"].mae))
        pvals[cell] = float(max(m1[f"p1_{r}"] for r in ("C0", "C1", "C2")))
    h = holm(pvals)
    for c in res:
        res[c]["weakest_p"] = pvals[c]; res[c]["holm_survives"] = bool(h[c]); res[c]["gate_established"] = bool(res[c]["gate_pass_unadjusted"] and h[c])
    # secondary models: 98.3 % intervals, reported
    sec = {}
    for cell in CELLS:
        s = sG[sG.cell == cell].set_index("model")
        for m in ("M2", "M3"):
            sec[f"{cell}/{m}"] = bool(all(s.loc[m][f"lo98_{r}"] > 0 for r in ("C0", "C1", "C2")))
    # margins & transfer coverage (R-WS6-3, R-WS6-4)
    mrows = []
    for cell in ("visc_253", "visc_233", "freezing"):
        for m in ("C1", "C2", "M1"):
            g = pf[(pf.cell == cell) & (pf.model == m)]
            direction = CELLS[cell][1]
            e = g.err.to_numpy()
            mall = one_sided_margin(e, g.group.to_numpy(), 0.10, direction)
            rr = dict(cell=cell, model=m, n_groups=len(g), margin90_all=mall)
            for a, b in (("av3322", "ajftd"), ("ajftd", "av3322")):
                ga, gb = g[g.family == a], g[g.family == b]
                ma = one_sided_margin(ga.err.to_numpy(), ga.group.to_numpy(), 0.10, direction) if len(ga) else np.inf
                if np.isfinite(ma) and len(gb):
                    k, nn = coverage(gb.err.to_numpy(), gb.group.to_numpy(), ma, direction)
                    rr[f"margin_{a}"] = ma; rr[f"cover_{b}_by_{a}"] = k / nn; rr[f"n_{b}"] = nn
                else:
                    rr[f"margin_{a}"] = ma; rr[f"cover_{b}_by_{a}"] = np.nan; rr[f"n_{b}"] = len(gb)
            mrows.append(rr)
    mg = pd.DataFrame(mrows)
    mg.to_csv(OUT / "margins.csv", index=False)
    r3 = {}
    for cell in ("visc_253", "visc_233"):
        for m in ("C2", "M1"):
            q = mg[(mg.cell == cell) & (mg.model == m)].iloc[0]
            cov = [q["cover_ajftd_by_av3322"], q["cover_av3322_by_ajftd"]]
            r3[f"{cell}/{m}"] = dict(margin_all=float(q.margin90_all), margin_av3322=float(q.margin_av3322), margin_ajftd=float(q.margin_ajftd), coverage=cov,
                                     holds=bool(np.isfinite(q.margin90_all) and all(np.isfinite(c) and c >= 0.90 for c in cov)))
    (OUT / "readings.json").write_text(json.dumps(dict(R1_R2=res, secondary_98pct_pass=sec, R3_cold_visc_margin=r3), indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))
    print(mg.to_string())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
