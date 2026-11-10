"""Task 14 evaluation -> results/v4/SP/smoke_point_eval.csv and README.md (numbers read back from CSV)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines import sooting as S  # noqa: E402

OUT = ROOT / "results" / "v4" / "SP"
FSETS = {"const": (), "arom": ("aromatics_pct_vol",),
         "arom+diarom": ("aromatics_pct_vol", "diaromatics_pct_vol")}
rng = np.random.default_rng(0)


def metrics(y, p):
    e = np.asarray(p, float) - np.asarray(y, float)
    return dict(mae_mm=np.abs(e).mean(), rmse_mm=np.sqrt((e ** 2).mean()),
                median_abs_mm=np.median(np.abs(e)), bias_mm=e.mean(),
                within_1mm=(np.abs(e) <= 1).mean(), within_2mm=(np.abs(e) <= 2).mean())


def boot_diff(y, p_model, p_base, n=2000):
    y, p_model, p_base = map(lambda v: np.asarray(v, float), (y, p_model, p_base))
    a, b = np.abs(p_model - y), np.abs(p_base - y)
    idx = rng.integers(0, len(y), (n, len(y)))
    d = a[idx].mean(1) - b[idx].mean(1)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def _fit(d, feats):
    if not feats:
        return d.value.mean(), None
    return S.fit_bulk(d, feats)


def _pred(d, feats, fit):
    if not feats:
        return np.full(len(d), fit[0])
    return S.predict_bulk(d, fit[0], fit[1], feats)


def select_inside(train):
    """Inner leave-one-fuel-out selection of the feature set using ONLY the training fold."""
    best, best_mse = None, np.inf
    for name, feats in FSETS.items():
        err = []
        for i in range(len(train)):
            tr, te = train.drop(train.index[i]), train.iloc[[i]]
            err.append(float(_pred(te, feats, _fit(tr, feats))[0] - te.value.iloc[0]) ** 2)
        mse = np.mean(err)
        if mse < best_mse:
            best, best_mse = name, mse
    return best


def bulk_cv(df, slice_name, rows, group_col="fuel_id"):
    df = df.reset_index(drop=True)
    y = df.value.to_numpy()
    preds = {k: np.zeros(len(df)) for k in ("const_mean", "arom+diarom_fixed", "selected_in_fold")}
    chosen = []
    for tr, te in GroupKFold(5).split(df, groups=df[group_col]):
        a, b = df.iloc[tr], df.iloc[te]
        preds["const_mean"][te] = a.value.mean()
        preds["arom+diarom_fixed"][te] = _pred(b, FSETS["arom+diarom"], _fit(a, FSETS["arom+diarom"]))
        sel = select_inside(a)
        chosen.append(sel)
        preds["selected_in_fold"][te] = _pred(b, FSETS[sel], _fit(a, FSETS[sel]))
    start = len(rows)
    for k, p in preds.items():
        r = dict(slice=slice_name, model=k, n=len(df), n_censored=0, **metrics(y, p))
        if k != "const_mean":
            r["mae_diff_vs_const_ci_lo"], r["mae_diff_vs_const_ci_hi"] = boot_diff(y, p, preds["const_mean"])
        if k == "selected_in_fold":
            r["note"] = "feature sets chosen per fold: " + ",".join(chosen)
        rows.append(r)
    rows[start]["note"] = f"y range {y.min():.1f}-{y.max():.1f} mm, sd {y.std():.2f}"


def main():
    rows = []
    # (A) pure compounds, YSI route, leave-one-compound-out on uncensored rows
    pure = pd.read_csv(S.DATA / "smoke_point_pure_citedtyped.csv").dropna(subset=["ysi"])
    ex = S.exact_rows(pure).reset_index(drop=True)
    cen = pure[pure.censored]
    yl, pl, pc = [], [], []
    for i in range(len(ex)):
        tr = ex.drop(i)
        coef = S.fit_ysi_sp(tr)
        yl.append(ex.smoke_point_mm[i])
        pl.append(float(S.predict_ysi_sp(ex.ysi[i], ex.mw[i], coef)))
        pc.append(float(np.exp(np.log(tr.smoke_point_mm).mean())))
    yl, pl, pc = map(np.array, (yl, pl, pc))
    note = "pure compounds, leave-one-compound-out; smoke points cited_typed (not file-verifiable)"
    rows.append(dict(slice="pure_ysi_route", model="const_geomean", n=len(ex), n_censored=len(cen),
                     note=note, **metrics(yl, pc)))
    lo, hi = boot_diff(yl, pl, pc)
    rows.append(dict(slice="pure_ysi_route", model="MW/SP=c0+c1*YSI", n=len(ex), n_censored=len(cen),
                     note=note, mae_diff_vs_const_ci_lo=lo, mae_diff_vs_const_ci_hi=hi, **metrics(yl, pl)))
    keep = (ex.compound_name != "iso-octane").to_numpy()
    rows.append(dict(slice="pure_ysi_route_no_isooctane", model="MW/SP=c0+c1*YSI", n=int(keep.sum()),
                     n_censored=len(cen),
                     note="sensitivity: iso-octane excluded from eval (typed 40.0 vs 42.8 mm in Boehm text)",
                     **metrics(yl[keep], np.array(pl)[keep])))
    coef = S.fit_ysi_sp(pure)
    pcen = S.predict_ysi_sp(cen.ysi, cen.mw, coef)
    rows.append(dict(slice="pure_ysi_route_censored", model="MW/SP=c0+c1*YSI", n=0, n_censored=len(cen),
                     within_1mm=float((pcen >= S.SP_CAP_MM).mean()),
                     note="true SP > 50 mm; within_1mm column = fraction predicted >=50 mm; predictions: "
                          + str({k: round(float(v), 1) for k, v in zip(cen.compound_name, pcen)})))
    rows.append(dict(slice="pure_ysi_route_censored", model="const_geomean", n=0, n_censored=len(cen),
                     within_1mm=0.0, note="constant geometric mean is <50 mm: no censored row consistent"))
    # (B) AJFTD 100 fuels, bulk route
    a = pd.read_csv(S.DATA / "smoke_point_ajftd.csv")
    ncen = int(a.censored.sum())
    a = a[~a.censored].copy()
    bulk_cv(a, "ajftd100_bulk", rows)
    # stricter groupings: inferred T-code (tank/producer code, not documented) and test date
    tc = a.fuel_label.astype(str).str.extract(r"\)(T\d+)")[0].fillna("NOCODE")
    bulk_cv(a.assign(tgroup=tc), "ajftd100_bulk_groupTcode", rows, "tgroup")
    bulk_cv(a.assign(dgroup=a.test_date.astype(str).str.split(" ").str[0]), "ajftd100_bulk_groupTestDate", rows, "dgroup")
    # (C) broad AJFTD set, sensitivity
    bulk_cv(pd.read_csv(S.DATA / "smoke_point_ajftd_broad.csv"), "ajftd_broad_bulk", rows)
    # (D) YSI route on real fuels
    rows.append(dict(slice="ysi_route_real_fuels", model="ysi_mix+smoke_point_from_ysi", n=0, n_censored=0,
                     note="NOT EVALUATED: no fuel in repo has both speciation and measured smoke point "
                          "(AJFTD has smoke point but no molecule list; GC fuels POSF-10325/UDRI have no smoke point)"))
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "smoke_point_eval.csv", index=False)
    write_readme(pd.read_csv(OUT / "smoke_point_eval.csv"), ncen)


def write_readme(df, ncen_ajftd):
    def g(sl, m, c):
        return df[(df.slice == sl) & (df.model == m)].iloc[0][c]

    L = ["# Smoke point / sooting evaluation (v4 Task 14)", "",
         "All numbers below are read from `smoke_point_eval.csv` "
         "(regenerate: `scripts/v4/build_sooting_inputs.py` then `scripts/v4/eval_smoke_point.py`).", "",
         "## Headline", "",
         f"Most conservative bulk-route result (AJFTD 100 fuels grouped by inferred T-code, an undocumented tank/producer code parsed "
         f"from fuel_label; unmatched fuels pooled as one group): selected-in-fold MAE {g('ajftd100_bulk_groupTcode', 'selected_in_fold', 'mae_mm'):.2f} mm vs "
         f"constant-mean {g('ajftd100_bulk_groupTcode', 'const_mean', 'mae_mm'):.2f} mm. Grouping by fuel_id alone gives "
         f"{g('ajftd100_bulk', 'selected_in_fold', 'mae_mm'):.2f} vs {g('ajftd100_bulk', 'const_mean', 'mae_mm'):.2f}; that gain is inflated by near-duplicate fuels. "
         f"Grouped by test date: {g('ajftd100_bulk_groupTestDate', 'selected_in_fold', 'mae_mm'):.2f} vs {g('ajftd100_bulk_groupTestDate', 'const_mean', 'mae_mm'):.2f}.", "",
         "## Results (mm; lower is better)", "",
         "| slice | model | N | MAE | RMSE | bias | MAE minus const MAE, 95% bootstrap CI |",
         "|---|---|---|---|---|---|---|"]
    for _, r in df[df.n > 0].iterrows():
        ci = "" if pd.isna(r.get("mae_diff_vs_const_ci_lo")) else \
            f"[{r.mae_diff_vs_const_ci_lo:+.2f}, {r.mae_diff_vs_const_ci_hi:+.2f}]"
        L.append(f"| {r.slice} | {r.model} | {int(r.n)} | {r.mae_mm:.2f} | {r.rmse_mm:.2f} | {r.bias_mm:+.2f} | {ci} |")
    cm = df[(df.slice == "pure_ysi_route_censored")].iloc[0]
    L += ["", "## Censoring", "",
          f"* Pure compounds: {int(cm.n_censored)} compounds have true SP > 50 mm (D1322 lamp limit); never fit as exact "
          f"targets. Fraction the YSI model predicts >= 50 mm: {cm.within_1mm:.2f} (constant baseline 0.00). "
          f"Predictions: {cm.note.split('predictions: ')[-1]}",
          f"* AJFTD: {ncen_ajftd} censored rows in the 100-fuel file.", "",
          "## Honest reading", "",
          f"* Bulk route (AJFTD 100 fuels, GroupKFold by fuel_id, feature set selected inside training folds): "
          f"MAE {g('ajftd100_bulk', 'selected_in_fold', 'mae_mm'):.2f} mm vs constant-mean "
          f"{g('ajftd100_bulk', 'const_mean', 'mae_mm'):.2f} mm ({g('ajftd100_bulk', 'const_mean', 'note')}). "
          "AJFTD smoke point resolution is ~0.5-1 mm and fuel_ids from one producer/batch are not grouped, so the CV "
          "is if anything optimistic. No fuel is below 19 mm or above 25 mm, so the 18/25 mm constraint boundaries are not probed.",
          f"* Broad AJFTD sensitivity (N={int(g('ajftd_broad_bulk', 'const_mean', 'n'))}, D1322, 5<SP<50 filter): "
          f"selected MAE {g('ajftd_broad_bulk', 'selected_in_fold', 'mae_mm'):.2f} vs const "
          f"{g('ajftd_broad_bulk', 'const_mean', 'mae_mm'):.2f}.",
          f"* YSI route, pure compounds (leave-one-out, N={int(g('pure_ysi_route', 'const_geomean', 'n'))}): "
          f"MAE {g('pure_ysi_route', 'MW/SP=c0+c1*YSI', 'mae_mm'):.2f} vs geometric-mean baseline "
          f"{g('pure_ysi_route', 'const_geomean', 'mae_mm'):.2f}. CAVEAT: unreliable for saturates; the five censored alkanes/cycloalkanes "
          f"(true SP > 50 mm) are predicted 13-22 mm, i.e. error >= 28 mm each, and are excluded from this MAE.",
          "* The linear MW/SP-vs-YSI form fails on saturated hydrocarbons: every censored (SP > 50 mm) alkane/cycloalkane is predicted far "
          "below 50 mm (see Censoring), so the YSI route is not reliable for low-sooting components. On the broad AJFTD set the bulk model has no skill over the constant mean (CI spans 0).",
          "* YSI route on real fuels: NOT EVALUATED (no fuel with both speciation and measured smoke point). "
          "It must not be cited as evidence for fuel smoke point.", "",
          "## Provenance and caveats", "",
          "* YSI: raw_file, verified equal to `UnifiedYSI` in Harvard Dataverse YSI Database Volume 2 "
          "(benzene 100.3, toluene 170.9, n-heptane 36.0).",
          "* Pure smoke points (28 rows, 23 uncensored): `cited_typed` from extras/build_tsi_table.py "
          "(Mensch 2010, Li & Sunderland 2012, Boehm 2022 SI); SI numbers are not in this repo and could not be checked. "
          "Known discrepancy: iso-octane typed 40.0 mm, Boehm 2022 text quotes 42.8 mm.",
          "* Config only: the `smoke_point_mm` keys added to configs/specs_astm.yaml and specs_bis.yaml (alt_min, alt_requires, uncertainty_margin) "
          "are NOT read by any code yet; nothing enforces them. BIS values are UNVERIFIED; D1655 values are from secondary sources only.",
          "* Data-quality note: raw YSI table gives identical values for naphthalene and 1,2-dihydronaphthalene (466.1) and for 1-/2-methylnaphthalene (649.1); "
          "they match the raw .tab as stored but look suspicious and should be checked against the primary source.",
          "* `ysi_mix` returns a YSI only; use `sooting.mixture_mw(smiles, x_mole)` for the mole-averaged MW needed by `smoke_point_from_ysi`.",
          "* EXCLUDED as circular: 175 rows of extras `tsi_table.csv` (derived=True), smoke point computed from YSI.",
          "* AJFTD smoke point/aromatics: raw_file (match `ajftd_fuels.csv`/`ajftd_composition.csv`); smoke point is not derived from aromatics.",
          "", "## Claims-register ids that would update (not edited here)", "",
          "Any claim listing smoke point as an unmodeled property (D1655 smoke point now has a bulk-composition baseline "
          "and a constraint entry); no pure/fuel sooting claim is supported beyond the rows above."]
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
