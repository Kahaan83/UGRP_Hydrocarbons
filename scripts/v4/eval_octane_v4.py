"""Task 15 evaluation -> results/v4/O/octane_eval_v4.csv and README.md (numbers read back from the CSV).

Targets: Co-Optima bRON blends only (ruling R17: per-property target_ok_ron / target_ok_mon, base RON/MON
as inputs; FACE rows are never targets).  CV: the stored v4 dev folds (blendstock-CAS groups, R19) plus
leave-one-blendstock-out; the locked set (4 CAS groups) is a separate thin test fit on all dev rows.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.models import octane_blend_v4 as O  # noqa: E402

CORPUS = ROOT / "data" / "real_fuels" / "v4" / "gasoline_octane_blends.csv"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
OUT = ROOT / "results" / "v4" / "O"
REFS = ["base_unchanged", "constant_mean", "linear_volume_pooled_bov"]


def run_view(table, folds, prop):
    """folds: list of (train_idx, test_idx).  Fresh models per fold; returns predictions DataFrame."""
    pred = {}
    idx_all = []
    for tr, te in folds:
        train, test = table.iloc[tr], table.iloc[te]
        for m in O.make_models(prop):
            pred.setdefault(m.name, {}).update(dict(zip(te, O.fit_predict(m, train, test))))
        idx_all += list(te)
    idx = sorted(idx_all)
    return pd.DataFrame({k: pd.Series(v).reindex(idx) for k, v in pred.items()}), idx


def rows_for(view, prop, table, P, idx, extra_note=""):
    y = table["y"].to_numpy()[idx]
    g = table["group"].to_numpy()[idx]
    common = np.all(np.isfinite(P.to_numpy()), axis=1)
    out = []
    err = {m: np.abs(P[m].to_numpy() - y) for m in P.columns}
    for m in P.columns:
        mt = O.metrics(y, P[m].to_numpy())
        pt, lo, hi, ng = O.cluster_boot(g, err[m])
        _, clo, chi, _ = O.cluster_boot(g[common], err[m][common]) if common.sum() else (0, np.nan, np.nan, 0)
        cm = O.metrics(y[common], P[m].to_numpy()[common])
        out.append(dict(record_type="metric", view=view, prop=prop, model=m, ref_model="", n=mt["n"], n_groups=ng,
                        mae=mt["mae"], rmse=mt["rmse"], bias=mt["bias"], ci_lo=lo, ci_hi=hi,
                        n_common=cm["n"], mae_common=cm["mae"], ci_lo_common=clo, ci_hi_common=chi,
                        excludes_zero="", note=extra_note))
        for r in REFS:
            if r == m or r not in P.columns:
                continue
            d, dlo, dhi, ngd = O.cluster_boot(g, err[m], err[r])
            out.append(dict(record_type="paired_diff", view=view, prop=prop, model=m, ref_model=r,
                            n=int((np.isfinite(err[m]) & np.isfinite(err[r])).sum()), n_groups=ngd,
                            mae=d, ci_lo=dlo, ci_hi=dhi, excludes_zero=bool(dlo > 0 or dhi < 0),
                            note="mae column = MAE(model) - MAE(ref); negative = model better" + ("; locked has 4 groups, CI NOT reliable, do not cite" if view.startswith("locked") else "")))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    corpus = pd.read_csv(CORPUS)
    sp = json.loads(SPLITS.read_text())
    locked = set(sp["locked"])
    fold_of = {b: k for k, v in sp["folds"].items() for b in v}
    out, reports = [], {}
    for prop in ("ron", "mon"):
        table, rep = O.prepare_blends(corpus, prop)
        reports[prop] = rep
        table["is_locked"] = table["blend_id"].isin(locked)
        table["fold"] = table["blend_id"].map(fold_of)
        dev = table[~table.is_locked].reset_index(drop=True)
        assert dev["fold"].notna().all()
        lock = table[table.is_locked].reset_index(drop=True)
        reports[prop].update(n_dev=len(dev), n_locked=len(lock), groups_dev=dev.group.nunique(),
                             groups_locked=lock.group.nunique())
        # (1) stored dev folds: blendstock-grouped (R19)
        folds = [(np.where(dev.fold != k)[0], np.where(dev.fold == k)[0]) for k in sorted(dev.fold.unique())]
        P, idx = run_view(dev, folds, prop)
        out += rows_for("dev_cv5_stored_folds", prop, dev, P, idx)
        # (2) leave-one-blendstock-out (leave-component-out)
        folds = list(LeaveOneGroupOut().split(dev, groups=dev.group))
        P, idx = run_view(dev, folds, prop)
        out += rows_for("dev_leave_one_blendstock_out", prop, dev, P, idx)
        # (3) locked: fit on all dev rows, thin test (4 CAS groups) -> point estimates only
        P, idx = run_view(pd.concat([dev, lock], ignore_index=True),
                          [(np.arange(len(dev)), np.arange(len(dev), len(dev) + len(lock)))], prop)
        rows = rows_for("locked_thin_test", prop, lock, P.reset_index(drop=True), list(range(len(lock))),
                        "locked: 4 CAS groups, no CI claim possible")
        out += rows
    df = pd.DataFrame(out)
    df.to_csv(OUT / "octane_eval_v4.csv", index=False, float_format="%.6g")
    (OUT / "prep_report.json").write_text(json.dumps(reports, indent=1))
    write_readme(pd.read_csv(OUT / "octane_eval_v4.csv"), reports, corpus)
    print("wrote", OUT / "octane_eval_v4.csv", len(df))


def fmt(v, d=2):
    return "nan" if pd.isna(v) else f"{v:.{d}f}"


def write_readme(df, reports, corpus):
    L = ["# Octane refit on real blends (v4 Task 15)", "",
         "Generated by `scripts/v4/eval_octane_v4.py`; every number below is read back from `octane_eval_v4.csv`.", "",
         "## Scope and honest caveats", "",
         "* Targets are **Co-Optima bRON blends only** (a base gasoline lump of measured RON/MON blended with one pure "
         "blendstock at a measured volume fraction). Ruling R17: base RON/MON enter as model **inputs**; the per-property flags "
         "`target_ok_ron` / `target_ok_mon` are used (not the OR column `target_ok`). FACE rows (144) are never targets "
         "(no fractions in the raw files) and are not evaluated.",
         "* Ruling R19: groups are blendstock CAS numbers, so every CV fold and the locked set test **unseen blendstocks**; "
         "the stored dev folds in `splits/v4/real_fuels_v4_splits.json` are reused. Hyper-parameters (Ghosh ridge alpha, GBM size) are "
         "selected by grouped CV inside each training fold (R6); locked rows are never used for selection.",
         "* Only two distinct base gasolines exist (RON 90.3 / 87.7). Base-fuel inputs therefore carry almost no information; "
         "all variation across blendstocks must be explained by blendstock descriptors (lump class, MW, density, O mass fraction, "
         "unsaturation, ring, aromatic) from 18 dev blendstock groups (16 pure + 2 named mixtures). The blendstock's own bRON cannot be an input (it is the target).",
         "* The locked set is thin: 12 RON / 12 MON blends from 4 CAS groups (ethyl butanoate, ethyl acetate, 2-pentanol, ethanol). "
         "No confidence statements are made from it.",
         "* The linear-by-mole model assumes a base gasoline of M = %.0f g/mol, rho = %.0f kg/m3 (typed typical values, not measured) "
         "and needs blendstock MW and density, so it is NaN for the two named mixtures (ketone, furan). Per-model N is shown; `*_common` "
         "columns restrict to rows every model predicts." % (O.BASE_MW_ASSUMED, O.BASE_RHO_ASSUMED),
         "* Lump priors for alcohol/ester/ketone/other oxygenate are the repo's `tail` lump (82/78) because the repo has no "
         "literature prior for them; they are learned (ridge-shrunk) by the Ghosh fit path.", ""]
    L += ["## Row accounting", "", "| prop | rows | flag false | flag true | dropped frac=None | dropped shape | dropped NaN | used | dev | locked | dev groups | locked groups |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for p, r in reports.items():
        L.append(f"| {p.upper()} | {r['n_rows']} | {r['n_flag_false']} | {r['n_flag_true']} | {r['n_dropped_frac_none']} | "
                 f"{r['n_dropped_shape']} | {r['n_dropped_nan']} | {r['n_used']} | {r['n_dev']} | {r['n_locked']} | "
                 f"{r['groups_dev']} | {r['groups_locked']} |")
    nfn = sum(any(c["frac"] is None for c in json.loads(x)) for x in corpus.components)
    L += ["", f"Rows in the corpus with a `frac=None` component: {nfn} (all FACE `TPRF 4+N`; handled by dropping, never imputed; "
          "none of them is a target row).", ""]
    m = df[df.record_type == "metric"]
    for view, title in (("dev_cv5_stored_folds", "Dev, stored 5 folds (blendstock-grouped)"),
                        ("dev_leave_one_blendstock_out", "Dev, leave-one-blendstock-out"),
                        ("locked_thin_test", "Locked (thin; point estimates, fit on all dev)")):
        for p in ("ron", "mon"):
            s = m[(m.view == view) & (m["prop"] == p)]
            L += [f"## {title} -- {p.upper()}", "",
                  "| model | N | groups | MAE | 95% cluster CI | RMSE | bias | N common | MAE common |", "|---|---|---|---|---|---|---|---|---|"]
            for r in s.itertuples():
                ci = "n/a" if view == "locked_thin_test" else f"[{fmt(r.ci_lo)}, {fmt(r.ci_hi)}]"
                L.append(f"| {r.model} | {r.n} | {r.n_groups} | {fmt(r.mae)} | {ci} | {fmt(r.rmse)} | {fmt(r.bias)} | "
                         f"{int(r.n_common)} | {fmt(r.mae_common)} |")
            L.append("")
    L += ["## Paired differences (MAE(model) - MAE(reference), cluster bootstrap over blendstocks; negative = model better)", ""]
    d = df[(df.record_type == "paired_diff") & (df.view != "locked_thin_test")]
    L += ["| view | prop | model | ref | dMAE | 95% CI | CI excludes 0 |", "|---|---|---|---|---|---|---|"]
    for r in d.itertuples():
        L.append(f"| {r.view} | {r.prop.upper()} | {r.model} | {r.ref_model} | {fmt(r.mae)} | [{fmt(r.ci_lo)}, {fmt(r.ci_hi)}] | {r.excludes_zero} |")
    L += ["", "## Reading the result", ""]
    for view in ("dev_cv5_stored_folds", "dev_leave_one_blendstock_out"):
        for p in ("ron", "mon"):
            s = m[(m.view == view) & (m["prop"] == p)].set_index("model")
            best = s.mae.idxmin()
            wins = d[(d.view == view) & (d["prop"] == p) & (d.model == best) & d.excludes_zero.astype(bool) & (d.mae < 0)]
            L.append(f"* {view} {p.upper()}: lowest MAE = `{best}` ({fmt(s.loc[best, 'mae'])}); references it beats with CI excluding 0: "
                     f"{', '.join(wins.ref_model) if len(wins) else 'none'}.")
    L += ["", "A model is called a winner over a reference only where the paired CI excludes 0 (above). With 18 dev blendstock groups "
          "such CIs are wide; absence of a CI-supported win means the models are not distinguished, not that they are equal.", "",
          "## Reproduce", "", "`.venv/Scripts/python.exe scripts/v4/eval_octane_v4.py`", ""]
    (OUT / "README.md").write_text("\n".join(L), encoding="utf-8")


if __name__ == "__main__":
    main()
