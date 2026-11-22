"""Task 15: heavy-tail GP and novelty flag for freezing, LOGO (52 dev groups) and leave-one-family-out (4 families). Dev rows only.
  python scripts/v7/freezing/lofo_gp.py [--quick]"""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts/v6/freezing"))
from fuelblend.v7 import freeze_gp as F  # noqa: E402

OUT = ROOT / "results/v7/freezing/gp"
LIVE = 13.8427
QUICK = "--quick" in sys.argv
warnings.filterwarnings("ignore")


def build():
    import common as C
    from fuelblend.api.predict_route import _tmap_freezing
    ds = C.dataset()
    feat = pd.read_csv(ROOT / "results/v6/freezing/features_preds.csv", dtype={"fuel_id": str})
    dev = feat[feat.role == "dev"].reset_index(drop=True)
    ents = C.entities(ds)
    tdes = _tmap_freezing()
    rows = []
    for e in dev.entity:
        for r in ents[e]["mapped"].itertuples():
            t = tdes.get(r.ik) or {}
            rows.append(dict(fuel=e, smiles=r.smiles, frac=r.frac, tm_K=float(t.get("tm_K", np.nan))))
    X = F.heavy_tail_features(pd.DataFrame(rows)).loc[dev.entity]
    fam = pd.read_csv(ROOT / "results/v7/freezing/families.csv").set_index("group").family
    ref = pd.read_csv(ROOT / "results/v6/freezing/dev_oof_predictions.csv")
    assert (ref.entity.to_numpy() == dev.entity.to_numpy()).all()
    return dev, X, fam.loc[dev.group].to_numpy(), ref


def folds_of(scheme, groups, fam):
    key = groups if scheme == "logo" else fam
    return [(np.where(key != k)[0], np.where(key == k)[0]) for k in pd.unique(key)]


def run_scheme(scheme, dev, X, fam, ref, lim):
    from sklearn.linear_model import RidgeCV
    y = dev.freezing_C.to_numpy()
    Xv = X.to_numpy(float)
    groups = dev.group.to_numpy()
    P = {k: np.full(len(y), np.nan) for k in ("constant", "ridge", "gp")}
    flag = np.zeros(len(y), bool); nov = np.full(len(y), np.nan)
    for tr, te in folds_of(scheme, groups, fam):
        med = np.nanmedian(Xv[tr], 0)
        A, B = np.where(np.isnan(Xv[tr]), med, Xv[tr]), np.where(np.isnan(Xv[te]), med, Xv[te])
        mu, sd = F.standardise(A)
        Za, Zb = (A - mu) / sd, (B - mu) / sd
        P["constant"][te] = y[tr].mean()
        inner = [(np.where(groups[tr] != g)[0], np.where(groups[tr] == g)[0]) for g in pd.unique(groups[tr])]
        rg = RidgeCV(alphas=np.logspace(-2, 3, 12), cv=inner).fit(Za, y[tr])
        P["ridge"][te] = rg.predict(Zb)
        gp = F.fit_gp(Za, y[tr], restarts=0 if QUICK else 5)
        P["gp"][te] = gp.predict(Zb)
        thr = F.novelty_threshold(A)
        nov[te] = F.novelty_score(A, B)
        flag[te] = nov[te] > thr
    P["baseline_b_unifac_svl"] = ref.b_unifac_svl.to_numpy()
    P["limiting_species"] = lim
    if scheme == "logo":
        P["v6_lgbm_residual"] = ref.f_lgbm_residual_on_e_isomer_tiles.to_numpy()
    return y, P, flag, nov, groups


def boot(d, groups, n=4000, seed=0):
    ug = pd.unique(groups)
    m = np.array([d[groups == g].mean() for g in ug])
    idx = np.random.default_rng(seed).integers(0, len(ug), (n, len(ug)))
    est = m[idx].mean(1)
    return float(m.mean()), float(np.percentile(est, 5)), float(np.percentile(est, 95))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dev, X, fam, ref = build()
    lim = pd.read_csv(ROOT / "results/v7/freezing/baseline_d/predictions.csv").pred_limiting.to_numpy()
    assert (pd.read_csv(ROOT / "results/v7/freezing/baseline_d/predictions.csv").entity.to_numpy() == dev.entity.to_numpy()).all()
    X.assign(entity=dev.entity.to_numpy(), family=fam).to_csv(OUT / "features.csv", index=False, float_format="%.5g")
    cov, allrows, nv, preds = [], [], [], dev[["entity", "group", "freezing_C"]].assign(family=fam)
    for scheme in ("logo", "lofo"):
        y, P, flag, nov, groups = run_scheme(scheme, dev, X, fam, ref, lim)
        for m, p in P.items():
            preds[f"{scheme}_{m}"] = p
        e = {m: np.abs(p - y) for m, p in P.items()}
        for m, p in P.items():
            row = dict(scheme=scheme, model=m, n=len(y), n_groups=len(np.unique(groups)), mae=float(np.nanmean(e[m])), bias=float(np.nanmean(p - y)),
                       mean_over=float(np.nanmean(np.clip(p - y, 0, None))), mean_under=float(np.nanmean(np.clip(y - p, 0, None))))
            for ref_m in ("constant", "ridge"):
                g, lo, hi = boot(e[ref_m] - e[m], groups)
                row.update({f"gain_vs_{ref_m}": g, f"gain_vs_{ref_m}_lo": lo, f"gain_vs_{ref_m}_hi": hi})
            allrows.append(row)
        preds[f"{scheme}_novelty"], preds[f"{scheme}_flag"] = nov, flag
        for f_ in np.unique(fam):
            mk = fam == f_
            for m in ("gp", "constant", "ridge", "baseline_b_unifac_svl", "limiting_species"):
                cov.append(dict(scheme=scheme, family=f_, model=m, n_groups=int(len(np.unique(groups[mk]))), n=int(mk.sum()),
                                coverage_live_margin=float(((y[mk] - P[m][mk]) <= LIVE).mean()), mae=float(np.nanmean(e[m][mk]))))
        nv.append(dict(scheme=scheme, n_flagged=int(flag.sum()), flagged_share=float(flag.mean()), n_groups_flagged=int(len(np.unique(groups[flag]))),
                       gp_mae_flagged=float(e["gp"][flag].mean()) if flag.any() else np.nan, gp_mae_unflagged=float(e["gp"][~flag].mean()) if (~flag).any() else np.nan,
                       flagged_families=";".join(sorted(set(fam[flag])))))
    pd.DataFrame(nv).to_csv(OUT / ("novelty_quick.csv" if QUICK else "novelty.csv"), index=False, float_format="%.5g")
    pd.DataFrame(allrows).to_csv(OUT / ("summary_quick.csv" if QUICK else "summary.csv"), index=False, float_format="%.5g")
    pd.DataFrame(cov).to_csv(OUT / ("coverage_by_family_quick.csv" if QUICK else "coverage_by_family.csv"), index=False, float_format="%.5g")
    preds.to_csv(OUT / ("predictions_quick.csv" if QUICK else "predictions.csv"), index=False, float_format="%.5g")
    print(pd.DataFrame(allrows).round(3).to_string())


if __name__ == "__main__":
    main()
