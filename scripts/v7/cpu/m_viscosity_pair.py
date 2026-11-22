"""Task 13 adapter: viscosity_pair (v5 delta_ml_viscosity). Reproduction gate = dev out-of-fold macro MAE (R-MLV-1, warm dev rows) of the v5 pair LightGBM."""
import sys, json, importlib.util
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[3]


def _redirect_removed_v3():
    """scripts/v3 was removed from the repo (commit b85a853); the v4 leaderboard script still loads one helper from it at import.
    Point that single load at a copy recovered from git (scratch/legacy), without touching the repo tree."""
    import importlib.util as iu
    orig = iu.spec_from_file_location

    def patched(name, location=None, *a, **k):
        if location is not None and str(location).replace("\\", "/").endswith("scripts/v3/build_viscosity_leaderboard.py"):
            location = ROOT / "scratch/legacy/build_viscosity_leaderboard.py"
        return orig(name, location, *a, **k)
    iu.spec_from_file_location = patched


def load_v5():
    _redirect_removed_v3()
    spec = importlib.util.spec_from_file_location("v5_visc", ROOT / "scripts/v5/base/delta_ml_viscosity.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["v5_visc"] = mod
    spec.loader.exec_module(mod)
    return mod


def run(nt, to, quick, out, tb, T):
    V = load_v5()
    from fuelblend.ml import viscosity_pair as vp
    ctx = vp.context()
    R = V.regenerate_rows(ctx)
    P = V.build_pairs(R, ctx)
    n = len(R)
    sysarr = R.system_id.to_numpy()
    resid_ref = R.ln_pred.to_numpy() - R.ln_meas.to_numpy()
    dev_rows = np.where(R.role.to_numpy() == "dev")[0]
    row_fold = np.full(n, -1)
    for f, (_, te) in enumerate(GroupKFold(n_splits=5).split(dev_rows, groups=sysarr[dev_rows])):
        row_fold[dev_rows[te]] = f
    # prediction units: all pairs of dev rows
    Pd = P[np.isin(P.row.to_numpy(), dev_rows)].reset_index(drop=True)
    prow = Pd.row.to_numpy()
    folds = row_fold[prow]
    grp = sysarr[prow]
    elig = (R.n_comp.to_numpy()[prow] == 2) & (Pd.xx.to_numpy() >= V.XMIN)
    r = (R.ln_meas.to_numpy() - R.ln_pred.to_numpy())[prow]
    y = np.clip(r, -V.WINS, V.WINS) / Pd.xx.to_numpy()
    w = Pd.xx.to_numpy() ** 2
    tot = pd.Series(np.where(elig, w, 0.0)).groupby(grp).transform("sum").to_numpy()
    w = np.where(elig, w / np.where(tot > 0, tot, 1.0), 0.0)
    # step 1: refit with the original parameters, exactly as v5 (per-fold train_set / targets / fit_lgb)
    oof_orig = np.zeros(n)
    for f in range(5):
        te_rows = np.where(row_fold == f)[0]
        tr_sys = set(sysarr[dev_rows[row_fold[dev_rows] != f]])
        Ptr = V.train_set(R, P, tr_sys)
        g, wt = V.targets(R, Ptr)
        model = V.fit_lgb(Ptr, g, wt)
        Pte = P[P.row.isin(te_rows)]
        oof_orig[te_rows] = V.corr_by_row(Pte, model.booster_.predict(Pte[vp.FEATURES]), n)[te_rows]
    dev = R.role.to_numpy() == "dev"
    warm = R.T_K.to_numpy() >= V.T_WARM
    E_orig = resid_ref + oof_orig
    mae_orig = float(V.macro(R, E_orig, dev & warm).mean())
    base = json.load(open(ROOT / "results/v5/BASE/delta_ml_viscosity/readings_dev.json"))["R-MLV-1"]["ml"]["mae_b"]
    gate = dict(model="viscosity_pair", repro_mae=mae_orig, base_mae=base, rel_diff=abs(mae_orig - base) / base, passed=bool(abs(mae_orig - base) / base < 1e-6))
    if not gate["passed"]:
        return gate, None
    # step 2: nested search, then OOF correction per row
    from _trees import nested_lgb, boot_gain
    X = Pd[vp.FEATURES]
    oof_g, res = nested_lgb(X, y, w, elig, grp, folds, cat=["class_pair"], prefix="v7_visc_pair" + ("_quick" if quick else ""), out_dir=out, nt=nt, to=to, tb=tb, quick=quick)
    oof_new = np.bincount(prow, weights=Pd.xx.to_numpy() * oof_g, minlength=n)
    E_new = resid_ref + oof_new
    m = dev & warm
    e_new, e_orig, e_ref = (V.macro(R, E, m) for E in (E_new, E_orig, resid_ref))
    s_ix = e_new.index
    gain, lo, hi = boot_gain(e_orig.reindex(s_ix).to_numpy(), e_new.to_numpy(), np.arange(len(s_ix)))
    summ = dict(model="viscosity_pair", metric="dev OOF macro MAE (ln, warm dev rows)", n_groups=int(len(s_ix)), mae_v7=float(e_new.mean()), mae_refit_original=mae_orig,
                mae_base_file=base, mae_physics_baseline=float(e_ref.mean()), gain_vs_original=gain, gain_lo=lo, gain_hi=hi, best_params=res.best_params.tolist())
    return gate, summ
