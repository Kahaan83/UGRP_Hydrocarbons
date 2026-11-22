"""Task 13 adapter: density_pair (v5 delta_ml_density). Reproduction gate = dev out-of-fold macro MAE (kg/m3) of the v5 pair LightGBM (R-MLD-1 ml)."""
import sys, json, importlib.util
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[3]


def load_v5():
    spec = importlib.util.spec_from_file_location("v5_dens", ROOT / "scripts/v5/base/delta_ml_density.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["v5_dens"] = mod
    spec.loader.exec_module(mod)
    return mod


def run(nt, to, quick, out, tb, T):
    V = load_v5()
    from fuelblend.ml import density_pair as dp
    R, desc = V.regenerate()
    P = V.build_pairs(R, desc)
    n = len(R)
    sysarr = R.system_id.to_numpy()
    dev = R.role.to_numpy() == "dev"
    dev_rows = np.where(dev)[0]
    row_fold = np.full(n, -1)
    for f, (_, te) in enumerate(GroupKFold(n_splits=5).split(dev_rows, groups=sysarr[dev_rows])):
        row_fold[dev_rows[te]] = f
    err_ideal = R.pred.to_numpy() - R.meas.to_numpy()

    def err_from(corr):
        return R.mm.to_numpy() / (R.v_ideal.to_numpy() + corr) * 1000.0 - R.meas.to_numpy()

    oof_orig = np.zeros(n)
    for f in range(5):
        te_rows = np.where(row_fold == f)[0]
        tr_sys = set(sysarr[dev_rows[row_fold[dev_rows] != f]])
        Ptr = V.train_set(R, P, tr_sys)
        g, wt = V.targets(R, Ptr)
        model = V.fit_lgb(Ptr, g, wt)
        Pte = P[P.row.isin(te_rows)]
        oof_orig[te_rows] = V.corr_by_row(Pte, model.booster_.predict(Pte[dp.FEATURES]), n)[te_rows]
    E_orig = err_from(oof_orig)
    mae_orig = float(V.macro(R, E_orig, dev).mean())
    base = json.load(open(ROOT / "results/v5/BASE/delta_ml_density/readings_dev.json"))["R-MLD-1"]["ml"]["mae_b"]
    gate = dict(model="density_pair", repro_mae=mae_orig, base_mae=base, rel_diff=abs(mae_orig - base) / base, passed=bool(abs(mae_orig - base) / base < 1e-6))
    if not gate["passed"]:
        return gate, None
    from _trees import nested_lgb, boot_gain
    Pd = P[np.isin(P.row.to_numpy(), dev_rows)].reset_index(drop=True)
    prow = Pd.row.to_numpy()
    folds, grp = row_fold[prow], sysarr[prow]
    elig = (R.n_comp.to_numpy()[prow] == 2) & (Pd.xx.to_numpy() >= V.XMIN)
    ve = (R.v_meas.to_numpy() - R.v_ideal.to_numpy())[prow]
    y = np.clip(ve, -V.WINS, V.WINS) / Pd.xx.to_numpy()
    w = Pd.xx.to_numpy() ** 2
    tot = pd.Series(np.where(elig, w, 0.0)).groupby(grp).transform("sum").to_numpy()
    w = np.where(elig, w / np.where(tot > 0, tot, 1.0), 0.0)
    oof_g, res = nested_lgb(Pd[dp.FEATURES], y, w, elig, grp, folds, cat=["class_pair"], prefix="v7_dens_pair" + ("_quick" if quick else ""),
                            out_dir=out, nt=nt, to=to, tb=tb, quick=quick)
    oof_new = np.bincount(prow, weights=Pd.xx.to_numpy() * oof_g, minlength=n)
    E_new = err_from(oof_new)
    e_new, e_orig, e_id = (V.macro(R, E, dev) for E in (E_new, E_orig, err_ideal))
    gain, lo, hi = boot_gain(e_orig.reindex(e_new.index).to_numpy(), e_new.to_numpy(), np.arange(len(e_new)))
    summ = dict(model="density_pair", metric="dev OOF macro MAE (kg/m3)", n_groups=int(len(e_new)), mae_v7=float(e_new.mean()), mae_refit_original=mae_orig,
                mae_base_file=base, mae_physics_baseline=float(e_id.mean()), gain_vs_original=gain, gain_lo=lo, gain_hi=hi, best_params=res.best_params.tolist())
    return gate, summ
