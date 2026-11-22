"""Task 13 adapter: ysi_lgbm (the LightGBM comparator in v5 delta_ml_ysi). Gate = OOF MAE of that LightGBM over all 222 v5 compounds
(R-MLY-1 mae.lgbm). The v7 comparison itself uses dev compounds only: 53 of the 222 are in the ysi_holdout lock and are excluded."""
import sys, json, importlib.util
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]


def load_v5():
    spec = importlib.util.spec_from_file_location("v5_ysi", ROOT / "scripts/v5/base/delta_ml_ysi.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["v5_ysi"] = mod
    spec.loader.exec_module(mod)
    return mod


def run(nt, to, quick, out, tb, T):
    from fuelblend.v7 import guard
    V = load_v5()
    YM = V.YM
    df = YM.load_training()
    smi = df.smiles.to_numpy()
    y = df.ysi.to_numpy(float)
    X = np.array([YM.features(s) for s in smi])
    groups = np.array([YM.group_key(s) for s in smi])
    # step 1: gate on all 222 rows, exactly as v5 (reproduces a published number; no tuning, not used for any v7 claim)
    p = np.full(len(y), np.nan)
    for tr, te in GroupKFold(5).split(X, y, groups):
        p[te] = V.lgbm_fit_predict(X[tr], y[tr], X[te])
    mae = float(np.abs(p - y).mean())
    base = json.load(open(ROOT / "results/v5/BASE/delta_ml_ysi/readings.json"))["R-MLY-1"]["mae"]["lgbm"]
    gate = dict(model="ysi_lgbm", repro_mae=mae, base_mae=base, rel_diff=abs(mae - base) / base, passed=bool(abs(mae - base) / base < 1e-6),
                note="gate uses all 222 v5 rows; v7 comparison below uses dev compounds only")
    if not gate["passed"]:
        return gate, None
    # step 2: dev-only comparison
    iks = [Chem.MolToInchiKey(Chem.MolFromSmiles(s)) for s in smi]
    keep = np.array([guard.canon(k) not in guard.locked_ids("ysi_holdout") for k in iks])
    Xd, yd, gd = X[keep], y[keep], groups[keep]
    folds = np.full(len(yd), -1)
    for f, (_, te) in enumerate(GroupKFold(5).split(Xd, yd, gd)):
        folds[te] = f
    p_orig = np.full(len(yd), np.nan)
    for f in range(5):
        tr, te = np.where(folds != f)[0], np.where(folds == f)[0]
        p_orig[te] = V.lgbm_fit_predict(Xd[tr], yd[tr], Xd[te])
    from _trees import nested_lgb, boot_gain
    w = YM.weights(yd)
    Xf = pd.DataFrame(Xd, columns=[f"f{i}" for i in range(Xd.shape[1])])
    oof, res = nested_lgb(Xf, yd, np.asarray(w, float), np.ones(len(yd), bool), gd, folds, cat=None, prefix="v7_ysi_lgbm" + ("_quick" if quick else ""),
                          out_dir=out, nt=nt, to=to, tb=tb, quick=quick, n_jobs=2)
    e_new, e_orig = np.abs(oof - yd), np.abs(p_orig - yd)
    gain, lo, hi = boot_gain(e_orig, e_new, gd)
    const = np.array([yd[folds != folds[i]].mean() for i in range(len(yd))])
    summ = dict(model="ysi_lgbm", metric="dev-only OOF MAE (YSI)", n_groups=int(len(np.unique(gd))), n=int(len(yd)), mae_v7=float(e_new.mean()),
                mae_refit_original=float(e_orig.mean()), mae_base_file_all222=base, mae_physics_baseline=float(np.abs(const - yd).mean()),
                gain_vs_original=gain, gain_lo=lo, gain_hi=hi, best_params=res.best_params.tolist())
    return gate, summ
