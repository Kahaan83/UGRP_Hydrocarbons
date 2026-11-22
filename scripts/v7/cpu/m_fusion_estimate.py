"""Task 13 adapter: fusion_estimate Tm model (v5 delta_ws3_fusion). Gate = leave-one-family-out Tm MAE, all hydrocarbons and C8 to C20
(readings_fusion.json R-W3-5 and R-W3-3). Only the Tm model is tuned; the dS_fus model is not part of this retrain."""
import sys, json, importlib.util
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]


def load_v5():
    spec = importlib.util.spec_from_file_location("v5_fusion", ROOT / "scripts/v5/base/delta_ws3_fusion.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["v5_fusion"] = mod
    spec.loader.exec_module(mod)
    return mod


def run(nt, to, quick, out, tb, T):
    V = load_v5()
    FE = V.FE
    # the measured table has since gained 10 compounds (119 vs 109 rows; Tm differs by up to 0.26 K), so the rows are pinned to the
    # compounds that produced the BASE number: fusion_oof.csv (read only) holds their inchikey, smiles and Tm
    base_rows = pd.read_csv(ROOT / "results/v5/BASE/delta_ws3_pure_gaps/fusion_oof.csv")
    train = base_rows[["inchikey", "smiles", "name", "tm_K", "tm_source", "dh_kJ", "dh_source"]].reset_index(drop=True)
    X = FE.features(train.smiles)
    fam = train.smiles.map(FE.family_of).to_numpy()
    y = train.tm_K.to_numpy(float)
    nC = X.n_carbon.to_numpy()
    c820 = (nC >= 8) & (nC <= 20)
    p_orig = FE.oof_predict(X, y, fam)
    r = json.load(open(ROOT / "results/v5/BASE/delta_ws3_pure_gaps/readings_fusion.json"))
    base_all, base_c = r["R-W3-5"]["mae_tm_lofo_all_K"], r["R-W3-3"]["mae_tm_lofo_K"]
    mae_all, mae_c = float(np.abs(p_orig - y).mean()), float(np.abs(p_orig - y)[c820].mean())
    ok = abs(mae_all - base_all) / base_all < 1e-6 and abs(mae_c - base_c) / base_c < 1e-6
    gate = dict(model="fusion_estimate", repro_mae=mae_all, base_mae=base_all, repro_mae_c8_c20=mae_c, base_mae_c8_c20=base_c,
                rel_diff=max(abs(mae_all - base_all) / base_all, abs(mae_c - base_c) / base_c), passed=bool(ok))
    if not gate["passed"]:
        return gate, None
    from _trees import nested_lgb, boot_gain
    codes = {f: i for i, f in enumerate(sorted(set(fam)))}
    folds = np.array([codes[f] for f in fam])
    oof, res = nested_lgb(X, y, np.ones(len(y)), np.ones(len(y), bool), fam, folds, cat=None, prefix="v7_fusion_tm" + ("_quick" if quick else ""),
                          out_dir=out, nt=nt, to=to, tb=tb, quick=quick, n_jobs=2)
    e_new, e_orig = np.abs(oof - y), np.abs(p_orig - y)
    gain, lo, hi = boot_gain(e_orig, e_new, fam)
    gc, lc, hc = boot_gain(e_orig[c820], e_new[c820], fam[c820])
    jo = X.joback_tm.to_numpy()
    summ = dict(model="fusion_estimate", metric="leave-one-family-out Tm MAE (K), all hydrocarbons", n_groups=int(len(np.unique(fam))), n=len(y),
                mae_v7=float(e_new.mean()), mae_refit_original=mae_all, mae_base_file=base_all, mae_physics_baseline=float(np.abs(jo - y).mean()),
                gain_vs_original=gain, gain_lo=lo, gain_hi=hi, mae_v7_c8_c20=float(e_new[c820].mean()), mae_refit_original_c8_c20=mae_c,
                gain_c8_c20=gc, gain_c8_c20_lo=lc, gain_c8_c20_hi=hc, best_params=res.best_params.tolist())
    return gate, summ
