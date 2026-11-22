"""Task 13 adapter: freezing_step3 monotone LightGBM residual on the isomer-tile base (v6 Track A step 3), leave-one-group-out on the 52 dev rows.
Gate = the v6 `dev_oof_predictions.csv` column and `dev_logo_summary.csv` MAE; those files store 4 decimals, so the gate tolerance is 1e-4 absolute
per row and 5e-5 on the MAE (the file precision), not 1e-6 relative."""
import sys, json
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/v6/freezing")
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]


def run(nt, to, quick, out, tb, T):
    import step3_models as S
    ds = pd.read_csv(ROOT / "results/v6/freezing/features_preds.csv", dtype={"fuel_id": str})
    dev = ds[ds.role == "dev"].reset_index(drop=True)
    y = dev.freezing_C.to_numpy()
    ref = pd.read_csv(ROOT / "results/v6/freezing/dev_oof_predictions.csv")
    col = "f_lgbm_residual_on_e_isomer_tiles"
    assert (ref.entity.to_numpy() == dev.entity.to_numpy()).all(), "dev row order differs from v6"
    p_orig = S.logo(dev, S.make_lgbm("pred_e"))
    d = float(np.abs(p_orig - ref[col].to_numpy()).max())
    mae = float(np.abs(p_orig - y).mean())
    summ_ref = pd.read_csv(ROOT / "results/v6/freezing/dev_logo_summary.csv")
    base = float(summ_ref[(summ_ref.model == col) & (summ_ref.subset == "all")].mae.iloc[0])
    gate = dict(model="freezing_step3", repro_mae=mae, base_mae=base, max_abs_row_diff=d, tolerance_row=1e-4, tolerance_mae=5e-5,
                rel_diff=abs(mae - base) / base, passed=bool(d <= 1e-4 and abs(mae - base) <= 5e-5))
    if not gate["passed"]:
        return gate, None
    from _trees import nested_lgb, boot_gain
    base_col = dev.pred_e.clip(*S.CLIP).to_numpy()
    X = dev[S.LGB_COLS].reset_index(drop=True)
    res_y = y - base_col
    folds = np.arange(len(dev))              # leave-one-out, as the original
    grp = dev.entity.to_numpy()
    mono = dict(monotone_constraints=S.LGB_MONO)
    oof, res = nested_lgb(X, res_y, np.ones(len(y)), np.ones(len(y), bool), grp, folds, cat=None, prefix="v7_freezing_step3" + ("_quick" if quick else ""),
                          out_dir=out, nt=nt, to=to, tb=tb, quick=quick, n_jobs=1, extra=mono)
    p_new = base_col + oof
    e_new, e_orig = np.abs(p_new - y), np.abs(p_orig - y)
    const = ref["a_constant"].to_numpy()
    gain, lo, hi = boot_gain(e_orig, e_new, grp)
    m64 = (dev.fuel_id != "64").to_numpy()
    summ = dict(model="freezing_step3", metric="dev LOGO MAE (degC), all 52", n_groups=int(len(np.unique(grp))), n=len(y), mae_v7=float(e_new.mean()),
                mae_refit_original=mae, mae_base_file=base, mae_physics_baseline=float(np.abs(const - y).mean()),
                gain_vs_original=gain, gain_lo=lo, gain_hi=hi, mae_v7_no_fuel64=float(e_new[m64].mean()), mae_refit_original_no_fuel64=float(e_orig[m64].mean()),
                bias_v7=float((p_new - y).mean()), over_pred_v7=float(np.clip(p_new - y, 0, None).mean()), over_pred_original=float(np.clip(p_orig - y, 0, None).mean()),
                best_params=res.best_params.tolist())
    return gate, summ
