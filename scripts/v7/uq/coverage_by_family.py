"""Task 19: coverage by family under global, Mondrian and shift-weighted one-sided margins (dev rows only, grouped leave-one-group-out calibration).
Unsafe side: flash = over-prediction (pred - measured); freezing and YSI = under-prediction (measured - pred).
For each held-out group the margin is calibrated on all other groups; coverage = share with unsafe-side score <= margin. Outputs results/v7/uq/."""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.uq import conformal_v7 as C  # noqa: E402
from fuelblend.v7 import guard  # noqa: E402

OUT = ROOT / "results/v7/uq"
OUT.mkdir(parents=True, exist_ok=True)
ALPHA = 0.10


def tasks():
    t = {}
    f = pd.read_csv(ROOT / "results/v5/BASE/rows_t2_flash.csv")
    f = f[~f.already_read_locked.astype(bool)].reset_index(drop=True)
    t["flash_designated"] = dict(score=(f.pred_designated - f.measured).to_numpy(), fam=f.method_family.to_numpy(), grp=f.group.to_numpy(), feat=f[["pred_designated"]].to_numpy())
    t["flash_lfl_overlay"] = dict(score=(f.pred_lfl_overlay - f.measured).to_numpy(), fam=f.method_family.to_numpy(), grp=f.group.to_numpy(), feat=f[["pred_lfl_overlay"]].to_numpy())
    g = pd.read_csv(ROOT / "results/v7/freezing/gp/predictions.csv")
    fam = pd.read_csv(ROOT / "results/v7/freezing/families.csv").set_index("group").family
    for m in ("logo_constant", "logo_ridge"):
        t[f"freezing_{m}"] = dict(score=(g.freezing_C - g[m]).to_numpy(), fam=g.group.map(fam).to_numpy(), grp=g.group.to_numpy(), feat=g[[m]].to_numpy())
    from rdkit import Chem
    y = pd.read_csv(ROOT / "results/v5/BASE/delta_ml_ysi/oof_predictions.csv")
    ik = [Chem.MolToInchiKey(Chem.MolFromSmiles(s)) for s in y.smiles]
    locked = guard.locked_ids("ysi_holdout")
    keep = np.array([guard.canon(k) not in locked for k in ik])
    y = y[keep].reset_index(drop=True)
    for m in ("pred_ridge", "pred_lgbm"):
        t[f"ysi_{m}"] = dict(score=(y.ysi - y[m]).to_numpy(), fam=y.cls.to_numpy(), grp=y.group.to_numpy(), feat=y[[m, "nC"]].to_numpy())
    return t


def main():
    rows = []
    for name, d in tasks().items():
        s, fam, grp, X = d["score"], d["fam"], d["grp"], d["feat"]
        cov = {k: np.full(len(s), np.nan) for k in ("global", "mondrian", "weighted")}
        for g in np.unique(grp):
            te, tr = grp == g, grp != g
            mg = C.one_sided_margin(s[tr], ALPHA)
            mm = C.group_conditional_margins(s[tr], fam[tr], ALPHA)
            w = C.shift_weights(X[tr], X[te]) if te.sum() >= 1 else None
            mw = C.one_sided_margin(s[tr], ALPHA, w)
            cov["global"][te] = s[te] <= mg
            cov["mondrian"][te] = [s[i] <= mm.get(str(fam[i]), mm["__global__"]) for i in np.where(te)[0]]
            cov["weighted"][te] = s[te] <= mw
        for f in ["__all__"] + sorted(set(fam)):
            mk = np.ones(len(s), bool) if f == "__all__" else fam == f
            r = dict(task=name, family=f, n=int(mk.sum()), n_groups=int(len(np.unique(grp[mk]))))
            for k, v in cov.items():
                r[f"coverage_{k}"] = float(np.mean(v[mk]))
            r["margin_global_full"] = C.one_sided_margin(s, ALPHA)
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "coverage_by_family.csv", index=False, float_format="%.4g")
    print(df.round(3).to_string())


if __name__ == "__main__":
    main()
