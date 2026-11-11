"""G1 freezing recomputed WITHOUT censored ('<') AJFTD values (correction 2026-11-12).  Uses G1's committed predictions (unchanged: the
model prediction does not depend on the target) and replaces `measured` by the median of uncensored rows only.  Writes
results/v5/S3/corrected_2026-11-12_g1_freezing_uncensored.csv.  G1 outputs are not touched."""
from __future__ import annotations
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("g1_eval_fz", ROOT / "scripts/v5/eval_ajftd_gcxgc_v5.py")
g1 = importlib.util.module_from_spec(spec); sys.modules["g1_eval_fz"] = g1; spec.loader.exec_module(g1)


def uncensored_freezing(pred: pd.DataFrame, tg: pd.DataFrame) -> pd.DataFrame:
    fz = pred[(pred.property == "freezing_point_C") & np.isfinite(pred.err)].copy()
    m = tg[tg.property == "freezing_point_C"][["fuel_id", "measured"]].rename(columns={"measured": "measured_unc"})
    fz = fz.merge(m, on="fuel_id", how="left")
    fz["n_censored_only"] = fz.measured_unc.isna()
    fz = fz[fz.measured_unc.notna()].copy()
    fz["err_unc"] = fz.pred - fz.measured_unc
    return fz


def summary(fz: pd.DataFrame) -> dict:
    e = fz.err_unc.to_numpy(float); y = fz.measured_unc.to_numpy(float); n = len(y)
    loo = np.abs(y - (y.sum() - y) / (n - 1))
    return dict(n=n, bias=e.mean(), mae=np.abs(e).mean(), loo_const_mae=loo.mean())


def main():
    pred = pd.read_csv(ROOT / "results/v5/G1/predictions.csv")
    old = pred[(pred.property == "freezing_point_C") & np.isfinite(pred.err)]
    fz = uncensored_freezing(pred, g1.targets(exclude_censored=True))
    sets = {"all eligible": fz, "not shared with another source": fz[~fz.shared_with_other_source.astype(bool)],
            "S4 role: dev": fz[fz.split_role == "dev"], "S4 role: locked": fz[fz.split_role == "locked"]}
    rows = [dict(set="all eligible (G1 as committed, incl. censored)", n=len(old), bias=old.err.mean(), mae=old.err.abs().mean(),
                 loo_const_mae=np.nan)]
    rows += [dict(set=k, **summary(v)) for k, v in sets.items()]
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/v5/S3/corrected_2026-11-12_g1_freezing_uncensored.csv", index=False, float_format="%.4g")
    print(out.to_string())
    print("dropped (censored only):", sorted(pred[(pred.property == "freezing_point_C") & np.isfinite(pred.err)].fuel_id.unique()
                                              .tolist() and set(old.fuel_id) - set(fz.fuel_id)))


if __name__ == "__main__":
    main()
