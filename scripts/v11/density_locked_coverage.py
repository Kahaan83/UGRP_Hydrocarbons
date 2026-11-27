"""Claim-risk addendum: independent locked coverage of the tuned density ML margin.
Margin from dev OOF only (GroupKFold on system_id); locked systems scored by the model fit on all dev systems.
Does not write to results/v5/BASE."""
import sys, json, warnings, importlib.util
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); warnings.simplefilter("ignore")
spec = importlib.util.spec_from_file_location("dmd", ROOT / "scripts/v5/base/delta_ml_density.py")
D = importlib.util.module_from_spec(spec); sys.modules["dmd"] = D; spec.loader.exec_module(D)
from fuelblend.uq import conformal as cf
R, desc = D.regenerate(); P = D.build_pairs(R, desc)
n = len(R); sysarr = R.system_id.to_numpy(); role = R.role.to_numpy()
dev, lock = role == "dev", role == "locked"; dev_rows = np.where(dev)[0]
def err_from(c): return R.mm.to_numpy() / (R.v_ideal.to_numpy() + c) * 1000.0 - R.meas.to_numpy()
def fitpred(train_sys, test_rows):
    Ptr = D.train_set(R, P, train_sys); g, w = D.targets(R, Ptr); m = D.fit_lgb(Ptr, g, w)
    Pte = P[P.row.isin(test_rows)]
    return D.corr_by_row(Pte, D.predict_all(m, D.fit_baselines(Ptr, g, w), Pte)["ml"] if False else None, n) if False else (m, Pte)
oof = np.zeros(n)
for tr, te in GroupKFold(n_splits=5).split(dev_rows, groups=sysarr[dev_rows]):
    m, Pte = fitpred(set(sysarr[dev_rows[tr]]), dev_rows[te])
    gv = D.predict_all(m, D.fit_baselines(*(lambda Pt: (Pt, *D.targets(R, Pt)))(D.train_set(R, P, set(sysarr[dev_rows[tr]])))), Pte)["ml"]
    oof[dev_rows[te]] = D.corr_by_row(Pte, gv, n)[dev_rows[te]]
m, Pl = fitpred(set(sysarr[dev]), np.where(lock)[0])
Ptr = D.train_set(R, P, set(sysarr[dev])); g, w = D.targets(R, Ptr)
gv = D.predict_all(m, D.fit_baselines(Ptr, g, w), Pl)["ml"]
lockc = D.corr_by_row(Pl, gv, n)
Ed, El = err_from(oof), err_from(lockc)
mg = cf.two_sided_margin(Ed[dev], sysarr[dev], 0.10)
cov = cf.coverage(El[lock], sysarr[lock], mg, two_sided=True)
res = dict(dev_oof_margin=float(mg), locked_coverage=list(map(int, cov)), locked_mae=float(np.abs(El[lock]).mean()))
Path(ROOT / "results/v11/claim_risk/density_locked_coverage.json").write_text(json.dumps(res, indent=1)); print(res)
