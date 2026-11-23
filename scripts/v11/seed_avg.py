"""v11c: does averaging 5 LightGBM seeds beat a single seed for the pair models? Same outer GroupKFold(5) dev folds as v5.
Rule (docs/v11c_predeclaration.md): adopt only if the paired system-bootstrap 95 % interval of the gain is above 0 and the mean gain is at least 1 % of the single-seed MAE."""
import sys, json
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
import nested_tune as N

out = {}
for which in ("visc", "dens"):
    S = N.Setup(which)
    n = len(S.R)
    c1, c5 = np.zeros(n), np.zeros(n)
    for a, b in GroupKFold(n_splits=5).split(S.dev_rows, groups=S.sysarr[S.dev_rows]):
        tr, te = S.dev_rows[a], S.dev_rows[b]
        trs = set(S.sysarr[tr])
        cs = [S.fit_predict({**S.default, "random_state": 42 + s, "n_jobs": 4}, 2, trs, te)[te] for s in range(5)]
        c1[te] = cs[0]; c5[te] = np.mean(cs, axis=0)
    a1, a5 = S.macro(S.err_from(c1), S.dev_rows), S.macro(S.err_from(c5), S.dev_rows)
    d = (a1 - a5.reindex(a1.index)).to_numpy()
    rng = np.random.default_rng(42)
    bs = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)])
    r = dict(n_systems=len(d), mae_single=float(a1.mean()), mae_avg5=float(a5.mean()), gain=float(d.mean()), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)))
    r["adopt"] = bool(r["lo"] > 0 and r["gain"] >= 0.01 * r["mae_single"])
    out[which] = r
    print(which, r, flush=True)
json.dump(out, open(N.ROOT / "results/v11/seed_avg.json", "w"), indent=1)
