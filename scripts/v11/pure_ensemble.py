"""v11c: does averaging the deployed pure-compound models (RF, LightGBM, D-MPNN) beat the deployed single model? Uses saved scaffold-CV OOF predictions.
Predeclared rule (docs/v11c_predeclaration.md): adopt an ensemble for a property only if the paired scaffold-group bootstrap interval of the MAE gain
(deployed minus ensemble) is above 0 and the Holm-adjusted one-sided p over the 10 properties is below 0.05.  Nothing is refit."""
import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]
a = pd.read_csv(ROOT / "results/v7/track_b/curves_v1_long/oof.csv")
b = pd.read_csv(ROOT / "results/v7/track_b/descriptor_gbm/oof.csv")
assert (a.inchikey == b.inchikey).all()
DEPLOYED = {"flash_C": "rf", "cetane": "rf", "dcn": "long", "density": "long", "dhfus_kJ_mol": "rf", "logvisc_mPas_293": "rf", "mon": "long", "ron": "long", "tm_K": "rf", "ysi": "rf"}
rng = np.random.default_rng(42)
rows = []
for t, dep in DEPLOYED.items():
    obs = b[f"{t}_obs"]
    P = {"rf": b[f"{t}_rf"], "gbm": b[f"{t}_gbm"], "long": a[f"{t}_long"]}
    ok = obs.notna() & pd.concat(P, axis=1).notna().all(axis=1)
    g = b.grp[ok].to_numpy()
    cands = {"rf+long": (P["rf"] + P["long"]) / 2, "rf+gbm+long": (P["rf"] + P["gbm"] + P["long"]) / 3}
    e_dep = (P[dep] - obs).abs()[ok].to_numpy()
    ug = np.unique(g)
    idx = [np.where(g == u)[0] for u in ug]
    for name, p in cands.items():
        e = (p - obs).abs()[ok].to_numpy()
        d = e_dep - e
        sd = np.array([d[i].sum() for i in idx]); nn = np.array([len(i) for i in idx])   # row-weighted statistic, resampling scaffold groups
        bs = np.array([(lambda k: sd[k].sum() / nn[k].sum())(rng.integers(0, len(sd), len(sd))) for _ in range(2000)])
        rows.append(dict(property=t, n=int(ok.sum()), deployed=dep, ensemble=name, mae_deployed=float(e_dep.mean()), mae_ensemble=float(e.mean()),
                         gain=float(d.mean()), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), p=float((np.sum(bs <= 0) + 1) / 2001)))
R = pd.DataFrame(rows)
best = R.loc[R.groupby("property").mae_ensemble.idxmin()].copy()   # one candidate per property (lowest ensemble MAE), selection disclosed
order = best.p.rank(method="first").astype(int)
best["p_holm"] = [min(1.0, p * (len(best) - r + 1)) for p, r in zip(best.p, order)]
best["adopt"] = (best.lo > 0) & (best.p_holm < 0.05)
R.to_csv(ROOT / "results/v11/ensemble/all_candidates.csv", index=False)
best.to_csv(ROOT / "results/v11/ensemble/best_per_property.csv", index=False)
print(best[["property", "deployed", "ensemble", "mae_deployed", "mae_ensemble", "gain", "lo", "hi", "p_holm", "adopt"]].round(4).to_string())
