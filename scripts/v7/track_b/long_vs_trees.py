"""Post-hoc: long-cap D-MPNN (curves_v1_long/oof.csv) vs better of GBM/RF (descriptor_gbm/oof.csv). Exploratory."""
import sys; sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
R = Path(__file__).resolve().parents[3]; B = R / "results/v7/track_b"
a = pd.read_csv(B / "curves_v1_long/oof.csv").set_index("inchikey"); t = pd.read_csv(B / "descriptor_gbm/oof.csv").set_index("inchikey").reindex(a.index)
tasks = [c[:-4] for c in a.columns if c.endswith("_obs")]
rows, ps = [], []
for k in tasks:
    y = a[k + "_obs"].values; gnn = a[k + "_long"].values; m = ~np.isnan(y) & ~np.isnan(gnn) & ~np.isnan(t[k + "_gbm"].values)
    en = np.abs(gnn[m] - y[m]); eg = np.abs(t[k + "_gbm"].values[m] - y[m]); er = np.abs(t[k + "_rf"].values[m] - y[m])
    best = eg if eg.mean() <= er.mean() else er; d = best - en; gg = a.grp.values[m]; ug = np.unique(gg)
    S = np.array([d[gg == u].sum() for u in ug]); N = np.array([(gg == u).sum() for u in ug], float)
    idx = np.random.default_rng(0).integers(0, len(ug), (4000, len(ug))); est = S[idx].sum(1) / N[idx].sum(1)
    p = 2 * min((est <= 0).mean(), (est >= 0).mean()); ps.append(p)
    rows.append(dict(target=k, n=int(m.sum()), mae_gnn_long=en.mean(), mae_gbm=eg.mean(), mae_rf=er.mean(), gain=d.mean(),
                     lo=np.percentile(est, 5), hi=np.percentile(est, 95), p_boot=p))
res = pd.DataFrame(rows); o = np.argsort(ps); adj = np.empty(len(ps)); run = 0
for r_, i in enumerate(o): run = max(run, min(1, (len(ps) - r_) * ps[i])); adj[i] = run
res["p_holm"] = adj; res["gnn_earns_place"] = (res.p_holm < 0.10) & (res.gain > 0)
res.to_csv(B / "curves_v1_long/long_vs_trees.csv", index=False); print(res.round(3).to_string())
