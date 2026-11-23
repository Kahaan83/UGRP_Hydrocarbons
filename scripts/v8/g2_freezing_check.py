"""G2 freezing-only check on real fuels (AJFTD, AV-33-22, FuelLib set of delta_freeze_t2_unifac) for the designated route."""
import math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np, pandas as pd
from fuelblend.uq.conformal import conformal_rank, group_scores, one_sided_margin

d = pd.read_csv(ROOT / "results/v5/BASE/delta_freeze_t2_unifac/rows.csv")
d = d[d.unifac_desig.notna() & d.freezing_C.notna() & d.C.notna()].copy()
d["e"] = d.unifac_desig - d.freezing_C           # pred - meas; unsafe = negative
d["c"] = d.C - d.freezing_C                      # leave-one-group-out constant
g = d.groupby("group").agg(ae=("e", lambda x: x.abs().mean()), ac=("c", lambda x: x.abs().mean())).reset_index()
gain = (g.ac - g.ae).to_numpy()
rng = np.random.default_rng(42)
bs = np.array([gain[rng.integers(0, len(gain), len(gain))].mean() for _ in range(5000)])
lo, hi = np.quantile(bs, [0.05, 0.95])
a_ok = lo > 0
n = d.group.nunique()
try:
    m = one_sided_margin(d.e.to_numpy(), d.group.to_numpy(), 0.10, "under")
except Exception as ex:
    m = float("nan")
b_ok = math.isfinite(m)
cov = {}
for tr, te in (("AJFTD", "AV-33-22"), ("AV-33-22", "AJFTD")):
    pass
src = d.source.where(d.source.isin(["AJFTD"]), "other")
# leave-source-out: margin from one source, coverage on the other (unsafe side = under-prediction)
cov = {}
for a, b in (("AJFTD", "other"), ("other", "AJFTD")):
    ta, tb = d[src == a], d[src == b]
    if ta.group.nunique() < 2 or tb.empty:
        cov[f"{a}->{b}"] = float("nan"); continue
    ma = one_sided_margin(ta.e.to_numpy(), ta.group.to_numpy(), 0.10, "under")
    gs = group_scores(tb.e.to_numpy(), tb.group.to_numpy(), "under")
    cov[f"{a}->{b}"] = float((gs <= ma).mean()) if math.isfinite(ma) else float("nan")
c_ok = all(np.isfinite(v) and v >= 0.80 for v in cov.values()) and True
out = dict(n_fuels=len(d), n_groups=n, mae_route=float(d.e.abs().mean()), mae_const=float(d.c.abs().mean()),
           bias_route=float(d.e.mean()), gain_group_mean=float(gain.mean()), gain_ci90_lo=float(lo), gain_ci90_hi=float(hi),
           a_beats_constant=bool(a_ok), margin_C=m, b_finite=bool(b_ok), cov=cov, c_transfer=bool(c_ok),
           reading="held" if (a_ok and b_ok and c_ok) else ("no_skill_shown" if not a_ok else ("inf_margin" if not b_ok else "transfer_not_held")))
print(out)
pd.DataFrame([dict(out, cov=str(cov))]).to_csv(ROOT / "results/v8/represented_margins/freezing_g2.csv", index=False)
