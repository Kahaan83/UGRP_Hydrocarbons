"""R-MLV-5 (descriptive): ML pair correction on the T2 real-fuel viscosity rows already scored in BASE (reproduces their Refutas values)."""
import sys, warnings, json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts/v5/base"))
warnings.simplefilter("ignore")
import cells_t2_bulk as c
from fuelblend.ml import viscosity_pair as vp
from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins, predict_properties
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
from fuelblend.data.reconstruction_r3 import SubReps
from fuelblend.baselines.flash_unifac import _lookup

reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
ctx = RouteCtx(reg, prov)
t = pd.read_csv(c.OUT / "t2_targets.csv"); t = t[t.property == "viscosity_kin_cSt"]
fu = pd.read_csv(c.OUT / "t2_fuels.csv").set_index("entity")
V5 = c.V5
comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"), "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
old = pd.read_csv(c.OUT / "rows_t2_bulk.csv")
old = old[(old.property == "viscosity_kin_cSt") & (old.model == "refutas_raw")]
rows, cache = [], {}
for ent, g in t.groupby("entity"):
    key = fu.loc[ent, "source"]; src = c.SRC[key]
    fid = int(g.fuel_id.iloc[0]) if key == "ajftd" else g.fuel_id.iloc[0]
    cc = comps[key]
    mf = map_bins(kind_bins(key, cc[cc.fuel_id == fid]), reps)
    keep = mf.w > 0
    idx, w = mf.idx[keep], mf.w[keep] / mf.w[keep].sum()
    smi = [reg.smiles[i] for i in idx]
    x = w / np.array([ctx.mw(s) for s in smi]); x = x / x.sum()
    for r in g.itertuples():
        o = predict_properties(ctx, mf, pd.DataFrame([dict(property=r.property, T_K=r.T_K)])).iloc[0]
        if o.status != "ok" or not np.isfinite(o.pred):
            continue
        try:
            corr = vp.ln_correction(smi, x, float(r.T_K)) if len(smi) > 1 else 0.0
        except Exception:
            continue
        rows.append(dict(entity=ent, source=src, group=fu.loc[ent, "group"], chosen=bool(fu.loc[ent, "chosen_for_group"]), T_K=c.snap(r.T_K),
                         measured=r.measured, ref=o.pred, ln_corr=corr, ml=o.pred * np.exp(corr)))
d = pd.DataFrame(rows)
# reproduction check against committed BASE refutas rows
m = d.merge(old[["entity", "T_K", "pred", "measured"]], on=["entity", "T_K", "measured"], how="left")
rep = float(np.nanmax(np.abs(m.ref - m.pred) / m.pred)) if m.pred.notna().any() else None
d = d[d.chosen]                       # one entity per double-count group, as BASE
d["e_ref"] = np.abs(np.log(d.ref / d.measured)); d["e_ml"] = np.abs(np.log(d.ml / d.measured))
d["s_ref"] = np.log(d.ref / d.measured); d["s_ml"] = np.log(d.ml / d.measured)
out = []
for (s, T), h in d.groupby(["source", "T_K"]):
    out.append(dict(source=s, T_K=T, n_groups=h.group.nunique(), mae_ln_refutas=h.e_ref.mean(), mae_ln_ml=h.e_ml.mean(), bias_refutas=h.s_ref.mean(), bias_ml=h.s_ml.mean(), mean_ln_corr=h.ln_corr.mean()))
o = pd.DataFrame(out)
o.to_csv(ROOT / "results/v5/BASE/delta_ml_viscosity/realfuel_transfer_R-MLV-5.csv", index=False, float_format="%.4f")
allT = pd.DataFrame([dict(scope="all", n_rows=len(d), mae_ln_refutas=d.e_ref.mean(), mae_ln_ml=d.e_ml.mean())])
json.dump(dict(refutas_reproduction_max_rel_diff=rep, n_rows=len(d), mae_ln_refutas=float(d.e_ref.mean()), mae_ln_ml=float(d.e_ml.mean()),
               temps_ml_not_worse=int((o.mae_ln_ml <= o.mae_ln_refutas).sum()), temps_total=int(len(o))),
          open(ROOT / "results/v5/BASE/delta_ml_viscosity/readings_R-MLV-5.json", "w"), indent=1)
print(o.round(3).to_string()); print(rep, allT)
