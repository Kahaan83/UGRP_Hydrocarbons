"""Debugging task: in-house D-MPNN schedule/budget variants on flash_C vs the Task 7 Chemprop OOF. See PREDECLARATION.md."""
import sys, json, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd
from fuelblend.ml import dmpnn as D
from fuelblend.ml.curves import CurveLogger
from fuelblend.v7 import runner

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/track_b" / ("debug_flash_chemprop_quick" if QUICK else "debug_flash_chemprop")
OUT.mkdir(parents=True, exist_ok=True)
NM = 2 if QUICK else 5
VARIANTS = {"b_cosine150": dict(scheduler="cosine", epochs=150, patience=25),
            "c_plateau150": dict(scheduler="plateau", epochs=150, patience=25),
            "d_const1000": dict(scheduler="none", epochs=1000, patience=100)}
if QUICK:
    VARIANTS = {k: dict(v, epochs=6, patience=3) for k, v in VARIANTS.items()}
T = "flash_C"
comp, folds, tasks = runner.load_track_b_pool()
g = comp.grp.values
y = comp[T].values
has = ~np.isnan(y)
feats = [D.featurize(s) for s in comp.smiles]
ref = pd.read_csv(R / "results/v7/track_b/chemprop_parity/oof.csv").set_index("inchikey").reindex(comp.inchikey)
cp = ref[T + "_chemprop"].values
pred = {k: np.full(len(comp), np.nan) for k in VARIANTS}
lg = CurveLogger(OUT / "curves.csv")
t0 = time.time()
for k in range(5):
    tr = np.where((folds != k) & has)[0]
    te = np.where((folds == k) & has)[0]
    ug = np.unique(g[tr])
    rng = np.random.default_rng(k)
    vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
    va = np.array([i for i in tr if g[i] in vg])
    tr2 = np.array([i for i in tr if g[i] not in vg])
    mu, sd = y[tr2].mean(), y[tr2].std() + 1e-9
    Z = (np.nan_to_num(y) - mu) / sd
    Mk = has[:, None].astype(float)
    for name, kw in VARIANTS.items():
        ms = D.fit_ensemble([feats[i] for i in tr2], Z[tr2, None], Mk[tr2], [feats[i] for i in va], Z[va, None], Mk[va],
                            n_models=NM, seed=k, hidden=200, depth=4, dropout=0.1, logger=lg.child(model=name, fold=k), **kw)
        pred[name][te] = D.predict(ms, [feats[i] for i in te])[0][:, 0] * sd + mu
        print(name, "fold", k, f"{time.time() - t0:.0f}s", flush=True)
lg.close()
rows = []
m0 = has & ~np.isnan(cp)
for name, p in pred.items():
    m = m0 & ~np.isnan(p)
    lo, hi = runner.paired_boot_gain(p[m] - y[m], cp[m] - y[m], g[m])
    rows.append(dict(variant=name, n=int(m.sum()), mae_variant=np.abs(p[m] - y[m]).mean(), mae_chemprop=np.abs(cp[m] - y[m]).mean(),
                     mae_inhouse_task7=np.abs(ref[T + "_inhouse"].values[m] - y[m]).mean(), gain_lo=lo, gain_hi=hi,
                     closes_gap=bool(hi >= 0)))
res = pd.DataFrame(rows)
oo = comp[["inchikey", "smiles", "grp"]].copy()
oo[T + "_obs"] = y
oo[T + "_chemprop"] = cp
for name, p in pred.items():
    oo[T + "_" + name] = p
runner.write_results(OUT, summary=res, oof=oo, meta=dict(task="debug_flash_vs_chemprop", n_models=NM, variants=VARIANTS,
                                                       total_s=time.time() - t0, n_locked_ysi=comp.attrs.get("n_locked_ysi")))
print(res.round(3).to_string())
