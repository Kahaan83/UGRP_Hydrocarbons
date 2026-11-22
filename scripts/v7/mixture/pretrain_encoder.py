"""v7 Task 11 step 0: refit the Track B pure-compound encoder (hidden 300, depth 3, dropout 0.2, the v6 Track C config) with the long
cap (epochs 1000, patience 100) and curve logging. Same recipe as scripts/v6/mixture/pretrain_encoder.py, which cannot be imported
(it runs on import); YSI labels of the ysi_holdout lock are masked. Saves results/v7/mixture/encoder.pt.
Usage: FUELBLEND_DEVICE=cuda python scripts/v7/mixture/pretrain_encoder.py [--quick]"""
import sys, time, os
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, torch
from fuelblend.ml import dmpnn as D
from fuelblend.ml.curves import CurveLogger
from fuelblend.v7 import runner

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/mixture"
OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs"))
CFG = dict(hidden=300, depth=3, dropout=0.2)
comp, folds, tasks = runner.load_track_b_pool()
g = comp.grp.values
Yraw = comp[tasks].values.astype(float)
Mk = ~np.isnan(Yraw)
rng = np.random.default_rng(0)
ug = np.unique(g)
vg = set(rng.choice(ug, int(.15 * len(ug)), replace=False))
va = np.where(np.isin(g, list(vg)))[0]
tr = np.where(~np.isin(g, list(vg)))[0]
Ym = np.where(Mk[tr], Yraw[tr], np.nan)
mu = np.nan_to_num(np.nanmean(Ym, 0))
sd = np.nanstd(Ym, 0) + 1e-9
sd = np.where(np.isnan(sd) | (Mk[tr].sum(0) < 2) | (sd < 1e-6), 1.0, sd)
Z = np.nan_to_num((Yraw - mu) / sd)
feats = [D.featurize(s) for s in comp.smiles]
lg = CurveLogger(OUT / ("curves_encoder_quick.csv" if QUICK else "curves_encoder.csv"), tb_dir=RUNS / "mixture" / "encoder")
t0 = time.time()
ms = D.fit_ensemble([feats[i] for i in tr], Z[tr], Mk[tr], [feats[i] for i in va], Z[va], Mk[va], n_models=2 if QUICK else 5,
                    epochs=6 if QUICK else 1000, patience=3 if QUICK else 100, seed=7, logger=lg.child(model="encoder"), **CFG)
lg.close()
torch.save({"cfg": CFG, "tasks": tasks, "mu": mu, "sd": sd, "states": [m.state_dict() for m in ms]},
           OUT / ("encoder_quick.pt" if QUICK else "encoder.pt"))
print("done", f"{time.time() - t0:.0f}s", len(comp), "compounds", "best epochs", [getattr(m, "best_epoch_", None) for m in ms], flush=True)
