"""v7 final pipeline: fit the deployed pure-compound models on all dev rows (locked ysi_holdout masked).
RF per property (fixed seed) and one multi-task D-MPNN ensemble (v1 config 0, long cap). Writes data/cache/pure_models/ (gitignored).
Usage: FUELBLEND_DEVICE=cuda python scripts/v7/final/fit_pure_models.py [--rf-only]"""
import sys, json
sys.path.insert(0, "src")
import numpy as np, torch
from fuelblend.api import best_methods as BM, predict_pure as PP
from fuelblend.ml import dmpnn as D
from fuelblend.v7 import runner

comp, folds, tasks = runner.load_track_b_pool()
PP.CACHE.mkdir(parents=True, exist_ok=True)
for p, e in BM.PURE.items():
    if e.method == "rf_rdkit2d":
        print("rf", p, PP.fit_rf(p, comp), flush=True)
if "--rf-only" in sys.argv:
    sys.exit()
Y = comp[tasks].values.astype(float); M = ~np.isnan(Y); g = comp.grp.values
ug = np.unique(g); vg = set(np.random.default_rng(0).choice(ug, max(1, int(.15 * len(ug))), replace=False))
va = np.array([i for i in range(len(comp)) if g[i] in vg]); tr = np.array([i for i in range(len(comp)) if g[i] not in vg])
mu = np.nan_to_num(np.nanmean(np.where(M[tr], Y[tr], np.nan), 0)); sd = np.nanstd(np.where(M[tr], Y[tr], np.nan), 0) + 1e-9
sd = np.where(np.isnan(sd), 1.0, sd); Z = np.nan_to_num((Y - mu) / sd)
arch = dict(hidden=200, depth=4, dropout=0.1)
feats = [D.featurize(s) for s in comp.smiles]
ms = D.fit_ensemble([feats[i] for i in tr], Z[tr], M[tr], [feats[i] for i in va], Z[va], M[va], n_models=5, seed=0,
                    epochs=1000, patience=100, **arch)
torch.save([m.state_dict() for m in ms], PP.CACHE / "dmpnn_long.pt")
(PP.CACHE / "dmpnn_long.json").write_text(json.dumps(dict(tasks=tasks, mu=mu.tolist(), sd=sd.tolist(), arch=arch,
    n_train=int(len(tr)), n_val=int(len(va)), best_epochs=[m.best_epoch_ for m in ms], n_locked_ysi=comp.attrs.get("n_locked_ysi"))))
print("dmpnn saved", [m.best_epoch_ for m in ms])
