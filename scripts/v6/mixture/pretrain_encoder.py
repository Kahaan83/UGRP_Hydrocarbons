"""Track C step 0: fit the Track B D-MPNN once on the full pure pool (config hidden 300, depth 3, dropout 0.2, the config chosen in
3 of 5 Track B outer folds) and save the weights. Pretraining on all pure labels is fair for the mixture task: the baselines
also receive measured pure values, and no mixture label enters the pool (pool = pure rows only, see results/v6/gnn/README.md).
The encoder is frozen afterwards (embeddings only)."""
import sys, json, time
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd, torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from fuelblend.ml import dmpnn as D

R = Path(__file__).resolve().parents[3]
OUT = R / "results/v6/mixture"
CFG = dict(hidden=300, depth=3, dropout=0.2)
QUICK = "--quick" in sys.argv

def group(s):
    m = Chem.MolFromSmiles(s)
    if m.GetRingInfo().NumRings() == 0: return f"acyc_C{sum(a.GetAtomicNum()==6 for a in m.GetAtoms())}"
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))

pool = pd.read_csv(R / "results/v6/gnn/pool.csv")
comp = pool.drop_duplicates("inchikey")[["inchikey", "smiles"]].reset_index(drop=True)
comp["grp"] = comp.smiles.map(group)
tasks = sorted(pool.target.unique())
Yraw = pd.pivot_table(pool, index="inchikey", columns="target", values="value").reindex(comp.inchikey)[tasks].values
Mask = ~np.isnan(Yraw)
rng = np.random.default_rng(0); ug = comp.grp.unique(); vg = set(rng.choice(ug, int(.15 * len(ug)), replace=False))
va = np.where(comp.grp.isin(vg))[0]; tr = np.where(~comp.grp.isin(vg))[0]
mu = np.nan_to_num(np.nanmean(np.where(Mask[tr], Yraw[tr], np.nan), 0)); sd = np.nanstd(np.where(Mask[tr], Yraw[tr], np.nan), 0) + 1e-9
Z = np.nan_to_num((Yraw - mu) / sd)
feats = [D.featurize(s) for s in comp.smiles]
t0 = time.time()
ms = D.fit_ensemble([feats[i] for i in tr], Z[tr], Mask[tr], [feats[i] for i in va], Z[va], Mask[va],
                    n_models=2 if QUICK else 5, epochs=6 if QUICK else 150, seed=7, **CFG)
torch.save({"cfg": CFG, "tasks": tasks, "mu": mu, "sd": sd, "states": [m.state_dict() for m in ms]},
           OUT / ("encoder_quick.pt" if QUICK else "encoder.pt"))
print("done", f"{time.time()-t0:.0f}s", len(comp), "compounds", flush=True)
