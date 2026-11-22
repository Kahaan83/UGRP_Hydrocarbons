"""v7 Task 7: in-house D-MPNN vs Chemprop 2 on flash_C and ysi (GPU venv). TensorBoard logs go to C:/fuelblend-gpu/runs (scripts/v7/monitor.ps1)."""
import sys, json, time, tempfile, warnings, logging
warnings.filterwarnings("ignore"); logging.getLogger("lightning.pytorch").setLevel(logging.ERROR)
sys.path.insert(0, "src")
from pathlib import Path
import numpy as np, pandas as pd, torch
import lightning.pytorch as pl
from lightning.pytorch.callbacks import EarlyStopping, ModelCheckpoint
from lightning.pytorch.loggers import TensorBoardLogger
from chemprop import data, featurizers, models, nn
from fuelblend.ml import dmpnn as D
from fuelblend.ml.curves import CurveLogger
from fuelblend.v7 import runner

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/track_b" / ("chemprop_parity_quick" if QUICK else "chemprop_parity"); OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path("C:/fuelblend-gpu/runs/chemprop_parity" + ("_quick" if QUICK else ""))
NM, EP, PAT = (2, 6, 3) if QUICK else (5, 150, 25)
TARGETS = ["flash_C", "ysi"]
comp, folds, tasks = runner.load_track_b_pool()
g = comp.grp.values; smi = comp.smiles.values
feat = featurizers.SimpleMoleculeMolGraphFeaturizer()
dfeats = {}
def df_(i):
    if i not in dfeats: dfeats[i] = D.featurize(smi[i])
    return dfeats[i]

def chemprop_fit_predict(tr2, va, te, y, seed):
    pl.seed_everything(seed, workers=True)
    mk = lambda idx: data.MoleculeDataset([data.MoleculeDatapoint.from_smi(smi[i], [y[i]]) for i in idx], feat)
    dtr, dva, dte = mk(tr2), mk(va), mk(te)
    sc = dtr.normalize_targets(); dva.normalize_targets(sc)
    ltr = data.build_dataloader(dtr, batch_size=32, num_workers=0, shuffle=True, seed=seed)
    lva = data.build_dataloader(dva, batch_size=64, num_workers=0, shuffle=False)
    lte = data.build_dataloader(dte, batch_size=64, num_workers=0, shuffle=False)
    mp = nn.BondMessagePassing(d_h=200, depth=4, dropout=0.1)
    ffn = nn.RegressionFFN(input_dim=mp.output_dim, hidden_dim=200, n_layers=1, output_transform=nn.UnscaleTransform.from_standard_scaler(sc))
    m = models.MPNN(mp, nn.MeanAggregation(), ffn, batch_norm=False, metrics=[nn.metrics.MAE()])
    with tempfile.TemporaryDirectory() as td:
        ck = ModelCheckpoint(dirpath=td, monitor="val_loss", mode="min", save_top_k=1)
        tr = pl.Trainer(logger=TensorBoardLogger(str(RUNS), name=f"cp_{TARGET}", version=f"f{FOLD}_s{seed}"), enable_checkpointing=True,
                        callbacks=[ck, EarlyStopping("val_loss", patience=PAT)], accelerator="gpu", devices=1, max_epochs=EP,
                        enable_progress_bar=False, enable_model_summary=False)
        tr.fit(m, ltr, lva)
        best = models.MPNN.load_from_checkpoint(ck.best_model_path)
        p = np.concatenate([b.numpy().ravel() for b in tr.predict(best, lte)])
    return p

rows, oofs = [], {}
lg = CurveLogger(OUT / "inhouse_curves.csv", tb_dir=RUNS / "inhouse")
t0 = time.time()
for TARGET in TARGETS:
    y = comp[TARGET].values; has = ~np.isnan(y)
    p_in = np.full(len(comp), np.nan); p_cp = np.full(len(comp), np.nan)
    for FOLD in range(5):
        tr = np.where((folds != FOLD) & has)[0]; te = np.where((folds == FOLD) & has)[0]
        ug = np.unique(g[tr]); rng = np.random.default_rng(FOLD)
        vg = set(rng.choice(ug, max(1, int(.15 * len(ug))), replace=False))
        va = np.array([i for i in tr if g[i] in vg]); tr2 = np.array([i for i in tr if g[i] not in vg])
        mu, sd = y[tr2].mean(), y[tr2].std() + 1e-9
        Z = ((y - mu) / sd)[:, None]; Mk = has[:, None].astype(float)
        ms = D.fit_ensemble([df_(i) for i in tr2], np.nan_to_num(Z[tr2]), Mk[tr2], [df_(i) for i in va], np.nan_to_num(Z[va]), Mk[va],
                            n_models=NM, epochs=EP, patience=PAT, hidden=200, depth=4, dropout=0.1, seed=FOLD,
                            logger=lg.child(model=f"inhouse_{TARGET}", fold=FOLD))
        p_in[te] = D.predict(ms, [df_(i) for i in te])[0][:, 0] * sd + mu
        p_cp[te] = np.mean([chemprop_fit_predict(tr2, va, te, y, FOLD * 100 + s) for s in range(NM)], 0)
        print(TARGET, "fold", FOLD, f"{time.time()-t0:.0f}s", flush=True)
    m = has & ~np.isnan(p_in) & ~np.isnan(p_cp)
    lo, hi = runner.paired_boot_gain(p_in[m] - y[m], p_cp[m] - y[m], g[m])
    rows.append(dict(target=TARGET, n=int(m.sum()), mae_inhouse=np.abs(p_in[m] - y[m]).mean(), mae_chemprop=np.abs(p_cp[m] - y[m]).mean(),
                     gain_lo=lo, gain_hi=hi, parity=bool(lo <= 0 <= hi), inhouse_worse=bool(hi < 0)))
    oofs[TARGET + "_obs"] = y; oofs[TARGET + "_inhouse"] = p_in; oofs[TARGET + "_chemprop"] = p_cp
lg.close()
res = pd.DataFrame(rows)
runner.write_results(OUT, summary=res, oof=pd.concat([comp[["inchikey", "smiles", "grp"]], pd.DataFrame(oofs)], axis=1),
                     meta=dict(task="7_chemprop_parity", n_models=NM, epochs=EP, patience=PAT, total_s=time.time() - t0,
                               n_locked_ysi=comp.attrs.get("n_locked_ysi")))
print(res.round(3).to_string())
