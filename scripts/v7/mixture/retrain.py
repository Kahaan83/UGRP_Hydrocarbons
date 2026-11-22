"""v7 Task 11: MolPool mixture heads (Track C) with nested Optuna on the GPU. See PREDECLARATION.md.
Usage: FUELBLEND_DEVICE=cuda python scripts/v7/mixture/retrain.py [--quick]   (run pretrain_encoder.py first)"""
import sys, json, time, os, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/v7/mixture")
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from fuelblend.ml import molpool as MP
from fuelblend.v7 import guard, runner
import _mp

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/mixture" / ("retrain_quick" if QUICK else "retrain")
OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs"))
NT, TO = (2, 120) if QUICK else (60, 1800)
V6 = R / "results/v6/mixture"
rows = pd.read_csv(V6 / "rows.csv")
dev_all = rows[rows.role == "dev"]
guard.assert_dev_only(dev_all.system_id.unique(), locks=("mixture_lock", "mixture_holdout51"))
v6 = pd.read_csv(V6 / "dev_oof.csv")
enc = MP.load_encoder(R / "results/v7/mixture" / ("encoder_quick.pt" if QUICK else "encoder.pt"))
COMP = {"octane": ["lin_molar", "corrected_pure", "desc_gbm"], "boehm": ["physics", "t_bias", "desc_gbm"]}
t0 = time.time()
out_rows, summ = [], []
for (head, prop), d in dev_all.groupby(["head", "prop"]):
    d = d.reset_index(drop=True)
    kind = "boehm" if head == "boehm" else "octane"
    ref = v6[(v6["head"] == head) & (v6["prop"] == prop)].reset_index(drop=True)
    assert len(ref) == len(d) and np.allclose(ref.y, d.y), "v6 OOF rows do not match"
    sm = [json.loads(s) for s in d.smiles]
    x = [np.array(json.loads(s), float) for s in d.x]
    allsmi = sorted({s for r in sm for s in r})
    E = dict(zip(allsmi, MP.embed(enc, allsmi)))
    H, X = MP.pad_sets([np.stack([E[s] for s in r]) for r in sm], x)
    g = d.group.to_numpy()
    base, y = d.base.to_numpy(float), d.y.to_numpy(float)
    A = (1000.0 / d.T_K.to_numpy(float))[:, None] if kind == "boehm" else None
    folds = np.full(len(d), -1)
    for f, (_, te) in enumerate(GroupKFold(min(5, len(np.unique(g)))).split(np.arange(len(d)), groups=g)):
        folds[te] = f
    oof, res = _mp.tune(H, X, A, y - base, g, folds, prefix=f"v7_mix_{head}_{prop}" + ("_quick" if QUICK else ""), out_dir=OUT / f"{head}_{prop}",
                        nt=NT, to=TO, tb=RUNS / "mixture" / "retrain" / f"{head}_{prop}", quick=QUICK)
    p7 = base + oof
    ev = np.abs(p7 - y)
    comp_mae = {k: float(np.abs(ref[f"pred_{k}"] - y).mean()) for k in COMP[kind]}
    best = min(comp_mae, key=comp_mae.get)
    e_best, e_v1 = np.abs(ref[f"pred_{best}"].to_numpy() - y), np.abs(ref.pred_molpool.to_numpy() - y)
    gb, glo, ghi, pb = _mp.boot_p(e_best - ev, g)
    g1, l1, h1, p1 = _mp.boot_p(e_v1 - ev, g)
    summ.append(dict(head=head, prop=prop, n=len(d), groups=int(len(np.unique(g))), mae_v7=ev.mean(), mae_v1=e_v1.mean(), best_comparator=best,
                     mae_best=e_best.mean(), **{f"mae_{k}": v for k, v in comp_mae.items()}, gain_vs_best=gb, best_lo=glo, best_hi=ghi, p_best=pb,
                     gain_vs_v1=g1, v1_lo=l1, v1_hi=h1, p_v1=p1, best_params=res.best_params.tolist()))
    o = d[["head", "prop", "group", "system_id", "n_comp", "T_K", "y", "base"]].copy()
    o["pred_v7"] = p7
    o["pred_v1"] = ref.pred_molpool.to_numpy()
    o[f"pred_{best}"] = ref[f"pred_{best}"].to_numpy()
    out_rows.append(o)
    print(head, prop, f"{time.time() - t0:.0f}s", flush=True)
S = pd.DataFrame(summ)
S["p_holm_best"] = _mp.holm(S.p_best.to_numpy())
S["p_holm_v1"] = _mp.holm(S.p_v1.to_numpy())
S["beats_best"] = (S.p_holm_best < 0.10) & (S.gain_vs_best > 0)
S["beats_v1"] = (S.p_holm_v1 < 0.10) & (S.gain_vs_v1 > 0)
S["best_params"] = S.best_params.map(json.dumps)
runner.write_results(OUT, summary=S, oof=pd.concat(out_rows, ignore_index=True),
                     meta=dict(task="11_mixture_molpool", n_trials=NT, timeout_s=TO, total_s=time.time() - t0))
print(S.drop(columns="best_params").round(3).to_string())
