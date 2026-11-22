"""Task 12: v5 set models (pds, gnn, st) retrained on dev rows with a 20 % group-held-out validation split, early stopping and curves.
Imports scripts/v5/base/delta_ws4_mixture_ml.py by path (unchanged) and replaces its fit_predict; its dev stage then runs as in v5.
  python scripts/v7/mixture/retrain_set_models.py [--cands pds,gnn,st] [--epochs 500] [--seeds 5] [--out results/v7/set_models] [--tag _quick]"""
import sys, os, hashlib, importlib.util, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts/v7/cpu"))
import numpy as np
import torch
import _trees
from fuelblend.ml import set_models as sm
from fuelblend.ml.curves import CurveLogger

VAL_FRAC, PATIENCE = 0.2, 50
spec = importlib.util.spec_from_file_location("v5_ws4", ROOT / "scripts/v5/base/delta_ws4_mixture_ml.py")
ws4 = importlib.util.module_from_spec(spec)
sys.modules["v5_ws4"] = ws4
spec.loader.exec_module(ws4)

OUT = ROOT / "results/v7/set_models"
for i, a in enumerate(sys.argv):
    if a == "--out":
        OUT = Path(sys.argv[i + 1])
OUT.mkdir(parents=True, exist_ok=True)
CL = CurveLogger(OUT / "curves.csv", tb_dir=Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs")) / "set_models")
FOLD_SEED = {"n": 0}


def fit_predict(cand, DS, train_rows, pred_rows, args, tag):
    sc = ws4.Scaler(DS, train_rows)
    bs, ys, sysl, props, ystd = [], [], [], [], []
    for p, (D, tr) in enumerate(zip(DS, train_rows)):
        sd = float(np.std(D["y"][tr])) + 1e-9
        ystd.append(sd)
        bs.append(sc.batch(D, tr, p))
        ys.append(D["y"][tr] / sd)
        sysl.append(D["sys"][tr])
        props.append(np.full(len(tr), p))
    cat = lambda name: torch.cat([getattr(b, name) for b in bs])  # noqa: E731
    sysarr, proparr = np.concatenate(sysl), np.concatenate(props)
    y = np.concatenate(ys)
    # 20 % of the training systems (shared across properties) held out for early stopping
    usys = np.unique(sysarr)
    seed = int(hashlib.md5(tag.encode()).hexdigest()[:8], 16)
    vs = set(np.random.default_rng(seed).choice(usys, max(1, int(VAL_FRAC * len(usys))), replace=False))
    isv = np.array([s in vs for s in sysarr])
    tri, vai = np.where(~isv)[0], np.where(isv)[0]
    full = sm.SetBatch(cat("F"), cat("X"), cat("M"), cat("invT"), cat("mol"), cat("prop"))
    data, vdata = full.take(torch.as_tensor(tri)), full.take(torch.as_tensor(vai))
    w_tr = sm.system_macro_weights(sysarr[tri], proparr[tri])
    w_va = sm.system_macro_weights(sysarr[vai], proparr[vai])
    t0 = time.time()
    ens = sm.train_ensemble(ws4.make_factory(cand, DS[0], args), data, y[tri], w_tr, seeds=list(range(args.seeds)), epochs=args.epochs,
                            batch_size=args.batch, lr=args.lr, wd=args.wd, device=args.device,
                            curve=CL.child(cand=cand, tag=tag), val=(vdata, y[vai], w_va), patience=PATIENCE)
    out = []
    for p, (D, pr) in enumerate(zip(DS, pred_rows)):
        if len(pr) == 0:
            out.append((np.zeros(0), np.zeros(0)))
            continue
        mu, sdv = sm.predict_ensemble(ens, sc.batch(D, pr, p), device=args.device)
        out.append((mu * ystd[p], sdv * ystd[p]))
    ws4.log(f"{tag} {cand}: train {len(tri)} val {len(vai)} rows, {time.time() - t0:.0f}s")
    return out


ws4.fit_predict = fit_predict
from m_viscosity_pair import _redirect_removed_v3   # v4 leaderboard helper lives in the removed scripts/v3 (git copy in scratch/legacy)
_redirect_removed_v3()
if "--stage" not in sys.argv:
    sys.argv += ["--stage", "dev"]
assert "locked" not in sys.argv, "locked systems are not scored in v7 Task 12"
for flag, val in (("--device", "cuda"), ("--epochs", "500"), ("--out", str(OUT))):
    if flag not in sys.argv:
        sys.argv += [flag, val]
h0 = _trees.tree_hash("src", "results/v5/BASE")
ws4.main()
CL.close()
assert _trees.tree_hash("src", "results/v5/BASE") == h0, "src or results/v5/BASE changed during the run"
