"""Task 16: joint pure-Tm plus fuel MolPool for freezing (GPU), MC-dropout uncertainty, LOGO and leave-one-family-out on the 52 dev fuels.
  FUELBLEND_DEVICE=cuda python scripts/v7/freezing/molpool_joint.py [--quick] [--schemes logo,lofo]"""
import sys, time, importlib.util, warnings, os
from pathlib import Path
import numpy as np, pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts/v6/freezing"))
from fuelblend.ml import molpool as MP  # noqa: E402
from fuelblend.ml.curves import CurveLogger  # noqa: E402
from fuelblend.v7 import runner  # noqa: E402

QUICK = "--quick" in sys.argv
SCHEMES = (sys.argv[sys.argv.index("--schemes") + 1] if "--schemes" in sys.argv else "logo,lofo").split(",")
OUT = ROOT / ("scratch/molpool_joint_quick" if QUICK else "results/v7/freezing/molpool_joint")
EPOCHS, SEEDS, MC = (6, (0,), 3) if QUICK else (400, (0, 1, 2), 30)
PATIENCE, BATCH, LR, WD = 50, 64, 3e-3, 1e-2
LIVE = 13.8427
DEV = torch.device("cuda" if torch.cuda.is_available() and os.environ.get("FUELBLEND_DEVICE", "cpu") == "cuda" else "cpu")
warnings.filterwarnings("ignore")
spec = importlib.util.spec_from_file_location("v7_lofo_gp", ROOT / "scripts/v7/freezing/lofo_gp.py")
LG = importlib.util.module_from_spec(spec)
spec.loader.exec_module(LG)


def fuel_sets(dev):
    import common as C
    ds = C.dataset()
    ents = C.entities(ds)
    g1 = C.load("g1_v7mj", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    d5 = C.load("d5_v7mj", "scripts/v5/base/delta_freeze_t2_unifac.py")
    out = []
    for e in dev.entity:
        x, iks, smi, _ = d5.s5_inputs(ents[e], g1)
        out.append((list(smi), np.asarray(x, np.float32)))
    return out


def build_arrays(dev):
    sets = fuel_sets(dev)
    pool, _, _ = runner.load_track_b_pool()
    pure = pool[pool.tm_K.notna() & pool.is_hc.astype(bool)].reset_index(drop=True)
    enc = MP.load_encoder(ROOT / "results/v7/mixture/encoder.pt")
    smis = sorted({s for ss, _ in sets for s in ss} | set(pure.smiles))
    E = {}
    for i in range(0, len(smis), 256):
        chunk = smis[i:i + 256]
        for s, h in zip(chunk, MP.embed(enc, chunk)):
            E[s] = h
    Hf, Xf = MP.pad_sets([np.stack([E[s] for s in ss]) for ss, _ in sets], [x for _, x in sets])
    K = Hf.shape[1]
    emb = Hf.shape[2]
    Hp = np.zeros((len(pure), K, emb), np.float32)
    Xp = np.zeros((len(pure), K), np.float32)
    for i, s in enumerate(pure.smiles):
        Hp[i, 0] = E[s]
        Xp[i, 0] = 1.0
    return Hf, Xf, Hp, Xp, (pure.tm_K.to_numpy(float) - 273.15)


def train_one(H, X, A, r, w, val, seed, lg):
    torch.manual_seed(1000 + seed)
    net = MP.MolPool(H.shape[2], 16, 8, 1, hidden=32, dropout=0.1).to(DEV)
    opt = torch.optim.AdamW(net.parameters(), lr=LR, weight_decay=WD)
    tr = torch.as_tensor(np.where(~val)[0], device=DEV)
    va = torch.as_tensor(np.where(val)[0], device=DEV)
    g = torch.Generator().manual_seed(seed)
    wsum = lambda ww: ww.sum().clamp(min=1e-12)  # noqa: E731
    best, best_sd, bad = float("inf"), None, 0
    for ep in range(EPOCHS):
        t0 = time.time()
        net.train()
        perm = tr[torch.randperm(len(tr), generator=g).to(DEV)]
        tl, tn = 0.0, 0
        for i in range(0, len(perm), BATCH):
            idx = perm[i:i + BATCH]
            if len(idx) < 2:
                continue
            e = net(H[idx], X[idx], A[idx]) - r[idx]
            loss = (w[idx] * e ** 2).sum() / wsum(w[idx])
            opt.zero_grad()
            loss.backward()
            opt.step()
            tl += float(loss.detach())
            tn += 1
        net.eval()
        with torch.no_grad():
            ev = net(H[va], X[va], A[va]) - r[va]
            vl = float((w[va] * ev ** 2).sum() / wsum(w[va]))
        if lg is not None:
            lg.log(ep + 1, train_loss=tl / max(tn, 1), val_loss=vl, lr=LR, sec=time.time() - t0)
        if vl < best - 1e-6:
            best, bad, best_sd = vl, 0, {a: b.clone() for a, b in net.state_dict().items()}
        else:
            bad += 1
            if bad >= PATIENCE:
                break
    net.load_state_dict(best_sd)
    return net.eval()


def mc_predict(net, H, X, A):
    net.eval()
    net.phi[2].train()                                  # dropout active, nothing else changes
    with torch.no_grad():
        P = torch.stack([net(H, X, A) for _ in range(MC)]).cpu().numpy()
    net.eval()
    return P


def run_scheme(scheme, dev, fam, base, Hf, Xf, Hp, Xp, tmp, use_pure, lg):
    y = dev.freezing_C.to_numpy()
    groups = dev.group.to_numpy()
    mu, sd = np.full(len(y), np.nan), np.full(len(y), np.nan)
    Ht, Xt = torch.tensor(Hf, device=DEV), torch.tensor(Xf, device=DEV)
    Hpt, Xpt = torch.tensor(Hp, device=DEV), torch.tensor(Xp, device=DEV)
    for k, (tr, te) in enumerate(LG.folds_of(scheme, groups, fam)):
        rng = np.random.default_rng(100 + k)
        r_f = y[tr] - base[tr]
        sf = float(r_f.std()) or 1.0
        nv = max(3, int(0.15 * len(tr)))
        vf = np.zeros(len(tr), bool)
        vf[rng.choice(len(tr), nv, replace=False)] = True
        rows_H, rows_X = [Ht[tr]], [Xt[tr]]
        rows_r = [torch.tensor(r_f / sf, dtype=torch.float32, device=DEV)]
        rows_w = [torch.full((len(tr),), (0.5 if use_pure else 1.0) / len(tr), device=DEV)]
        rows_a = [torch.zeros(len(tr), 1, device=DEV)]
        val = [vf]
        if use_pure:
            sp = float(tmp.std())
            vp = rng.random(len(tmp)) < 0.15
            rows_H.append(Hpt)
            rows_X.append(Xpt)
            rows_r.append(torch.tensor((tmp - tmp.mean()) / sp, dtype=torch.float32, device=DEV))
            rows_w.append(torch.full((len(tmp),), 0.5 / len(tmp), device=DEV))
            rows_a.append(torch.ones(len(tmp), 1, device=DEV))
            val.append(vp)
        H, X, r, w, A = torch.cat(rows_H), torch.cat(rows_X), torch.cat(rows_r), torch.cat(rows_w), torch.cat(rows_a)
        v = np.concatenate(val)
        P = []
        for s in SEEDS:
            net = train_one(H, X, A, r, w, v, s, None if lg is None else lg.child(fold=f"{scheme}_{k}", seed=s))
            P.append(mc_predict(net, Ht[te], Xt[te], torch.zeros(len(te), 1, device=DEV)))
        P = np.concatenate(P) * sf
        mu[te], sd[te] = base[te] + P.mean(0), P.std(0)
        print(f"{scheme} fold {k + 1} done", flush=True)
    return mu, sd


def main():
    from scipy.stats import spearmanr
    OUT.mkdir(parents=True, exist_ok=True)
    dev, X, fam, ref = LG.build()
    base = dev.pred_e.clip(-100, -20).to_numpy()
    Hf, Xf, Hp, Xp, tmp = build_arrays(dev)
    gp = pd.read_csv(ROOT / "results/v7/freezing/gp/predictions.csv")
    assert (gp.entity.to_numpy() == dev.entity.to_numpy()).all()
    y = dev.freezing_C.to_numpy()
    lg = CurveLogger(OUT / "curves.csv", tb_dir=Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs")) / "freezing_joint")
    rows, cov, preds = [], [], dev[["entity", "group", "freezing_C"]].assign(family=fam)
    for scheme in SCHEMES:
        groups = dev.group.to_numpy()
        res = {}
        for name, up in (("molpool_joint", True), ("molpool_fuel_only", False)):
            mu, sd = run_scheme(scheme, dev, fam, base, Hf, Xf, Hp, Xp, tmp, up, lg.child(model=name))
            res[name] = (mu, sd)
            preds[f"{scheme}_{name}"], preds[f"{scheme}_{name}_mcsd"] = mu, sd
        comps = {m: gp[f"{scheme}_{m}"].to_numpy() for m in ("constant", "ridge", "gp", "baseline_b_unifac_svl", "limiting_species")}
        if scheme == "logo":
            comps["v6_lgbm_residual"] = gp["logo_v6_lgbm_residual"].to_numpy()
        allm = {**{k: v[0] for k, v in res.items()}, **comps}
        e = {m: np.abs(p - y) for m, p in allm.items()}
        for m, p in allm.items():
            row = dict(scheme=scheme, model=m, n=len(y), mae=float(e[m].mean()), bias=float((p - y).mean()),
                       mean_over=float(np.clip(p - y, 0, None).mean()), mean_under=float(np.clip(y - p, 0, None).mean()))
            for rm in ("constant", "ridge"):
                g_, lo, hi = LG.boot(e[rm] - e[m], groups)
                row.update({f"gain_vs_{rm}": g_, f"gain_vs_{rm}_lo": lo, f"gain_vs_{rm}_hi": hi})
            if m == "molpool_joint":
                g_, lo, hi = LG.boot(e["molpool_fuel_only"] - e[m], groups)
                row.update(gain_vs_fuel_only=g_, gain_vs_fuel_only_lo=lo, gain_vs_fuel_only_hi=hi)
            if m in res:
                row["spearman_mcsd_abs_err"] = float(spearmanr(res[m][1], e[m])[0])
            rows.append(row)
        for f_ in np.unique(fam):
            mk = fam == f_
            for m in ("molpool_joint", "molpool_fuel_only", "constant", "ridge"):
                cov.append(dict(scheme=scheme, family=f_, model=m, n_groups=int(len(np.unique(groups[mk]))),
                                coverage_live_margin=float(((y[mk] - allm[m][mk]) <= LIVE).mean()), mae=float(e[m][mk].mean())))
    lg.close()
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False, float_format="%.5g")
    pd.DataFrame(cov).to_csv(OUT / "coverage_by_family.csv", index=False, float_format="%.5g")
    preds.to_csv(OUT / "predictions.csv", index=False, float_format="%.5g")
    print(pd.DataFrame(rows).round(3).to_string())


if __name__ == "__main__":
    main()
