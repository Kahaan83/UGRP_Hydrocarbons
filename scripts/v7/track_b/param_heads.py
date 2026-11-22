"""Task 18: D-MPNN equation-parameter heads (Antoine, Andrade), scaffold GroupKFold on molecules, comparators and flash recompute.
  FUELBLEND_DEVICE=cuda python scripts/v7/track_b/param_heads.py [--quick]"""
import sys, json, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import torch
from scipy.optimize import curve_fit
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")
from fuelblend.ml import dmpnn as D  # noqa: E402
from fuelblend.ml import param_heads as PH  # noqa: E402
from fuelblend.ml.curves import CurveLogger  # noqa: E402
from fuelblend.v7 import guard, runner  # noqa: E402
from fuelblend.baselines.flash import predict_flash_le_chatelier  # noqa: E402

QUICK = "--quick" in sys.argv
OUT = ROOT / ("scratch/param_heads_quick" if QUICK else "results/v7/track_b/param_heads")
EPOCHS, SEEDS, PATIENCE = (4, (0,), 4) if QUICK else (300, (0, 1, 2), 40)
rng = np.random.default_rng(11)


def load_points():
    d = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=["property", "kind", "is_hc", "T_K", "comp_inchikeys", "value", "system_id"])
    d = d[(d.kind == "pure") & d.is_hc.astype(bool)].copy()
    d["ik"] = d.comp_inchikeys.map(lambda s: json.loads(s)[0])
    comp = pd.read_csv(ROOT / "data/unified/compounds.csv").set_index("inchikey")
    out = {}
    for prop, key in (("vapor_pressure_kPa", "P"), ("viscosity_dyn_mPas", "mu")):
        p = d[(d.property == prop) & d.T_K.notna() & (d.value > 0)].copy()
        p["y"] = np.log10(p.value) if key == "P" else np.log(p.value)
        n = p.groupby("ik").T_K.nunique()
        out[key] = p[p.ik.isin(n[n >= 3].index)][["ik", "T_K", "y"]].reset_index(drop=True)
    mols = sorted(set(out["P"].ik) | set(out["mu"].ik))
    mols = [m for m in mols if m in comp.index]
    df = pd.DataFrame({"inchikey": mols, "smiles": [comp.smiles[m] for m in mols]})
    df["grp"] = df.smiles.map(runner._group)
    for k in out:
        out[k] = out[k][out[k].ik.isin(mols)].reset_index(drop=True)
    return df, out, comp, d


def descr(smi):
    m = Chem.MolFromSmiles(smi)
    return [sum(a.GetAtomicNum() == 6 for a in m.GetAtoms()), m.GetNumHeavyAtoms(), m.GetRingInfo().NumRings(), sum(a.GetIsAromatic() for a in m.GetAtoms())]


def pts(out, key, mol_index):
    p = out[key]
    i = p.ik.map(mol_index)
    ok = i.notna()
    return i[ok].astype(int).to_numpy(), p.T_K[ok].to_numpy(float), p.y[ok].to_numpy(float)


def train_one(feats, tr, va, PT, seed, lg):
    torch.manual_seed(seed)
    net = PH.ParamHeadDMPNN().to(D.DEVICE)
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-5)
    r = np.random.default_rng(seed)
    dev = D.DEVICE

    def loss_for(par, mols_sel, mol_pos):
        tot, cnt = 0.0, 0
        for key, (mi, T, y) in PT.items():
            m = np.isin(mi, mols_sel)
            if not m.any():
                continue
            pos = torch.as_tensor(np.array([mol_pos[j] for j in mi[m]]), device=dev)
            Tt, yt = torch.tensor(T[m], dtype=torch.float32, device=dev), torch.tensor(y[m], dtype=torch.float32, device=dev)
            pred = PH.log10_p(par[pos], Tt) if key == "P" else PH.ln_mu(par[pos], Tt)
            tot = tot + ((pred - yt) ** 2).sum()
            cnt += int(m.sum())
        return tot / max(cnt, 1), cnt

    best, best_sd, bad = 1e9, None, 0
    vb = D.collate([feats[j] for j in va])
    vpos = {int(j): k for k, j in enumerate(va)}
    for ep in range(EPOCHS):
        t0 = time.time()
        net.train()
        perm = r.permutation(tr)
        tl, tn = 0.0, 0
        for i in range(0, len(perm), 32):
            idx = perm[i:i + 32]
            if len(idx) < 2:
                continue
            par = net.params(D.collate([feats[j] for j in idx]))
            l, c = loss_for(par, idx, {int(j): k for k, j in enumerate(idx)})
            if c == 0:
                continue
            opt.zero_grad(); l.backward(); opt.step()
            tl += float(l.detach()); tn += 1
        net.eval()
        with torch.no_grad():
            vl = float(loss_for(net.params(vb), va, vpos)[0])
        if lg is not None:
            lg.log(ep + 1, train_loss=tl / max(tn, 1), val_loss=vl, lr=1e-3, sec=time.time() - t0)
        if vl < best - 1e-6:
            best, bad, best_sd = vl, 0, {a: b.clone() for a, b in net.state_dict().items()}
        else:
            bad += 1
            if bad >= PATIENCE:
                break
    net.load_state_dict(best_sd)
    return net.eval()


def fit_pooled(PT, tr_mask):
    res = {}
    mi, T, y = PT["P"]
    m = tr_mask[mi]
    f = lambda T_, A, B, C: A - B / (T_ + C)  # noqa: E731
    res["P"] = curve_fit(f, T[m], y[m], p0=[6, 1500, -50], bounds=([-5, 10, -150], [20, 10000, 0]), maxfev=20000)[0]
    mi, T, y = PT["mu"]
    m = tr_mask[mi]
    b, a = np.polyfit(1 / T[m], y[m], 1)
    res["mu"] = (a, b)
    return res


def ridge_fit(F, T, y, alpha=1.0):
    X = np.column_stack([np.ones(len(T)), 1 / T, F, F / T[:, None]])
    mu, sd = X[:, 1:].mean(0), X[:, 1:].std(0) + 1e-9
    Z = np.column_stack([np.ones(len(T)), (X[:, 1:] - mu) / sd])
    w = np.linalg.solve(Z.T @ Z + alpha * np.diag([0] + [1] * (Z.shape[1] - 1)), Z.T @ y)
    return lambda F_, T_: np.column_stack([np.ones(len(T_)), (np.column_stack([1 / T_, F_, F_ / T_[:, None]]) - mu) / sd]) @ w


def boot(d, g, n=4000):
    ug = np.unique(g)
    pm = np.array([d[g == k].mean() for k in ug])
    b = np.array([pm[rng.integers(0, len(ug), len(ug))].mean() for _ in range(n)])
    return float(pm.mean()), float(np.percentile(b, 5)), float(np.percentile(b, 95))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, out, comp, dd = load_points()
    mol_index = {m: i for i, m in enumerate(df.inchikey)}
    feats = [D.featurize(s) for s in df.smiles]
    Fd = np.array([descr(s) for s in df.smiles], float)
    folds = runner.scaffold_folds(df, 5)
    PT = {k: pts(out, k, mol_index) for k in out}
    lg = CurveLogger(OUT / "curves.csv", tb_dir=Path(__import__("os").environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs")) / "param_heads")
    n = len(df)
    head_par = np.full((n, 5), np.nan)
    pred = {k: {m: np.full(len(PT[k][0]), np.nan) for m in ("head", "pooled", "ridge", "refit")} for k in PT}
    pooled_par = {}
    for f in range(5):
        te = np.where(folds == f)[0]
        trm = np.where(folds != f)[0]
        nv = max(3, int(0.15 * len(trm)))
        va = r_choice = rng.choice(trm, nv, replace=False)
        tr = np.setdiff1d(trm, va)
        P = []
        for s in SEEDS:
            net = train_one(feats, tr, va, {k: (v[0], v[1], v[2]) for k, v in PT.items()}, 100 * f + s, lg.child(model="param_head", fold=f, seed=s))
            with torch.no_grad():
                P.append(net.params(D.collate([feats[j] for j in te])).cpu().numpy())
        head_par[te] = np.mean(P, axis=0)
        trmask = np.zeros(n, bool); trmask[trm] = True
        pool = fit_pooled(PT, trmask)
        pooled_par[f] = pool
        for key in PT:
            mi, T, y = PT[key]
            ptr = trmask[mi]
            pte = np.isin(mi, te)
            rf = ridge_fit(Fd[mi[ptr]], T[ptr], y[ptr])
            pos = {int(j): k for k, j in enumerate(te)}
            for j in np.where(pte)[0]:
                par = head_par[mi[j]]
                pred[key]["head"][j] = (par[0] - par[1] / (T[j] + par[2])) if key == "P" else par[3] + par[4] / T[j]
                pred[key]["pooled"][j] = (pool["P"][0] - pool["P"][1] / (T[j] + pool["P"][2])) if key == "P" else pool["mu"][0] + pool["mu"][1] / T[j]
            pred[key]["ridge"][pte] = rf(Fd[mi[pte]], T[pte])
    # per-molecule refit reference
    for key in PT:
        mi, T, y = PT[key]
        for j in range(n):
            m = mi == j
            if m.sum() < (4 if key == "P" else 3):
                continue
            try:
                if key == "P":
                    p, _ = curve_fit(lambda T_, A, B, C: A - B / (T_ + C), T[m], y[m], p0=[6, 1500, -50], bounds=([-5, 10, -150], [20, 10000, 0]), maxfev=20000)
                    pred[key]["refit"][m] = p[0] - p[1] / (T[m] + p[2])
                else:
                    b, a = np.polyfit(1 / T[m], y[m], 1)
                    pred[key]["refit"][m] = a + b / T[m]
            except Exception:
                pass
    rows = []
    for key in PT:
        mi, T, y = PT[key]
        g = mi
        err = {m: np.abs(pred[key][m] - y) for m in pred[key]}
        for m in ("head", "pooled", "ridge", "refit"):
            row = dict(prop="log10_P" if key == "P" else "ln_mu", model=m, n_points=len(y), n_molecules=int(len(np.unique(g))),
                       mae=float(np.nanmean(err[m])))
            if m != "head":
                continue
            for comp_ in ("pooled", "ridge"):
                gm = boot(err[comp_] - err["head"], g)
                row.update({f"gain_vs_{comp_}": gm[0], f"gain_vs_{comp_}_lo": gm[1], f"gain_vs_{comp_}_hi": gm[2]})
            rows.append(row)
        for m in ("pooled", "ridge", "refit"):
            rows.append(dict(prop="log10_P" if key == "P" else "ln_mu", model=m, n_points=len(y), n_molecules=int(len(np.unique(g))), mae=float(np.nanmean(err[m]))))
    pd.DataFrame(rows).to_csv(OUT / "summary.csv", index=False, float_format="%.5g")
    pd.DataFrame(dict(inchikey=df.inchikey, smiles=df.smiles, fold=folds, A=head_par[:, 0], B=head_par[:, 1], C=head_par[:, 2], a=head_par[:, 3], b=head_par[:, 4])).to_csv(OUT / "head_params.csv", index=False, float_format="%.5g")
    # flash recompute
    fl = dd[(dd.property == "flash_point_C")].copy()
    fl = fl[fl.ik.isin(df.inchikey)]
    loc = set()
    for nm, sp in guard.LOCKS.items():
        if sp.id_kind == "system_id":
            loc |= {guard.canon(i) for i in guard.locked_ids(nm)}
    fl = fl[~fl.system_id.isin(loc)]
    frows, n_c10 = [], 0
    for r in fl.itertuples():
        c = comp.loc[r.ik]
        if c.cls in ("n_paraffin", "iso_paraffin") and c.n_carbon >= 10:
            n_c10 += 1
            continue
        if not np.isfinite(c.lfl_vol_pct) or not np.isfinite(c.antoine_A):
            continue
        j = mol_index[r.ik]
        f = int(folds[j])
        pl = pooled_par[f]["P"]

        def fp(ant):
            T_, st = predict_flash_le_chatelier([1.0], [tuple(map(float, ant))], [float(c.lfl_vol_pct)], lfl_temp=True)
            return T_ if st == "ok" else np.nan
        frows.append(dict(ik=r.ik, group=df.grp[j], measured=r.value, head=fp(head_par[j, :3]), pooled=fp(pl), registry=fp([c.antoine_A, c.antoine_B, c.antoine_C])))
    F_ = pd.DataFrame(frows)
    F_.to_csv(OUT / "flash_rows.csv", index=False, float_format="%.4f")
    ok = F_[["head", "pooled", "registry"]].notna().all(axis=1)
    F2 = F_[ok]
    fs = dict(n_rows=int(len(F2)), n_molecules=int(F2.ik.nunique()), n_excluded_c10_paraffin=int(n_c10), n_failed=int((~ok).sum()))
    for m in ("head", "pooled", "registry"):
        e = F2[m] - F2.measured
        fs[m] = dict(mae=float(e.abs().mean()), mean_over=float(e.clip(lower=0).mean()), max_over=float(e.max()))
    gm = boot(np.abs(F2.pooled - F2.measured).to_numpy() - np.abs(F2["head"] - F2.measured).to_numpy(), F2.ik.to_numpy())
    fs["head_gain_vs_pooled"] = dict(mean=gm[0], lo=gm[1], hi=gm[2])
    (OUT / "flash_summary.json").write_text(json.dumps(fs, indent=1))
    lg.close()
    print(pd.DataFrame(rows).round(3).to_string()); print(json.dumps(fs, indent=1))


if __name__ == "__main__":
    main()
