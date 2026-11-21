"""Track C: score the MolPool mixture head and the comparators on the declared rows (results/v6/mixture/README.md).

  python scripts/v6/mixture/run_cv.py --stage dev      dev groups only, GroupKFold out-of-fold  -> dev_oof.csv
  python scripts/v6/mixture/run_cv.py --stage locked   fit on all dev, score locked groups ONCE -> locked_pred.csv
  python scripts/v6/mixture/run_cv.py --stage summary  paired group bootstrap + gate           -> readings.json, summary.csv
Add --quick for a smoke run (writes *_quick files)."""
import sys, json, time, argparse
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold
RDLogger.DisableLog("rdApp.*")
R = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(R / "src")); sys.path.insert(0, str(R / "scripts/external_ingest"))
from fuelblend.ml import molpool as MP  # noqa: E402

OUT = R / "results/v6/mixture"
CONFIGS = [None, dict(d=8, wd=1e-2), dict(d=8, wd=1e-1), dict(d=24, wd=1e-2), dict(d=24, wd=1e-1)]   # None = baseline (no correction)
ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["dev", "locked", "summary"])
ap.add_argument("--quick", action="store_true")
ap.add_argument("--v2", action="store_true", help="post-hoc dev-only optimisation variant (locked groups were already scored once)")
args = ap.parse_args()
Q = "_quick" if args.quick else ("_v2" if args.v2 else "")
if args.v2:
    assert args.stage != "locked", "locked groups are spent; v2 is dev-only"
    CONFIGS = [None] + [dict(d=d, wd=wd, epochs=ep, lr=lr) for d in (8, 24) for wd in (1e-2, 1e-1, 1.0) for ep, lr in ((60, 3e-3), (300, 1e-3))]
SHRINKS = (0.5, 1.0) if args.v2 else (1.0,)
EPOCHS, SEEDS, INNER_EP = (15, (0,), 15) if args.quick else (150, (0, 1, 2, 3, 4), 120)
COMPARATORS = {"octane": ["lin_molar", "corrected_pure", "desc_gbm"], "boehm": ["physics", "t_bias", "desc_gbm"]}


def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    return np.array([Descriptors.MolWt(m) / 100, rdMolDescriptors.CalcNumAromaticRings(m), rdMolDescriptors.CalcFractionCSP3(m),
                     sum(a.GetAtomicNum() == 8 for a in m.GetAtoms()), rdMolDescriptors.CalcNumRings(m),
                     sum(a.GetAtomicNum() == 6 for a in m.GetAtoms()) / 10], float)


def pure_values(head, prop):
    """inchikey -> measured pure value used by the linear molar baseline of this head (same source as build_tables.py)."""
    if head == "kaust":
        from kaust_registry import COMPONENTS, pure_props
        from rdkit.Chem import MolToInchiKey
        pp = pure_props()
        return {MolToInchiKey(Chem.MolFromSmiles(c["smiles"])): pp[k][prop] for k, c in COMPONENTS.items()}
    if head == "leenhouts":
        pool = pd.read_csv(R / "results/v6/gnn/pool.csv")
        return pool[pool.target == "mon"].set_index("inchikey").value.to_dict()
    return {}


class Head:
    def __init__(self, df, enc):
        self.df = df.reset_index(drop=True)
        self.ks = [json.loads(s) for s in self.df.comps]
        self.sm = [json.loads(s) for s in self.df.smiles]
        self.x = [np.array(json.loads(s), float) for s in self.df.x]
        allsmi = sorted({s for r in self.sm for s in r})
        E = dict(zip(allsmi, MP.embed(enc, allsmi)))
        Dsc = {s: descriptors(s) for s in allsmi}
        self.H, self.X = MP.pad_sets([np.stack([E[s] for s in r]) for r in self.sm], self.x)
        self.D = np.stack([(x[:, None] * np.stack([Dsc[s] for s in r])).sum(0) for r, x in zip(self.sm, self.x)])
        self.T = self.df.T_K.to_numpy(float)
        self.base = self.df.base.to_numpy(float)
        self.y = self.df.y.to_numpy(float)
        self.g = self.df.group.to_numpy()


def corrected_pure(h, tr, te, pure, include_te=False):
    """Linear molar with per-component effective values c = p + (X'X + lam I)^-1 X'(y - Xp); lam by inner GroupKFold(3)."""
    comps = sorted({k for i in (list(tr) + list(te) if include_te else tr) for k in h.ks[i] if k in pure})
    pos = {k: j for j, k in enumerate(comps)}
    p = np.array([pure[k] for k in comps])

    def design(idx):
        M = np.zeros((len(idx), len(comps)))
        cov = np.zeros(len(idx))
        for r, i in enumerate(idx):
            for k, xi in zip(h.ks[i], h.x[i]):
                if k in pos:
                    M[r, pos[k]] += xi
                    cov[r] += xi
        return M, cov

    def solve(idx, lam):
        M, cov = design(idx)
        Mn = M / np.maximum(cov, 1e-9)[:, None]
        return p + np.linalg.solve(Mn.T @ Mn + lam * np.eye(len(comps)), Mn.T @ (h.y[idx] - Mn @ p))

    def predict(idx, c):
        M, cov = design(idx)
        return (M / np.maximum(cov, 1e-9)[:, None]) @ c

    lams = [0.1, 1, 10, 100, 1000]
    ng = len(np.unique(h.g[tr]))
    best = (1e18, 10)
    if ng >= 3:
        for lam in lams:
            e = []
            for a, b in GroupKFold(min(3, ng)).split(tr, groups=h.g[tr]):
                c = solve(tr[a], lam)
                e.append(np.abs(predict(tr[b], c) - h.y[tr[b]]))
            m = np.concatenate(e).mean()
            if m < best[0]:
                best = (m, lam)
    c = solve(tr, best[1])
    return predict(te, c)


def gbm_residual(h, tr, te, extra_T):
    def f(idx):
        cols = [h.D[idx], h.base[idx, None]] + ([h.T[idx, None]] if extra_T else [])
        return np.hstack(cols)
    m = GradientBoostingRegressor(n_estimators=45, max_depth=3, learning_rate=0.08, random_state=42)
    m.fit(f(tr), h.y[tr] - h.base[tr])
    return h.base[te] + m.predict(f(te))


def candidate(h, tr, te, use_T, rec):
    r = h.y[tr] - h.base[tr]
    s = float(r.std()) if r.std() > 1e-9 else 1.0
    aux_all = (1000.0 / h.T)[:, None] if use_T else None
    if aux_all is not None:
        am, asd = aux_all[tr].mean(0), aux_all[tr].std(0) + 1e-9
        aux_all = (aux_all - am) / asd
    sub = lambda idx: (h.H[idx], h.X[idx], None if aux_all is None else aux_all[idx])
    ng = len(np.unique(h.g[tr]))
    scores = {}
    for ci, cfg in enumerate(CONFIGS):
        if cfg is None:
            scores[(ci, 1.0)] = float(np.abs(r).mean() / s)   # no correction: held-out error is |residual|
            continue
        if ng < 3:
            scores[(ci, 1.0)] = 1e9
            continue
        kw = dict(cfg)
        kw.setdefault("epochs", INNER_EP)
        e = {sh: [] for sh in SHRINKS}
        for a, b in GroupKFold(min(3, ng)).split(tr, groups=h.g[tr]):
            Ha, Xa, Aa = sub(tr[a]); Hb, Xb, Ab = sub(tr[b])
            pr = MP.fit_predict(Ha, Xa, (h.y[tr[a]] - h.base[tr[a]]) / s, Hb, Xb, Aa, Ab, seeds=(0,), **kw)
            for sh in SHRINKS:
                e[sh].append(np.abs(sh * pr - (h.y[tr[b]] - h.base[tr[b]]) / s))
        for sh in SHRINKS:
            scores[(ci, sh)] = float(np.concatenate(e[sh]).mean())
    best, bsh = min(scores, key=scores.get)
    rec.append(best)
    if CONFIGS[best] is None or ng < 3:
        return h.base[te].copy()
    Ha, Xa, Aa = sub(tr); Hb, Xb, Ab = sub(te)
    kw = dict(CONFIGS[best])
    kw.setdefault("epochs", EPOCHS)
    pr = MP.fit_predict(Ha, Xa, r / s, Hb, Xb, Aa, Ab, seeds=SEEDS, **kw)
    return h.base[te] + bsh * pr * s


def fit_all(h, tr, te, kind, pure, rec, include_te=False):
    out = {}
    if kind == "octane":
        out["lin_molar"] = h.base[te].copy()
        out["corrected_pure"] = corrected_pure(h, tr, te, pure, include_te)
        out["desc_gbm"] = gbm_residual(h, tr, te, False)
        out["molpool"] = candidate(h, tr, te, False, rec)
    else:
        out["physics"] = h.base[te].copy()
        res = h.y[tr] - h.base[tr]
        tb = pd.Series(res).groupby(np.round(h.T[tr])).mean()
        out["t_bias"] = h.base[te] + np.array([tb.get(round(t), res.mean()) for t in h.T[te]])
        out["desc_gbm"] = gbm_residual(h, tr, te, True)
        out["molpool"] = candidate(h, tr, te, True, rec)
    return out


def run(stage):
    df = pd.read_csv(OUT / "rows.csv")
    enc = MP.load_encoder(OUT / ("encoder_quick.pt" if args.quick else "encoder.pt"))
    rows, cfgs = [], {}
    t0 = time.time()
    for (head, prop), sub in df.groupby(["head", "prop"]):
        kind = "boehm" if head == "boehm" else "octane"
        pure = pure_values(head, prop)
        dev = sub[sub.role == "dev"]
        h_dev = Head(dev, enc)
        if stage == "dev":
            n_splits = min(5, len(np.unique(h_dev.g)))
            preds = {k: np.full(len(dev), np.nan) for k in COMPARATORS[kind] + ["molpool"]}
            rec = []
            for tr, te in GroupKFold(n_splits).split(np.arange(len(dev)), groups=h_dev.g):
                o = fit_all(h_dev, tr, te, kind, pure, rec)
                for k, v in o.items():
                    preds[k][te] = v
            out = dev.reset_index(drop=True).copy()
        else:
            lk = sub[sub.role == "locked"].reset_index(drop=True)
            full = pd.concat([dev, lk]).reset_index(drop=True)
            hf = Head(full, enc)
            tr, te = np.arange(len(dev)), np.arange(len(dev), len(full))
            rec = []
            preds = fit_all(hf, tr, te, kind, pure, rec)
            out = lk.copy()
        for k, v in preds.items():
            out[f"pred_{k}"] = v
        rows.append(out)
        cfgs[f"{head}_{prop}"] = [str(CONFIGS[i]) for i in rec]
        print(stage, head, prop, f"{time.time() - t0:.0f}s", flush=True)
    res = pd.concat(rows, ignore_index=True)
    res.to_csv(OUT / (("dev_oof" if stage == "dev" else "locked_pred") + Q + ".csv"), index=False)
    json.dump(cfgs, open(OUT / (f"chosen_cfg_{stage}{Q}.json"), "w"), indent=1)


def boot(d, groups, rng, B=2000):
    ug = np.unique(groups)
    idx = {g: np.where(groups == g)[0] for g in ug}
    out = np.empty(B)
    for b in range(B):
        pick = rng.choice(len(ug), len(ug))
        out[b] = d[np.concatenate([idx[ug[i]] for i in pick])].mean()
    return np.percentile(out, [5, 95])


def summary():
    dev = pd.read_csv(OUT / f"dev_oof{Q}.csv")
    lk = None if args.v2 else pd.read_csv(OUT / "locked_pred.csv")
    rng = np.random.default_rng(0)
    rows, readings = [], {}
    for (head, prop), d in dev.groupby(["head", "prop"]):
        kind = "boehm" if head == "boehm" else "octane"
        l = None if args.v2 else lk[(lk["head"] == head) & (lk["prop"] == prop)]
        mae = {k: float(np.abs(d[f"pred_{k}"] - d.y).mean()) for k in COMPARATORS[kind] + ["molpool"]}
        best = min(COMPARATORS[kind], key=lambda k: mae[k])
        r = dict(head=head, prop=prop, best_comparator=best)
        gate = []
        for name, part in ((("dev", d),) if args.v2 else (("dev", d), ("locked", l))):
            ec = np.abs(part[f"pred_{best}"] - part.y).to_numpy()
            em = np.abs(part["pred_molpool"] - part.y).to_numpy()
            lo, hi = boot(ec - em, part.group.to_numpy(), rng)
            r.update({f"{name}_n": len(part), f"{name}_groups": part.group.nunique(), f"{name}_mae_best": ec.mean(), f"{name}_mae_molpool": em.mean(),
                      f"{name}_gain_lo": lo, f"{name}_gain_hi": hi})
            for k in COMPARATORS[kind]:
                r[f"{name}_mae_{k}"] = float(np.abs(part[f"pred_{k}"] - part.y).mean())
            mc = (part.n_comp >= 3).to_numpy()
            r[f"{name}_n_multi"] = int(mc.sum())
            r[f"{name}_mae_best_multi"] = float(ec[mc].mean()) if mc.any() else np.nan
            r[f"{name}_mae_molpool_multi"] = float(em[mc].mean()) if mc.any() else np.nan
            gate.append(lo > 0)
            gate.append(True if not mc.any() else bool(em[mc].mean() <= ec[mc].mean()))
        r["gate_pass"] = None if args.v2 else bool(all(gate))
        rows.append(r)
        readings[f"{head}_{prop}"] = r
    pd.DataFrame(rows).to_csv(OUT / f"summary{Q}.csv", index=False)
    json.dump(readings, open(OUT / f"readings{Q}.json", "w"), indent=1, default=float)
    show = ["head", "prop", "best_comparator", "dev_n", "dev_mae_best", "dev_mae_molpool", "dev_gain_lo", "dev_gain_hi"]
    if not args.v2:
        show += ["locked_n", "locked_mae_best", "locked_mae_molpool", "locked_gain_lo", "locked_gain_hi", "gate_pass"]
    print(pd.DataFrame(rows)[show].round(3).to_string())


if __name__ == "__main__":
    summary() if args.stage == "summary" else run(args.stage)
