"""v6 flash arm: oxygenate mixtures (results/v6/flash_oxy/PREDECLARATION.md).
  --stage dev      GroupKFold(5) over dev systems (hc+oxy jointly); arm B trains on hc+oxy, arm A on hc only  -> dev_oof.csv
  --stage locked   fit on all dev, score locked systems ONCE                                               -> locked_pred.csv
  --stage summary  paired system bootstrap, gate, non-inferiority, safety                                  -> readings.json
--quick is a smoke run (writes *_quick files, never touches the locked file)."""
import sys, json, time, argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GroupKFold
RDLogger.DisableLog("rdApp.*"); warnings.simplefilter("ignore")
R = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(R / "src"))
from fuelblend.ml import molpool as MP   # noqa: E402
from fuelblend.baselines.flash import predict_flash_le_chatelier   # noqa: E402
from fuelblend.baselines.flash_unifac import _prepare   # noqa: E402

OUT = R / "results/v6/flash_oxy"
ENC = R / "results/v6/mixture/encoder.pt"
CONFIGS = [None, dict(d=8, wd=1e-2), dict(d=8, wd=1e-1), dict(d=24, wd=1e-2), dict(d=24, wd=1e-1)]   # v1 configs, None = no correction
THETAS = [10, 20, 40, 80, 160, 1e4]   # blending-index scale (degC); 1e4 is linear molar
ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["dev", "locked", "summary"])
ap.add_argument("--quick", action="store_true")
args = ap.parse_args()
Q = "_quick" if args.quick else ""
EPOCHS, SEEDS, INNER_EP = (15, (0,), 15) if args.quick else (150, (0, 1, 2, 3, 4), 120)
COMP = ["lin_molar", "blend_law", "le_chatelier", "desc_gbm"]


def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    return np.array([Descriptors.MolWt(m) / 100, rdMolDescriptors.CalcNumAromaticRings(m), rdMolDescriptors.CalcFractionCSP3(m),
                     sum(a.GetAtomicNum() == 8 for a in m.GetAtoms()), rdMolDescriptors.CalcNumRings(m),
                     sum(a.GetAtomicNum() == 6 for a in m.GetAtoms()) / 10], float)


class Data:
    def __init__(self, df, enc):
        self.df = df.reset_index(drop=True)
        self.sm = [json.loads(s) for s in self.df.smiles]
        self.x = [np.array(json.loads(s), float) for s in self.df.x]
        self.pf = [np.array([np.nan if v is None else v for v in json.loads(s)], float) for s in self.df.pure_flash]
        allsmi = sorted({s for r in self.sm for s in r})
        E = dict(zip(allsmi, MP.embed(enc, allsmi)))
        Dsc = {s: descriptors(s) for s in allsmi}
        self.H, self.X = MP.pad_sets([np.stack([E[s] for s in r]) for r in self.sm], self.x)
        self.D = np.stack([(x[:, None] * np.stack([Dsc[s] for s in r])).sum(0) for r, x in zip(self.sm, self.x)])
        self.y = self.df.y.to_numpy(float)
        self.g = self.df.system_id.to_numpy()
        self.cls = self.df.chem_class.to_numpy()
        self.n = self.df.n_comp.to_numpy()
        self._b = {}

    def blend(self, theta):
        if theta in self._b:
            return self._b[theta]
        out = np.empty(len(self.y))
        for i, (x, t) in enumerate(zip(self.x, self.pf)):
            h = ~np.isnan(t)
            xx = x[h] / x[h].sum()
            out[i] = (xx * t[h]).sum() if theta >= 1e4 else -theta * np.log((xx * np.exp(-t[h] / theta)).sum())
        self._b[theta] = out
        return out

    def antoine(self):
        out = np.full(len(self.y), np.nan)
        for i, (s, x) in enumerate(zip(self.sm, self.x)):
            xn, ant, lfl, missing = _prepare(s, x, None)
            if missing:
                continue
            v, st = predict_flash_le_chatelier(list(xn), ant, list(lfl), lfl_temp=True)
            out[i] = v if st == "ok" else np.nan
        return out


def fit_theta(d, tr):
    return min(THETAS, key=lambda th: np.abs(d.blend(th)[tr] - d.y[tr]).mean())


def gbm_resid(d, base, tr, te):
    f = lambda i: np.hstack([d.D[i], base[i, None]])
    m = GradientBoostingRegressor(n_estimators=45, max_depth=3, learning_rate=0.08, random_state=42)
    m.fit(f(tr), d.y[tr] - base[tr])
    return base[te] + m.predict(f(te))


def candidate(d, base, tr, te):
    r = d.y[tr] - base[tr]
    s = float(r.std()) if r.std() > 1e-9 else 1.0
    ng = len(np.unique(d.g[tr]))
    sub = lambda idx: (d.H[idx], d.X[idx], None)
    scores = {}
    for ci, cfg in enumerate(CONFIGS):
        if cfg is None:
            scores[ci] = float(np.abs(r).mean() / s)
            continue
        if ng < 3:
            scores[ci] = 1e9
            continue
        e = []
        for a, b in GroupKFold(min(3, ng)).split(tr, groups=d.g[tr]):
            Ha, Xa, Aa = sub(tr[a]); Hb, Xb, Ab = sub(tr[b])
            pr = MP.fit_predict(Ha, Xa, (d.y[tr[a]] - base[tr[a]]) / s, Hb, Xb, Aa, Ab, seeds=(0,), epochs=INNER_EP, **cfg)
            e.append(np.abs(pr - (d.y[tr[b]] - base[tr[b]]) / s))
        scores[ci] = float(np.concatenate(e).mean())
    best = min(scores, key=scores.get)
    if CONFIGS[best] is None or ng < 3:
        return base[te].copy(), best
    Ha, Xa, Aa = sub(tr); Hb, Xb, Ab = sub(te)
    pr = MP.fit_predict(Ha, Xa, r / s, Hb, Xb, Aa, Ab, seeds=SEEDS, epochs=EPOCHS, **CONFIGS[best])
    return base[te] + pr * s, best


def fit_arm(d, tr, te, base_kind, ant):
    """One arm in one fold: comparators and MolPool, all trained on rows tr only; returns predictions for te."""
    th = fit_theta(d, tr)
    lin, bl = d.blend(1e4), d.blend(th)
    base = lin if base_kind == "lin_molar" else bl
    cand, ci = candidate(d, base, tr, te)
    return dict(lin_molar=lin[te], blend_law=bl[te], le_chatelier=ant[te], desc_gbm=gbm_resid(d, base, tr, te), molpool=cand,
                _theta=th, _cfg=ci)


def run(stage):
    df = pd.read_csv(OUT / "rows.csv")
    enc = MP.load_encoder(ENC)
    dev = df[df.role == "dev"].reset_index(drop=True)
    t0 = time.time()
    if stage == "dev":
        d = Data(dev, enc)
        ant = d.antoine()
        hc_i = np.where(d.cls == "hc")[0]
        folds = list(GroupKFold(5).split(np.arange(len(d.y)), groups=d.g))
        # pass 1: base selection by dev OOF MAE of the two blending baselines (arm B: all dev rows; arm A: hc rows)
        oof_bl = np.full(len(d.y), np.nan)
        for tr, te in folds:
            oof_bl[te] = d.blend(fit_theta(d, tr))[te]
        lin = d.blend(1e4)
        kinds = {}
        for nm, m in (("B", np.ones(len(d.y), bool)), ("A", d.cls == "hc")):
            kinds[nm] = min(("lin_molar", "blend_law"), key=lambda k: np.abs((lin if k == "lin_molar" else oof_bl)[m] - d.y[m]).mean())
        preds = {a: {k: np.full(len(d.y), np.nan) for k in COMP + ["molpool"]} for a in "AB"}
        thetas, cfgs = [], []
        for tr, te in folds:
            for arm in "BA":
                trr = tr if arm == "B" else np.intersect1d(tr, hc_i)
                tee = te if arm == "B" else np.intersect1d(te, hc_i)
                if len(tee) == 0:
                    continue
                o = fit_arm(d, trr, tee, kinds[arm], ant)
                for k in COMP + ["molpool"]:
                    preds[arm][k][tee] = o[k]
                thetas.append((arm, o["_theta"])); cfgs.append((arm, int(o["_cfg"])))
            print("fold", f"{time.time() - t0:.0f}s", flush=True)
        out = dev.copy()
        for a in "AB":
            for k, v in preds[a].items():
                out[f"pred{a}_{k}"] = v
        out.to_csv(OUT / f"dev_oof{Q}.csv", index=False)
        json.dump(dict(base_kind=kinds, thetas=thetas, cfgs=cfgs, antoine_coverage=float(np.isfinite(ant).mean()),
                       antoine_coverage_oxy=float(np.isfinite(ant[d.cls == "oxy"]).mean())), open(OUT / f"chosen_dev{Q}.json", "w"), indent=1)
    else:
        assert not args.quick, "locked is scored once, never quick"
        assert not (OUT / "locked_pred.csv").exists(), "locked set already scored"
        lk = df[df.role == "locked"].reset_index(drop=True)
        full = pd.concat([dev, lk]).reset_index(drop=True)
        d = Data(full, enc)
        ant = d.antoine()
        nd = len(dev)
        tr_all, te = np.arange(nd), np.arange(nd, len(full))
        kinds = json.load(open(OUT / "chosen_dev.json"))["base_kind"]
        out = lk.copy()
        for arm in "BA":
            trr = tr_all if arm == "B" else tr_all[d.cls[tr_all] == "hc"]
            tee = te if arm == "B" else te[d.cls[te] == "hc"]
            o = fit_arm(d, trr, tee, kinds[arm], ant)
            for k in COMP + ["molpool"]:
                col = np.full(len(lk), np.nan)
                col[tee - nd] = o[k]
                out[f"pred{arm}_{k}"] = col
            print("locked arm", arm, f"{time.time() - t0:.0f}s", o["_theta"], o["_cfg"], flush=True)
        out.to_csv(OUT / "locked_pred.csv", index=False)


def boot(v, groups, rng, B=2000):
    ug = np.unique(groups)
    idx = {g: np.where(groups == g)[0] for g in ug}
    o = np.empty(B)
    for b in range(B):
        o[b] = v[np.concatenate([idx[ug[i]] for i in rng.choice(len(ug), len(ug))])].mean()
    return [float(x) for x in np.percentile(o, [5, 95])]


def safety(e):
    return dict(mean_signed=float(e.mean()), mean_over=float(np.clip(e, 0, None).mean()), p90_signed=float(np.percentile(e, 90)))


def summary():
    dev = pd.read_csv(OUT / f"dev_oof{Q}.csv")
    lk = None if args.quick else pd.read_csv(OUT / "locked_pred.csv")
    rng = np.random.default_rng(0)
    ant_cov = float(dev[dev.chem_class == "oxy"].predB_le_chatelier.notna().mean())
    cands = [k for k in COMP if k != "le_chatelier" or ant_cov >= 0.9]
    do = dev[dev.chem_class == "oxy"]
    mae = {k: float(np.abs(do[f"predB_{k}"] - do.y).mean()) for k in cands}
    best = min(cands, key=lambda k: mae[k])
    res = dict(gate_comparator=best, dev_oxy_mae_comparators=mae, antoine_coverage_oxy_dev=ant_cov)
    gate = []
    parts = (("dev", dev), ("locked", lk)) if lk is not None else (("dev", dev),)
    for name, part in parts:
        o = part[part.chem_class == "oxy"]
        ec = (o[f"predB_{best}"] - o.y).to_numpy()
        em = (o.predB_molpool - o.y).to_numpy()
        lo, hi = boot(np.abs(ec) - np.abs(em), o.system_id.to_numpy(), rng)
        m3 = (o.n_comp >= 3).to_numpy()
        ns3 = int(o[m3].system_id.nunique())
        r = dict(n=len(o), systems=int(o.system_id.nunique()), mae_comparator=float(np.abs(ec).mean()), mae_molpool=float(np.abs(em).mean()),
                 gain_ci90=[lo, hi], n_multi_systems=ns3,
                 mae_multi_comparator=float(np.abs(ec[m3]).mean()) if m3.any() else None,
                 mae_multi_molpool=float(np.abs(em[m3]).mean()) if m3.any() else None,
                 safety_comparator=safety(ec), safety_molpool=safety(em),
                 all_comparators_mae={k: float(np.abs(o[f"predB_{k}"] - o.y).mean()) for k in COMP if o[f"predB_{k}"].notna().any()})
        gate.append(lo > 0)
        gate.append(True if ns3 < 5 else bool(np.abs(em[m3]).mean() <= np.abs(ec[m3]).mean()))
        r["veto_overprediction"] = bool(r["safety_molpool"]["p90_signed"] - r["safety_comparator"]["p90_signed"] > 1.0)
        h = part[part.chem_class == "hc"]
        eb, ea = (h.predB_molpool - h.y).to_numpy(), (h.predA_molpool - h.y).to_numpy()
        lo2, hi2 = boot(np.abs(eb) - np.abs(ea), h.system_id.to_numpy(), rng)
        r["hc"] = dict(n=len(h), systems=int(h.system_id.nunique()), mae_armA=float(np.abs(ea).mean()), mae_armB=float(np.abs(eb).mean()),
                       diff_B_minus_A_ci90=[lo2, hi2], noninferior=bool(hi2 <= 1.0), safety_armA=safety(ea), safety_armB=safety(eb),
                       mae_lin=float(np.abs(h.predB_lin_molar - h.y).mean()), mae_blend=float(np.abs(h.predB_blend_law - h.y).mean()))
        res[name] = r
    res["primary_gate_pass"] = None if lk is None else bool(all(gate))
    json.dump(res, open(OUT / f"readings{Q}.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    summary() if args.stage == "summary" else run(args.stage)
