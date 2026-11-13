"""v5 WS4 delta_ws4_mixture_ml: pairwise DeepSets / MolSets-style GNN / multi-task set transformer versus the LightGBM pair residuals.
Declared text: results/v5/BASE/delta_ws4_mixture_ml/README.md (gates fixed before any fit).

Usage:  python scripts/v5/base/delta_ws4_mixture_ml.py --stage dev      (dev out-of-fold only; locked systems NOT scored)
        python scripts/v5/base/delta_ws4_mixture_ml.py --stage locked   (final 5-seed ensembles; locked scored ONCE; reads the dev OOF file)
Flags for the GPU pack: --device cuda --epochs 300 --out <dir>.  The declared CPU reading uses the defaults.
The baseline errors are READ from the committed delta_ml_viscosity / delta_ml_density rows_scored.csv; nothing under models/v5 is touched.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT_DEFAULT = ROOT / "results" / "v5" / "BASE" / "delta_ws4_mixture_ml"
BASE_V = ROOT / "results" / "v5" / "BASE" / "delta_ml_viscosity" / "rows_scored.csv"
BASE_D = ROOT / "results" / "v5" / "BASE" / "delta_ml_density" / "rows_scored.csv"
CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]
N_CONT = 8                   # continuous feature columns standardised (6 static + ln_eta + v_molar)
K = 3
T_WARM = 283.15
WINS = {"viscosity": 0.5, "density": 2.0}
XMIN = 0.01
SEED = 42
CANDS = ("pds", "gnn", "st")


def load_by_path(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


# ------------------------------------------------------------------------------------------------ data
def _onehot(c):
    return [1.0 if c == k else 0.0 for k in CLASSES]


def viscosity_data(mv, ctx) -> dict:
    R = mv.regenerate_rows(ctx)
    B = pd.read_csv(BASE_V)
    assert len(B) == len(R) and (B.system_id.to_numpy() == R.system_id.to_numpy()).all(), "viscosity rows differ from baseline rows_scored"
    assert np.abs(B.T_K.to_numpy() - R.T_K.to_numpy()).max() < 1e-3 and (B.n_comp.to_numpy() == R.n_comp.to_numpy()).all()
    assert (B.role.to_numpy() == R.role.to_numpy()).all()
    phys = R.ln_pred.to_numpy() - R.ln_meas.to_numpy()
    assert np.abs(phys - B.resid_refutas.to_numpy()).max() < 1e-5, "Refutas residual differs from baseline file"
    n = len(R)
    F = np.zeros((n, K, 14), np.float32)
    X = np.zeros((n, K), np.float32)
    M = np.zeros((n, K), bool)
    mol = np.zeros((n, K), np.int64)
    uniq = sorted({i for ii in R.ii for i in ii})
    pos = {i: k for k, i in enumerate(uniq)}
    rr, ss, ci, TT = [], [], [], []
    for r, (ii, x, T) in enumerate(zip(R.ii, R.x, R.T_K)):
        assert len(ii) <= K
        for s, (i, xi) in enumerate(zip(ii, x)):
            rr.append(r); ss.append(s); ci.append(i); TT.append(T)
            X[r, s], M[r, s], mol[r, s] = xi, True, pos[i]
    rr, ss, ci, TT = map(np.array, (rr, ss, ci, TT))
    le, rho = ctx.prov.lookup(ci, TT)
    vm = ctx.reg.mw[ci] / rho * 1000.0
    F[rr, ss, :6] = ctx.static[ci]
    F[rr, ss, 6] = le
    F[rr, ss, 7] = vm
    F[rr, ss, 8:13] = np.array([_onehot(ctx.cls[i]) for i in ci])
    assert np.isfinite(F).all(), "non-finite viscosity component feature"
    X = X / X.sum(1, keepdims=True)
    from fuelblend.ml.set_models import GraphBank
    bank = GraphBank.from_smiles([ctx.reg.smiles[i] for i in uniq])
    xx = np.array([x[0] * x[1] if len(x) == 2 else np.nan for x in R.x])
    r_target = np.clip(R.ln_meas.to_numpy() - R.ln_pred.to_numpy(), -WINS["viscosity"], WINS["viscosity"])
    return dict(name="viscosity", R=R, F=F, X=X, M=M, invT=1000.0 / R.T_K.to_numpy(), mol=mol, bank=bank,
                sys=R.system_id.to_numpy(), role=R.role.to_numpy(), n_comp=R.n_comp.to_numpy(), t_bin=R.t_bin.to_numpy(),
                T=R.T_K.to_numpy(), eligible=(R.n_comp.to_numpy() == 2) & (xx >= XMIN), y=r_target,
                phys_err=phys, base_err=B.resid_ml.to_numpy(), err_fn=lambda c, phys=phys: phys + c,
                primary=R.T_K.to_numpy() >= T_WARM, unit="ln")


def density_data(md, ctx) -> dict:
    from rdkit import Chem
    R, desc = md.regenerate()
    B = pd.read_csv(BASE_D)
    assert len(B) == len(R) and (B.system_id.to_numpy() == R.system_id.to_numpy()).all(), "density rows differ from baseline rows_scored"
    assert np.abs(B.T_K.to_numpy() - R.T_K.to_numpy()).max() < 1e-3 and (B.n_comp.to_numpy() == R.n_comp.to_numpy()).all()
    assert (B.role.to_numpy() == R.role.to_numpy()).all()
    phys = R.pred.to_numpy() - R.meas.to_numpy()
    assert np.abs(phys - B.err_ideal.to_numpy()).max() < 1e-3, "ideal-mixing error differs from baseline file"
    lab = sys.modules["eval_lab_known_composition"]
    thermo = pd.read_csv(md.THERMOML, low_memory=False)
    inchi_of = {}
    for s in thermo[thermo.prop.str.contains("Mass density") & (thermo.n_comp >= 2)].inchis.unique():
        for inchi in json.loads(s):
            info = lab._inchi_info(inchi)
            if info and info[0] not in inchi_of:
                inchi_of[info[0]] = inchi
    n = len(R)
    F = np.zeros((n, K, 14), np.float32)
    X = np.zeros((n, K), np.float32)
    M = np.zeros((n, K), bool)
    mol = np.zeros((n, K), np.int64)
    uniq = sorted({k for iks in R.iks for k in iks})
    pos = {k: j for j, k in enumerate(uniq)}
    reg_i = {k: ctx.reg.ik2i.get(k) for k in uniq}
    for r, (iks, x, rho, T) in enumerate(zip(R.iks, R.x, R.rho_i, R.T_K)):
        assert len(iks) <= K
        for s, (k, xi, rh) in enumerate(zip(iks, x, rho)):
            st, cl = desc[k]
            F[r, s, :6] = st
            F[r, s, 7] = st[1] / rh * 1000.0
            F[r, s, 8:13] = _onehot(cl)
            i = reg_i[k]
            le = np.nan
            if i is not None and ctx.usable(i):
                le = float(np.asarray(ctx.prov.lookup(np.array([i]), np.array([T]))[0])[0])
            if np.isfinite(le):
                F[r, s, 6] = le
            else:
                F[r, s, 13] = 1.0
            X[r, s], M[r, s], mol[r, s] = xi, True, pos[k]
    assert np.isfinite(F).all()
    from fuelblend.ml.set_models import GraphBank
    bank = GraphBank.from_mols([Chem.MolFromInchi(inchi_of[k]) if k in inchi_of else None for k in uniq])
    missing_graph = [k for k in uniq if k not in inchi_of]
    xx = np.array([x[0] * x[1] if len(x) == 2 else np.nan for x in R.x])
    ve = np.clip(R.v_meas.to_numpy() - R.v_ideal.to_numpy(), -WINS["density"], WINS["density"])
    mm, vid, meas = R.mm.to_numpy(), R.v_ideal.to_numpy(), R.meas.to_numpy()
    return dict(name="density", R=R, F=F, X=X, M=M, invT=1000.0 / R.T_K.to_numpy(), mol=mol, bank=bank,
                sys=R.system_id.to_numpy(), role=R.role.to_numpy(), n_comp=R.n_comp.to_numpy(), t_bin=R.t_bin.to_numpy(),
                T=R.T_K.to_numpy(), eligible=(R.n_comp.to_numpy() == 2) & (xx >= XMIN), y=ve,
                phys_err=phys, base_err=B.err_ml.to_numpy(), err_fn=lambda c: mm / (vid + c) * 1000.0 - meas,
                primary=np.ones(n, bool), unit="kg/m3", missing_graph=missing_graph,
                ln_eta_missing_share=float(F[M][:, 13].mean()))


# ------------------------------------------------------------------------------------------------ fitting
class Scaler:
    def __init__(self, D_list, rows_list):
        live = np.concatenate([D["F"][r][D["M"][r]][:, :N_CONT] for D, r in zip(D_list, rows_list)])
        self.mu, self.sd = live.mean(0), live.std(0) + 1e-6
        t = np.concatenate([D["invT"][r] for D, r in zip(D_list, rows_list)])
        self.tmu, self.tsd = t.mean(), t.std() + 1e-6

    def batch(self, D, rows, prop):
        from fuelblend.ml.set_models import SetBatch
        F = D["F"][rows].copy()
        F[..., :N_CONT] = (F[..., :N_CONT] - self.mu) / self.sd
        F[~D["M"][rows]] = 0.0
        return SetBatch.from_numpy(F, D["X"][rows], D["M"][rows], (D["invT"][rows] - self.tmu) / self.tsd, D["mol"][rows],
                                   np.full(len(rows), prop))


def make_factory(cand, D, args):
    from fuelblend.ml import set_models as sm
    if cand == "pds":
        return lambda: sm.PairDeepSets(width=args.width)
    if cand == "gnn":
        return lambda: sm.MolSetsGNN(D["bank"], width=args.width)
    return lambda: sm.SetTransformer(width=args.width, n_props=2)


def fit_predict(cand, DS, train_rows, pred_rows, args, tag):
    """DS: list of property dicts (1 for pds/gnn, 2 for st). train_rows / pred_rows: per-property index arrays.
    Returns per-property (mean, std) corrections in raw target units for pred_rows."""
    from fuelblend.ml import set_models as sm
    import torch
    sc = Scaler(DS, train_rows)
    bs, ys, sysl, props, ystd = [], [], [], [], []
    for p, (D, tr) in enumerate(zip(DS, train_rows)):
        sd = float(np.std(D["y"][tr])) + 1e-9
        ystd.append(sd)
        bs.append(sc.batch(D, tr, p))
        ys.append(D["y"][tr] / sd)
        sysl.append(D["sys"][tr])
        props.append(np.full(len(tr), p))
    cat = lambda name: torch.cat([getattr(b, name) for b in bs])  # noqa: E731
    data = sm.SetBatch(cat("F"), cat("X"), cat("M"), cat("invT"), cat("mol"), cat("prop"))
    y = np.concatenate(ys)
    w = sm.system_macro_weights(np.concatenate(sysl), np.concatenate(props))
    t0 = time.time()
    ens = sm.train_ensemble(make_factory(cand, DS[0], args), data, y, w, seeds=list(range(args.seeds)), epochs=args.epochs,
                            batch_size=args.batch, lr=args.lr, wd=args.wd, device=args.device)
    out = []
    for p, (D, pr) in enumerate(zip(DS, pred_rows)):
        if len(pr) == 0:
            out.append((np.zeros(0), np.zeros(0)))
            continue
        mu, sd = sm.predict_ensemble(ens, sc.batch(D, pr, p), device=args.device)
        out.append((mu * ystd[p], sd * ystd[p]))
    log(f"{tag} {cand}: train rows {[len(t) for t in train_rows]}, {time.time() - t0:.0f}s")
    return out


def dev_oof(cand, Dp, Dq, args):
    """Out-of-fold corrections for property Dp (all its dev rows) with the baseline's exact folds.
    For the set transformer, Dq rows enter training except those of systems in Dp's test fold."""
    from fuelblend.ml.set_models import fold_test_systems
    dev_rows = np.where(Dp["role"] == "dev")[0]
    folds = fold_test_systems(Dp["sys"], dev_rows)
    mu, sd = np.full(len(Dp["sys"]), np.nan), np.full(len(Dp["sys"]), np.nan)
    for f, test_sys in enumerate(folds):
        test_sys = set(test_sys)
        trp = np.where((Dp["role"] == "dev") & Dp["eligible"] & ~np.isin(Dp["sys"], list(test_sys)))[0]
        tep = np.where((Dp["role"] == "dev") & np.isin(Dp["sys"], list(test_sys)))[0]
        if cand == "st":
            trq = np.where((Dq["role"] == "dev") & Dq["eligible"] & ~np.isin(Dq["sys"], list(test_sys)))[0]
            DS, tr, pr = ([Dp, Dq], [trp, trq], [tep, np.zeros(0, int)]) if Dp["name"] == "viscosity" else \
                ([Dq, Dp], [trq, trp], [np.zeros(0, int), tep])
            res = fit_predict(cand, DS, tr, pr, args, f"{Dp['name']} fold {f}")
            m_, s_ = res[0] if Dp["name"] == "viscosity" else res[1]
        else:
            (m_, s_), = fit_predict(cand, [Dp], [trp], [tep], args, f"{Dp['name']} fold {f}")
        mu[tep], sd[tep] = m_, s_
    return mu, sd, folds


# ------------------------------------------------------------------------------------------------ readings
def boot(mv, D, ea, eb, mask):
    R = D["R"]
    return mv.paired_boot(mv.macro(R, ea, mask), mv.macro(R, eb, mask), np.random.default_rng(SEED))


def readings_for(mv, D, err_c, mask_role, label):
    from fuelblend.ml.set_models import gate_pass
    pm = mask_role & D["primary"]
    vs_pair = boot(mv, D, D["base_err"], err_c, pm)
    vs_phys = boot(mv, D, D["phys_err"], err_c, pm)
    pair_vs_phys = boot(mv, D, D["phys_err"], D["base_err"], pm)
    return {"mask": label, "vs_pair": vs_pair, "vs_physics": vs_phys, "pair_vs_physics_check": pair_vs_phys, "pass": gate_pass(vs_pair)}


def ternary_reading(D, E, mask_role):
    tm = (D["n_comp"] == 3) & mask_role
    out = dict(rows=int(tm.sum()), systems=int(len(set(D["sys"][tm]))))
    for k, e in E.items():
        out[f"mae_{k}"] = float(np.abs(e[tm]).mean()) if tm.any() else None
    return out


def cold_reading(mv, D, err_c, mask_role):
    cm = mask_role & (D["T"] < T_WARM)
    if not cm.any():
        return None
    return dict(rows=int(cm.sum()), systems=int(len(set(D["sys"][cm]))), vs_pair=boot(mv, D, D["base_err"], err_c, cm),
                vs_physics=boot(mv, D, D["phys_err"], err_c, cm),
                under_pred_rate=dict(physics=float((D["phys_err"][cm] < 0).mean()), pair=float((D["base_err"][cm] < 0).mean()),
                                     candidate=float((err_c[cm] < 0).mean())))


def uq_reading(mv, D, err_c, std_c, dev, lock):
    from scipy.stats import spearmanr
    out = {"spearman_std_vs_abs_err_dev": float(spearmanr(std_c[dev], np.abs(err_c[dev])).statistic)}
    sysarr = D["sys"]
    if D["name"] == "viscosity":
        for b in ("T_lt_263K", "T_263_283K", "T_ge_283K"):
            msk = dev & (D["t_bin"] == b)
            d = pd.DataFrame({"s": sysarr[msk], "e": np.abs(err_c[msk])}).groupby("s").e.max()
            v, ng, k = mv.conformal(d.to_numpy()) if len(d) >= 9 else (float("inf"), len(d), None)
            entry = dict(margin_ln=v, n_dev_systems=int(ng))
            if lock is not None:
                ml = lock & (D["t_bin"] == b)
                dl = pd.DataFrame({"s": sysarr[ml], "e": np.abs(err_c[ml])}).groupby("s").e.max()
                entry.update(locked_systems=int(len(dl)), locked_covered=int((dl <= v).sum()) if np.isfinite(v) else None)
            out[b] = entry
    else:
        from fuelblend.uq import conformal as cf
        v = cf.two_sided_margin(err_c[dev], sysarr[dev], 0.10)
        out["dev_oof_margin_kg_m3"] = v
        if lock is not None:
            out["locked_coverage_by_dev_margin"] = list(cf.coverage(err_c[lock], sysarr[lock], v, two_sided=True))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=("dev", "locked"), required=True)
    ap.add_argument("--cands", default=",".join(CANDS))
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--width", type=int, default=64)
    ap.add_argument("--batch", type=int, default=512)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--wd", type=float, default=1e-4)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    args = ap.parse_args()
    warnings.simplefilter("ignore")
    import torch
    torch.set_num_threads(args.threads)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cands = [c for c in args.cands.split(",") if c]
    config = dict(epochs=args.epochs, seeds=args.seeds, width=args.width, batch=args.batch, lr=args.lr, wd=args.wd, threads=args.threads,
                  device=args.device, torch=torch.__version__, cands=cands)
    log(f"config {config}")

    mv = load_by_path("delta_ml_viscosity", "scripts/v5/base/delta_ml_viscosity.py")
    md = load_by_path("delta_ml_density", "scripts/v5/base/delta_ml_density.py")
    from fuelblend.ml import viscosity_pair as vp
    ctx = vp.context()
    DV = viscosity_data(mv, ctx)
    log(f"viscosity rows {len(DV['sys'])}, systems {len(set(DV['sys']))}, molecules {DV['bank'].n_mol}; aligned with baseline rows_scored")
    DD = density_data(md, ctx)
    log(f"density rows {len(DD['sys'])}, systems {len(set(DD['sys']))}, molecules {DD['bank'].n_mol} (no graph: {len(DD['missing_graph'])}),"
        f" ln_eta missing share {DD['ln_eta_missing_share']:.3f}; aligned with baseline rows_scored")
    roles_v = dict(zip(DV["sys"], DV["role"]))
    clash = [s for s, r in zip(DD["sys"], DD["role"]) if s in roles_v and roles_v[s] != r]
    assert not clash, f"role disagreement between tables: {sorted(set(clash))[:5]}"
    shared = len(set(DV["sys"]) & set(DD["sys"]))
    props = {"viscosity": DV, "density": DD}
    oof_path = out / "oof_corrections.csv.gz"

    if args.stage == "dev":
        cols, folds_all = {}, {}
        for pname, Dp in props.items():
            Dq = DD if pname == "viscosity" else DV
            for c in cands:
                mu, sd, folds = dev_oof(c, Dp, Dq, args)
                cols[(pname, c)] = (mu, sd)
                folds_all[pname] = folds
        for pname, folds in folds_all.items():
            (out / f"folds_{pname}.json").write_text(json.dumps({f"fold_{i}": f for i, f in enumerate(folds)}, indent=0))
        frames = []
        for pname, Dp in props.items():
            fr = pd.DataFrame({"property": pname, "row": np.arange(len(Dp["sys"])), "system_id": Dp["sys"], "role": Dp["role"]})
            for c in cands:
                if (pname, c) in cols:
                    fr[f"corr_{c}"], fr[f"std_{c}"] = cols[(pname, c)]
            frames.append(fr)
        pd.concat(frames).to_csv(oof_path, index=False, float_format="%.7g")
        readings = dict(config=config, shared_systems_between_tables=shared, density_ln_eta_missing_share=DD["ln_eta_missing_share"],
                        density_molecules_without_graph=len(DD["missing_graph"]))
        for pname, Dp in props.items():
            dev = Dp["role"] == "dev"
            for c in cands:
                mu, sd = cols[(pname, c)]
                e = np.where(dev, Dp["err_fn"](np.nan_to_num(mu)), np.nan)
                r = {"R-WS4-1": readings_for(mv, Dp, e, dev, "dev OOF, " + ("T >= 283.15 K" if pname == "viscosity" else "all T"))}
                r["ternary_dev_oof"] = ternary_reading(Dp, {"physics": Dp["phys_err"], "pair": Dp["base_err"], c: e}, dev)
                if pname == "viscosity":
                    r["R-WS4-4_cold_dev"] = cold_reading(mv, Dp, e, dev)
                r["R-WS4-5_dev"] = uq_reading(mv, Dp, e, sd, dev, None)
                readings[f"{pname}/{c}"] = r
                log(f"R-WS4-1 {pname}/{c}: pair {r['R-WS4-1']['vs_pair']['mae_a']:.4g} -> cand {r['R-WS4-1']['vs_pair']['mae_b']:.4g}, "
                    f"gain [{r['R-WS4-1']['vs_pair']['lo']:.4g}, {r['R-WS4-1']['vs_pair']['hi']:.4g}] pass {r['R-WS4-1']['pass']}")
        (out / "readings_dev.json").write_text(json.dumps(readings, indent=1, default=float))
        log("dev stage done; locked systems not scored")
        return

    # ---------------- locked stage: final ensembles on all dev training rows, locked scored once
    oof = pd.read_csv(oof_path)
    dev_readings = json.loads((out / "readings_dev.json").read_text())
    readings = dict(config=config, dev_config=dev_readings["config"], shared_systems_between_tables=shared)
    final = {}
    trv = np.where((DV["role"] == "dev") & DV["eligible"])[0]
    trd = np.where((DD["role"] == "dev") & DD["eligible"])[0]
    lv, ld = np.where(DV["role"] == "locked")[0], np.where(DD["role"] == "locked")[0]
    for c in cands:
        if c == "st":
            (a, b) = fit_predict(c, [DV, DD], [trv, trd], [lv, ld], args, "final")
            final[("viscosity", c)], final[("density", c)] = a, b
        else:
            final[("viscosity", c)], = fit_predict(c, [DV], [trv], [lv], args, "final")
            final[("density", c)], = fit_predict(c, [DD], [trd], [ld], args, "final")
    rows_out = []
    for pname, Dp in props.items():
        dev, lock = Dp["role"] == "dev", Dp["role"] == "locked"
        o = oof[oof.property == pname].sort_values("row")
        assert len(o) == len(Dp["sys"]) and (o.system_id.to_numpy() == Dp["sys"]).all()
        lr = np.where(lock)[0]
        E = {"physics": Dp["phys_err"], "pair": Dp["base_err"]}
        S = {}
        for c in cands:
            corr = o[f"corr_{c}"].to_numpy().copy()
            sd = o[f"std_{c}"].to_numpy().copy()
            corr[lr], sd[lr] = final[(pname, c)]
            assert np.isfinite(corr[dev | lock]).all()
            E[c], S[c] = Dp["err_fn"](corr), sd
        for c in cands:
            r = {"R-WS4-1": dev_readings[f"{pname}/{c}"]["R-WS4-1"],
                 "R-WS4-2": readings_for(mv, Dp, E[c], lock, "locked, " + ("T >= 283.15 K" if pname == "viscosity" else "all T"))}
            tern = ternary_reading(Dp, {k: E[k] for k in ("physics", "pair", c)}, dev | lock)
            tern["dev_oof"] = ternary_reading(Dp, {k: E[k] for k in ("physics", "pair", c)}, dev)
            tern["locked"] = ternary_reading(Dp, {k: E[k] for k in ("physics", "pair", c)}, lock)
            tern["pass"] = bool(tern["rows"] > 0 and tern[f"mae_{c}"] <= tern["mae_pair"])
            r["R-WS4-3"] = tern
            r["gate_pass"] = bool(r["R-WS4-1"]["pass"] and r["R-WS4-2"]["pass"] and tern["pass"])
            if pname == "viscosity":
                r["R-WS4-4_cold_dev_plus_locked"] = cold_reading(mv, Dp, E[c], dev | lock)
            r["R-WS4-5"] = uq_reading(mv, Dp, E[c], S[c], dev, lock)
            readings[f"{pname}/{c}"] = r
            log(f"{pname}/{c}: R1 {r['R-WS4-1']['pass']} R2 {r['R-WS4-2']['pass']} "
                f"(locked pair {r['R-WS4-2']['vs_pair']['mae_a']:.4g} -> {r['R-WS4-2']['vs_pair']['mae_b']:.4g}, "
                f"[{r['R-WS4-2']['vs_pair']['lo']:.4g}, {r['R-WS4-2']['vs_pair']['hi']:.4g}]) R3 {tern['pass']} "
                f"(tern pair {tern['mae_pair']:.4g} cand {tern[f'mae_{c}']:.4g}) GATE {r['gate_pass']}")
        fr = pd.DataFrame({"property": pname, "system_id": Dp["sys"], "T_K": Dp["T"], "n_comp": Dp["n_comp"], "role": Dp["role"],
                           "err_physics": E["physics"], "err_pair": E["pair"]})
        for c in cands:
            fr[f"err_{c}"], fr[f"std_{c}"] = E[c], S[c]
        rows_out.append(fr)
    readings["any_gate_pass"] = bool(any(v["gate_pass"] for k, v in readings.items() if isinstance(v, dict) and "gate_pass" in v))
    readings["designation"] = "not made; user decision (no variant exposed in predict_fuel)"
    pd.concat(rows_out).to_csv(out / "rows_scored.csv.gz", index=False, float_format="%.6g")
    (out / "readings.json").write_text(json.dumps(readings, indent=1, default=float))
    log("locked stage done")


if __name__ == "__main__":
    main()
