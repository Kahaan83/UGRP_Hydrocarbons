"""
T9 v2 - Inverse blend design on the composition simplex
=======================================================
Replaces fuelblend/mixture/t9_inverse_design.py (needs the new cold_flow.py next to it; common.py and
t8_deepsets.py are unchanged).

NEW in v2
  * class_wt_pct   constraints: bound the mass % of aromatic / n_paraffin / iso_paraffin / cycloalkane / alkene_alkyne
  * wat_C          constraint: ideal-solubility wax appearance temperature (proxy for freezing / cloud point),
                   see cold_flow.py for its assumptions and limits
  * --reference    CSV of a real fuel's composition -> penalty on class marginals + carbon-number distribution
  * --min-neff     floor on the effective number of components (use for realistic 40+ component blends)
  * --drop-aux     reproduce v1 behaviour (only viscosity / density / flash constraints)
  * fixes: duplicate-keyword lam bug, .detach() before scalar conversion

Realistic diesel-like run (many components, reference-guided):
    python -m fuelblend.mixture.t9_inverse_design --bundle results/t8/<bundle>.pt --pool palette --preset diesel \
        --reference diesel_ref.csv --lam-ref 2 --sparsity 0 --min-frac 0.002 --min-neff 15

Aux constraint defaults (class windows, WAT limits, 5 K proxy margin) are ASSUMPTIONS - edit PRESETS or pass --spec.
Model-based constraints are tightened by k_sigma * ensemble std (+ conformal q for flash); aux constraints use their
own fixed `margin`.
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from scipy.optimize import minimize

from fuelblend.mixture.cold_flow import (CLASSES, carbon_matrix, class_matrix, compound_meta, load_reference,
                                         reference_distance, thermal_arrays, wax_appearance_C)
from fuelblend.mixture.common import (CompoundRegistry, PurePropertyProvider, RegistryTensors, assemble_batch,
                                      ctx_np, mole_fractions, physics_baselines)
from fuelblend.mixture.t8_pairwise import load_any_bundle as load_bundle   # DeepSets or pairwise bundles

SURROGATES = [("n-octane", "CCCCCCCC"), ("n-decane", "CCCCCCCCCC"), ("n-dodecane", "CCCCCCCCCCCC"),
              ("n-tetradecane", "CCCCCCCCCCCCCC"), ("isooctane", "CC(C)CC(C)(C)C"),
              ("isododecane", "CC(C)(C)CC(C)CC(C)(C)C"), ("cyclohexane", "C1CCCCC1"),
              ("methylcyclohexane", "CC1CCCCC1"), ("decalin", "C1CCC2CCCCC2C1"),
              ("bicyclohexyl", "C1CCC(CC1)C1CCCCC1"), ("toluene", "Cc1ccccc1"),
              ("ethylbenzene", "CCc1ccccc1"), ("m-xylene", "Cc1cccc(C)c1"),
              ("mesitylene", "Cc1cc(C)cc(C)c1"), ("1-methylnaphthalene", "Cc1cccc2ccccc12")]
COL = {"kin_visc_cSt": 0, "density_kg_m3": 1, "flash_C": 2}
AUX = {"class_wt_pct", "wat_C"}
DEFAULT_SCALE = {"kin_visc_cSt": 0.05, "density_kg_m3": 3.0, "flash_C": 1.5, "class_wt_pct": 2.0, "wat_C": 2.0}


@dataclass
class Constraint:
    prop: str                      # kin_visc_cSt | density_kg_m3 | flash_C | class_wt_pct | wat_C
    T_K: float | None = None       # None for flash / aux
    lo: float | None = None
    hi: float | None = None
    scale: float | None = None
    weight: float = 1.0
    group: str | None = None       # class name for class_wt_pct
    margin: float = 0.0            # fixed safety margin for aux constraints (metric units)

    def __post_init__(self):
        if self.scale is None:
            self.scale = DEFAULT_SCALE[self.prop]

    def bounds_m(self):            # bounds in metric space (ln for kinematic viscosity)
        f = math.log if self.prop == "kin_visc_cSt" else (lambda v: v)
        return (None if self.lo is None else f(self.lo)), (None if self.hi is None else f(self.hi))

    @property
    def key(self):
        if self.prop == "class_wt_pct":
            return f"class_wt_pct[{self.group}]"
        if self.prop == "wat_C":
            return "wat_C"
        return f"{self.prop}@{'-' if self.T_K is None else round(self.T_K - 273.15, 2)}C"


PRESETS = {
    "jet_a": [Constraint("kin_visc_cSt", 253.15, hi=8.0), Constraint("density_kg_m3", 288.15, lo=775.0, hi=840.0),
              Constraint("flash_C", None, lo=38.0),
              Constraint("class_wt_pct", group="aromatic", hi=25.0),            # D1655 limit is 25 vol%; ~wt% here
              Constraint("wat_C", hi=-40.0, margin=5.0)],                       # freezing point <= -40 C (proxy)
    "diesel": [Constraint("kin_visc_cSt", 313.15, lo=2.0, hi=4.5), Constraint("density_kg_m3", 288.15, lo=820.0, hi=845.0),
               Constraint("flash_C", None, lo=55.0),                            # EN 590-style limits
               Constraint("class_wt_pct", group="n_paraffin", lo=10.0, hi=35.0),
               Constraint("class_wt_pct", group="iso_paraffin", lo=10.0, hi=35.0),
               Constraint("class_wt_pct", group="cycloalkane", lo=15.0, hi=40.0),
               Constraint("class_wt_pct", group="aromatic", lo=15.0, hi=35.0),  # typical GC x GC windows (assumption)
               Constraint("wat_C", hi=0.0, margin=5.0)],                        # cloud-point-like limit (assumption)
}


def load_spec(preset, path, drop_aux=False):
    if path:
        out = []
        for d in json.load(open(path)):
            d = dict(d)
            T = d.pop("T_K", None)
            if "T_C" in d:
                T = d.pop("T_C") + 273.15
            out.append(Constraint(T_K=T, **d))
    else:
        out = list(PRESETS[preset])
    return [c for c in out if not (drop_aux and c.prop in AUX)]


def resolve_pool(reg, provider, arg):
    names, cost, idx = {}, {}, []
    if arg == "surrogate":
        from rdkit import Chem
        first = {k[:14]: i for i, k in enumerate(reg.ik)}
        for nm, smi in SURROGATES:
            ik = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
            i = reg.ik2i.get(ik, first.get(ik[:14]))
            if i is None:
                print(f"[pool] {nm} not found in compounds.csv - skipped")
                continue
            idx.append(i)
            names[i] = nm
    elif arg == "palette":
        idx = [i for i in range(len(reg.ik)) if reg.usable[i] and reg.has_flash[i]]
    else:
        df = pd.read_csv(arg)
        for r in df.itertuples(index=False):
            ik = getattr(r, "inchikey", None)
            if ik is None or ik not in reg.ik2i:
                from rdkit import Chem
                ik = Chem.MolToInchiKey(Chem.MolFromSmiles(r.smiles))
            i = reg.ik2i.get(ik)
            if i is not None:
                idx.append(i)
                cost[i] = float(getattr(r, "cost", 0.0))
    idx = [i for i in dict.fromkeys(idx) if reg.usable[i] and reg.has_flash[i]]
    provider.ensure(idx)
    idx = [i for i in idx if provider.valid[i]]
    return np.array(idx), names, cost


class BlendEvaluator:
    def __init__(self, bundle_paths, reg, rt, provider, pool_idx, spec, dev, aux, k_sigma=1.0,
                 conformal_props=("flash_C",), cost=None, lam=None, min_neff=0.0, flash_margin=None, unc_inflate=0.0):
        self.rt, self.dev, self.spec, self.k_sigma = rt, dev, spec, k_sigma
        self.unc_inflate = unc_inflate            # uncertainty grows with effective #components beyond training
        self.flash_margin = flash_margin          # fixed engineering margin (degC) replacing the bundle's conformal q
        self.conformal_props = set(conformal_props)
        self.K = len(pool_idx)
        self.idx = torch.as_tensor(pool_idx, dtype=torch.long, device=dev)
        t = lambda a: torch.as_tensor(np.asarray(a), dtype=torch.float32, device=dev)
        self.cls_oh, self.nc_oh, self.tm, self.hf = t(aux["cls_oh"]), t(aux["nc_oh"]), t(aux["tm"]), t(aux["hf"])
        self.ref = None if aux.get("ref") is None else (t(aux["ref"][0]), t(aux["ref"][1]))
        self.min_neff = min_neff
        self.owner = {}
        for p in bundle_paths:
            b = load_bundle(p, rt, dev)
            b["owned"] = []
            for tk in b["tasks"]:
                if tk not in self.owner:
                    self.owner[tk] = b
                    b["owned"].append(tk)
        need = {c for s in spec if s.prop not in AUX for c in ((0, 1) if s.prop == "kin_visc_cSt" else (COL[s.prop],))}
        missing = need - set(self.owner)
        if missing:
            raise ValueError(f"No T8 bundle provides head(s) {sorted(missing)} (0 visc, 1 dens, 2 flash); pass more --bundle files.")
        self.uniq = list({id(b): b for b in self.owner.values()}.values())
        counts = np.asarray(self.uniq[0]["comp_count"])[pool_idx]
        self.cov = t(np.minimum(1.0, np.log1p(counts) / math.log1p(20.0)))
        self.n_train_max = max(b["n_max_train"] for b in self.uniq)
        self.cost = t(cost if cost is not None else np.zeros(self.K))
        self.lam = {**dict(unc=0.05, dom=0.1, cost=0.0, ref=0.0, neff=1.0), **(lam or {})}
        self.pure, self.ctx = {}, {}
        for key in {c.T_K for c in spec if c.prop not in AUX}:
            T = 298.15 if key is None else key
            le, rho = provider.lookup(pool_idx, T)
            self.pure[key] = (t(le), t(rho))
            self.ctx[key] = torch.as_tensor(ctx_np([T], [key is None])[0], device=dev)

    def forward_T(self, w, key):
        B = w.shape[0]
        idx = self.idx[None].expand(B, -1)
        ln_i, rho_i = (v[None].expand(B, -1) for v in self.pure[key])
        base, valid, _ = physics_baselines(self.rt, idx, w, ln_i, rho_i)
        batch = assemble_batch(self.rt, idx, w, ln_i, rho_i, self.ctx[key][None].expand(B, -1), base, valid)
        res = {}
        for b in self.uniq:
            outs = torch.stack([m(batch) for m in b["models"]])
            mu = outs.mean(0)
            sd = outs.std(0) if outs.shape[0] > 1 else torch.zeros_like(mu)
            for p in b["owned"]:
                res[p] = (mu[:, p], sd[:, p])
        return res

    def aux_metric(self, c, w):
        if c.prop == "class_wt_pct":
            return 100.0 * (w @ self.cls_oh[:, CLASSES.index(c.group)])
        x = mole_fractions(self.rt, self.idx[None].expand(w.shape[0], -1), w)
        return wax_appearance_C(x, self.tm, self.hf)

    def metric(self, c, cols, w):
        if c.prop in AUX:
            mu = self.aux_metric(c, w)
            return mu, torch.zeros_like(mu), c.margin
        if c.prop == "kin_visc_cSt":
            (me, se), (mr, sr) = cols[0], cols[1]
            mr = mr.clamp_min(1.0)
            mu, sd = math.log(1000.0) + me - torch.log(mr), torch.sqrt(se ** 2 + (sr / mr) ** 2 + 1e-12)
        else:
            mu, sd = cols[COL[c.prop]]
        sd = sd * self.inflate(w)
        return mu, sd, self.k_sigma * sd + self.q_for(c)

    def inflate(self, w):
        """1 + g * max(0, n_eff - n_train_max) / n_train_max ; set g from the binary->ternary growth factor."""
        if self.unc_inflate <= 0:
            return 1.0
        n_eff = torch.exp(-(w * torch.log(w.clamp_min(1e-12))).sum(-1))
        return 1.0 + self.unc_inflate * F.relu(n_eff - self.n_train_max) / self.n_train_max

    def q_for(self, c):
        if c.prop not in self.conformal_props:
            return 0.0
        if c.prop == "flash_C" and self.flash_margin is not None:
            return float(self.flash_margin)
        col = 0 if c.prop == "kin_visc_cSt" else COL[c.prop]
        q = self.owner[col]["conformal"][col]
        return 0.0 if q is None or math.isnan(q) else float(q)

    def _cols(self, w, cache, c):
        if c.prop not in AUX and c.T_K not in cache:
            cache[c.T_K] = self.forward_T(w, c.T_K)
        return cache.get(c.T_K)

    def objective(self, w, sparsity=0.0):
        B = w.shape[0]
        J, unc, cache = w.new_zeros(B), w.new_zeros(B), {}
        n_mod = max(1, sum(c.prop not in AUX for c in self.spec))
        for c in self.spec:
            mu, sd, margin = self.metric(c, self._cols(w, cache, c), w)
            lo, hi = c.bounds_m()
            v = w.new_zeros(B)
            if lo is not None:
                v = v + F.relu(lo - (mu - margin))
            if hi is not None:
                v = v + F.relu((mu + margin) - hi)
            J = J + c.weight * (v / c.scale) ** 2
            if c.prop not in AUX:
                unc = unc + sd / c.scale
        J = J + self.lam["unc"] * unc / n_mod
        J = J + self.lam["dom"] * (w * (1.0 - self.cov)).sum(-1)
        J = J + self.lam["cost"] * (w * self.cost).sum(-1)
        if self.ref is not None and self.lam["ref"] > 0:
            J = J + self.lam["ref"] * reference_distance(w, self.cls_oh, self.nc_oh, *self.ref)
        if self.min_neff > 0:
            n_eff = torch.exp(-(w * torch.log(w.clamp_min(1e-12))).sum(-1))
            J = J + self.lam["neff"] * F.relu(self.min_neff - n_eff) / self.min_neff
        if sparsity:
            J = J + sparsity * (torch.sqrt(w + 1e-8).sum(-1) - 1.0)
        return J

    @torch.no_grad()
    def report(self, w):
        wt = torch.as_tensor(w, dtype=torch.float32, device=self.dev)[None]
        cache, out = {}, {}
        for c in self.spec:
            mu, sd, m = self.metric(c, self._cols(wt, cache, c), wt)
            mu, sd, m = float(mu[0]), float(sd[0]), float(m)
            lo, hi = c.bounds_m()
            nat = math.exp if c.prop == "kin_visc_cSt" else (lambda v: v)
            out[c.key] = dict(value=nat(mu), lo=nat(mu - m), hi=nat(mu + m), std_metric=sd,
                              nominal_ok=bool((lo is None or mu >= lo) and (hi is None or mu <= hi)),
                              robust_ok=bool((lo is None or mu - m >= lo) and (hi is None or mu + m <= hi)))
        return out


def prune(w, k, min_frac):
    order = np.argsort(-w)
    act = order[:k] if k else order
    act = act[w[act] >= min_frac]
    return np.sort(act if len(act) else order[:1])


def polish(ev, w, active, iters):
    act = torch.as_tensor(active, device=ev.dev)

    def fg(v):
        vt = torch.tensor(v, dtype=torch.float32, device=ev.dev, requires_grad=True)
        wf = torch.zeros(1, ev.K, device=ev.dev)
        wf[0, act] = vt
        J = ev.objective(wf).sum()
        J.backward()
        return float(J.detach()), vt.grad.double().cpu().numpy()

    v0 = w[active] / w[active].sum()
    r = minimize(fg, v0, jac=True, method="SLSQP", bounds=[(0.0, 1.0)] * len(active),
                 constraints=[{"type": "eq", "fun": lambda v: v.sum() - 1.0, "jac": lambda v: np.ones_like(v)}],
                 options={"maxiter": iters, "ftol": 1e-9})
    out = np.zeros(ev.K)
    v = np.clip(r.x, 0, None)
    out[active] = v / v.sum()
    return out


def optimise(ev, a):
    torch.manual_seed(a.seed)
    K, dev = ev.K, ev.dev
    half = a.starts // 2
    init = torch.cat([torch.distributions.Dirichlet(torch.full((K,), a.alpha)).sample((half,)),
                      torch.softmax(0.3 * torch.randn(a.starts - half, K), -1)])
    theta = torch.log(init.clamp_min(1e-4)).to(dev).requires_grad_(True)
    opt = torch.optim.Adam([theta], lr=a.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, a.iters)
    for it in range(a.iters):
        J = ev.objective(torch.softmax(theta, -1), sparsity=a.sparsity * min(1.0, it / (0.6 * a.iters)))
        opt.zero_grad()
        J.sum().backward()
        opt.step()
        sched.step()
        if it % 100 == 0:
            print(f"[opt] iter {it:4d}  best J {float(J.detach().min()):.4f}")
    with torch.no_grad():
        W = torch.softmax(theta, -1)
        J0 = ev.objective(W).cpu().numpy()
    W = W.detach().cpu().numpy()
    finals, seen = [], set()
    for s in np.argsort(J0)[: a.polish * 4]:
        act = prune(W[s], a.max_components, a.min_frac)
        w = polish(ev, W[s], act, a.polish_iters)
        w[w < a.min_frac] = 0.0
        if a.max_components and (w > 0).sum() > a.max_components:
            w[np.argsort(-w)[a.max_components:]] = 0.0
        w = w / w.sum()
        key = tuple(np.round(w, 2))
        if key in seen:
            continue
        seen.add(key)
        finals.append(w)
        if len(finals) >= a.polish:
            break
    return finals


def verify_with_api(reg, w, pool_idx, spec):
    try:
        from fuelblend.api.predict_blend import predict
    except Exception as e:                                       # pragma: no cover
        return {"error": str(e)}
    comp = [(reg.ik[int(pool_idx[i])], float(w[i])) for i in np.where(w > 0)[0]]
    out = {}
    for c in spec:
        if c.prop in AUX:
            continue
        try:
            r = predict(comp, T_K=c.T_K or 298.15)
            if c.prop == "flash_C":
                out[c.key] = {"value": r["flash_point_C"]["value"], "interval": r["flash_point_C"].get("interval")}
            elif c.prop == "density_kg_m3":
                out[c.key] = {"value": r["density_kg_m3"]["value"]}
            else:
                out[c.key] = {"value": 1000.0 * r["viscosity_dyn_mPas"]["value"] / r["density_kg_m3"]["value"]}
        except Exception as e:
            out[c.key] = {"error": str(e)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--bundle", nargs="+", required=True)
    ap.add_argument("--pool", default="surrogate", help="surrogate | palette | path/to/pool.csv")
    ap.add_argument("--preset", default="jet_a", choices=list(PRESETS))
    ap.add_argument("--spec", default=None, help="JSON list of constraints (overrides --preset)")
    ap.add_argument("--drop-aux", action="store_true", help="ignore class-fraction and wax-appearance constraints")
    ap.add_argument("--reference", default=None, help="CSV: inchikey|smiles, mass_frac|wt_pct of a real fuel")
    ap.add_argument("--lam-ref", type=float, default=1.0)
    ap.add_argument("--min-neff", type=float, default=0.0)
    ap.add_argument("--unc-inflate", type=float, default=0.0,
                    help="g in sd *= 1 + g*max(0, n_eff - n_train_max)/n_train_max; take it from t8_pairwise's growth factor")
    ap.add_argument("--max-components", type=int, default=None)
    ap.add_argument("--min-frac", type=float, default=0.005)
    ap.add_argument("--k-sigma", type=float, default=1.0)
    ap.add_argument("--conformal-props", nargs="*", default=["flash_C"])
    ap.add_argument("--flash-margin", type=float, default=5.0,
                    help="fixed flash safety margin in degC (default 5). Negative = use the bundle's conformal q, "
                         "which is only trustworthy when it was calibrated on >= ~50 flash rows")
    ap.add_argument("--min-train-mix", type=int, default=0,
                    help="drop pool compounds that appear in fewer than N training MIXTURE rows of the T8 bundle")
    ap.add_argument("--lam-unc", type=float, default=0.05)
    ap.add_argument("--lam-dom", type=float, default=0.1)
    ap.add_argument("--lam-cost", type=float, default=0.0)
    ap.add_argument("--sparsity", type=float, default=0.02)
    ap.add_argument("--starts", type=int, default=128)
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--lr", type=float, default=0.08)
    ap.add_argument("--alpha", type=float, default=0.3)
    ap.add_argument("--polish", type=int, default=6)
    ap.add_argument("--polish-iters", type=int, default=150)
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/t9")
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    a = ap.parse_args()

    root, dev = Path(a.root), torch.device(a.device)
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    rt = RegistryTensors(reg, dev)
    provider = PurePropertyProvider(reg, root / "data/cache/pure_grid.npz")
    spec = load_spec(a.preset, a.spec, a.drop_aux)
    pool, names, cost = resolve_pool(reg, provider, a.pool)
    z = torch.load(a.bundle[0], map_location="cpu", weights_only=False)
    print("[conformal] bundle q (visc ln-units, dens kg/m3, flash degC) = "
          + ", ".join(f"{q:.3g}" for q in z["conformal"]) + f";  max components in training mixtures = {z['n_max_train']}")
    if z["conformal"][2] == z["conformal"][2] and z["conformal"][2] > 20:
        print("[conformal] WARNING: flash q > 20 degC - calibration set is too small; using --flash-margin instead "
              "(pass a negative value to override)")
    if a.min_train_mix > 0:
        keep = np.asarray(z["comp_count"])[pool] >= a.min_train_mix
        print(f"[pool] --min-train-mix {a.min_train_mix}: kept {int(keep.sum())}/{len(pool)} compounds")
        pool = pool[keep]

    cdf, cls, nc = compound_meta(root / "data/unified/compounds.csv")
    tm, hf = thermal_arrays(cdf, cls, nc, root / "data/unified/thermal.csv")
    aux = {"cls_oh": class_matrix([cls[i] for i in pool]), "nc_oh": carbon_matrix(nc[pool]),
           "tm": tm[pool], "hf": hf[pool],
           "ref": load_reference(a.reference, reg, cls, nc) if a.reference else None}
    n_th = int(np.isfinite(aux["tm"]).sum())
    print(f"[pool] {len(pool)} usable components ({n_th} with melting data for the wax proxy; the rest are "
          f"assumed never to precipitate). Spec: " + "; ".join(f"{c.key} in [{c.lo}, {c.hi}]" for c in spec))
    ev = BlendEvaluator(a.bundle, reg, rt, provider, pool, spec, dev, aux, a.k_sigma, a.conformal_props,
                        cost=np.array([cost.get(int(i), 0.0) for i in pool]) if cost else None,
                        lam={"unc": a.lam_unc, "dom": a.lam_dom, "cost": a.lam_cost, "ref": a.lam_ref},
                        min_neff=a.min_neff, flash_margin=None if a.flash_margin < 0 else a.flash_margin,
                        unc_inflate=a.unc_inflate)
    for c in spec:
        if c.prop not in AUX:
            print(f"[margin] {c.key}: k_sigma={a.k_sigma}, conformal q={ev.q_for(c):.3g}")

    finals = optimise(ev, a)
    rows = []
    for w in finals:
        rep = ev.report(w)
        J = float(ev.objective(torch.as_tensor(w, dtype=torch.float32, device=dev)[None]).detach())
        rows.append((not all(v["robust_ok"] for v in rep.values()), J, w, rep))
    rows.sort(key=lambda r: (r[0], r[1]))

    out = []
    for rank, (infeasible, J, w, rep) in enumerate(rows[: a.top], 1):
        act = np.where(w > 0)[0]
        n_eff = float(np.exp(-(w[act] * np.log(w[act])).sum()))
        comp = sorted(((names.get(int(pool[i]), reg.smiles[int(pool[i])]), 100 * w[i], int(pool[i])) for i in act),
                      key=lambda t: -t[1])
        cls_mass = {k: float(100 * (w @ aux["cls_oh"][:, j])) for j, k in enumerate(CLASSES)}
        print(f"\n=== Blend #{rank}  J={J:.4f}  robustly feasible={not infeasible}  "
              f"n_active={len(act)}  n_eff={n_eff:.1f} ===")
        for nm, pct, _ in comp:
            print(f"  {pct:6.2f} wt%  {nm}")
        print("  class wt%: " + ", ".join(f"{k} {v:.1f}" for k, v in cls_mass.items()))
        for k, v in rep.items():
            print(f"  {k:30s} {v['value']:9.3f}  safe-interval [{v['lo']:.3f}, {v['hi']:.3f}]  "
                  f"nominal_ok={v['nominal_ok']} robust_ok={v['robust_ok']}")
        ver = {} if a.no_verify else verify_with_api(reg, w, pool, spec)
        if ver:
            print("  consistency check (predict_blend API, same pure models/data - NOT independent):",
                  json.dumps(ver, default=float))
        warn = []
        if len(act) > ev.n_train_max:
            warn.append(f"{len(act)} active components > {ev.n_train_max} max seen in training mixtures (extrapolation)")
        unseen = [nm for nm, _, i in comp if ev.uniq[0]["comp_count"][i] == 0]
        if unseen:
            warn.append("never seen in a training mixture: " + ", ".join(unseen))
        no_th = 100 * sum(w[i] for i in act if not np.isfinite(aux["tm"][i]))
        if any(c.prop == "wat_C" for c in spec) and no_th > 0:
            warn.append(f"{no_th:.1f} wt% has no melting data (assumed never to precipitate in the wax proxy)")
        for m in warn:
            print("  WARNING:", m)
        out.append({"rank": rank, "J": J, "robust_feasible": not infeasible, "n_active": len(act), "n_eff": n_eff,
                    "class_wt_pct": cls_mass,
                    "composition": [{"name": nm, "mass_pct": pct, "inchikey": reg.ik[i]} for nm, pct, i in comp],
                    "predicted": rep, "api_check": ver, "warnings": warn})
    Path(a.out).mkdir(parents=True, exist_ok=True)
    tag = Path(a.spec).stem if a.spec else a.preset
    json.dump(out, open(Path(a.out) / f"t9_{tag}.json", "w"), indent=2, default=float)
    print(f"\n[done] wrote {Path(a.out) / f't9_{tag}.json'}")


if __name__ == "__main__":
    main()
