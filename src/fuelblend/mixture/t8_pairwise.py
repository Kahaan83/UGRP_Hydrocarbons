"""
t8_pairwise - does pairwise additivity extrapolate to more components than you trained on?
=========================================================================================
Install at fuelblend_data/src/fuelblend/mixture/t8_pairwise.py  (needs your existing t8_deepsets.py,
t8_vs_hgb.py, common.py, cold_flow.py)

Model  (PairwiseMixture, "regular-solution" form, O(N) in the number of components)
    p_i  = enc(static_i, pure props_i(T), T)  in R^{3 x R}        (one R-dim vector per property)
    out  = sum_r lam_r * Var_x(p_r)  +  sum_i x_i u_i               x = mole fractions
    Var_x(p) = sum_{i<j} x_i x_j (p_i - p_j)^2   -> exactly pairwise-additive, zero for a pure fluid,
    the learned analogue of Scatchard-Hildebrand / Grunberg-Nissan excess terms. Residual over the physics
    baseline, same inputs / loss / splits as DeepSets, so the two are directly comparable. Because it never
    forms an N x N tensor it also runs on K = 170 components inside T9.

Experiment (the evidence for blends larger than your training data)
    python -m fuelblend.mixture.t8_pairwise --root . --train-max-comp 2 --seeds 3
    trains DeepSets and Pairwise on pure + binary rows only, then scores both on ternaries (extrapolation)
    and on held-out binaries (in-distribution). Prints the error-growth factor extrap/in-dist for each model
    (feed it to T9 as --unc-inflate) and a cluster bootstrap of Pairwise - DeepSets over system_id.
Production bundles for T9 (train on everything):
    python -m fuelblend.mixture.t8_pairwise --root . --train-max-comp 99 --save --seeds 5 --epochs 80
    -> results/pairwise/extrap_{deepsets,pairwise}_bundle.pt   (load_any_bundle handles both)
"""
from __future__ import annotations

import argparse
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider, RegistryTensors
from fuelblend.mixture.t8_deepsets import (MLP, SHORT, Cfg, DeepSetsMixture, build_data, compute_norm,
                                           group_holdout, load_bundle, make_batch, predict_rows, task_loss,
                                           to_device)
from fuelblend.mixture.t8_vs_hgb import cluster_bootstrap, metric


class PairwiseMixture(nn.Module):
    def __init__(self, d_static, cfg: Cfg, norm, rank=32):
        super().__init__()
        self.cfg, self.rank = cfg, rank
        d_tok = d_static + 2 + 2                              # static, pure(ln eta, rho), T-context
        self.enc = MLP(d_tok, cfg.hidden, 3 * rank, 2, cfg.dropout)
        self.unary = MLP(d_tok, cfg.hidden, 3, 1, cfg.dropout)
        self.lam = nn.Parameter(torch.zeros(3, rank))
        if cfg.mode == "residual":
            nn.init.zeros_(self.unary.net[-1].weight)
            nn.init.zeros_(self.unary.net[-1].bias)           # starts exactly at the physics baseline
        for k in ("res_sd", "y_mu", "y_sd"):
            self.register_buffer(k, norm[k].detach().clone().float().cpu())

    def forward(self, b):
        x, mask = b["x"], b["mask"]
        N = x.shape[1]
        tok = torch.cat([b["static"], b["pure"], b["ctx"][:, None, :].expand(-1, N, -1)], -1)
        p = self.enc(tok).view(*x.shape, 3, self.rank)        # (B,N,3,R)
        xm = (x * mask)[..., None, None]
        m1 = (xm * p).sum(1)
        m2 = (xm * p * p).sum(1)
        pair = ((m2 - m1 ** 2) * self.lam).sum(-1)            # (B,3)   variance (regular-solution) form
        uni = ((x * mask)[..., None] * self.unary(tok)).sum(1)
        out = pair + uni
        if self.cfg.mode == "residual":
            return b["base_raw"] + self.res_sd * out
        return self.y_mu + self.y_sd * out


def fit(model_fn, t, rt, rows_tr, rows_va, cfg, norm, seed, dev):
    torch.manual_seed(seed)
    model = model_fn(rt.S.shape[1], cfg, norm).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.wd)
    tr, va = torch.as_tensor(rows_tr, device=dev), torch.as_tensor(rows_va, device=dev)
    steps = math.ceil(len(tr) / cfg.bs)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=cfg.lr, total_steps=cfg.epochs * steps, pct_start=0.1)
    scale, tasks = norm["res_sd"].to(dev), list(cfg.tasks)
    lab = t["lab"][tr][:, tasks]
    pw = (lab.float() * lab.sum(0).clamp_min(1).float().pow(-cfg.balance)[None]).sum(1)
    best, best_state, bad = float("inf"), None, 0
    for _ in range(cfg.epochs):
        model.train()
        order = tr[torch.multinomial(pw, len(tr), replacement=True)] if cfg.balance > 0 \
            else tr[torch.randperm(len(tr), device=dev)]
        for s in range(0, len(order), cfg.bs):
            b = make_batch(t, rt, order[s:s + cfg.bs])
            loss = task_loss(model(b), b, scale, cfg)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            sched.step()
        if len(va):
            model.eval()
            with torch.no_grad():
                vl = float(sum(task_loss(model(b), b, scale, cfg) for b in
                               (make_batch(t, rt, va[s:s + 4096]) for s in range(0, len(va), 4096))))
            if vl < best - 1e-5:
                best, bad, best_state = vl, 0, {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
                if bad >= cfg.patience:
                    break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model.eval()


def load_any_bundle(path, rt, dev):
    """Drop-in replacement for t8_deepsets.load_bundle that also understands pairwise bundles."""
    z = torch.load(path, map_location="cpu", weights_only=False)
    if z.get("arch", "deepsets") == "deepsets":
        return load_bundle(path, rt, dev)
    c = z["cfg"]
    cfg = Cfg(**{**c, "tasks": tuple(c["tasks"]), "lam": tuple(c["lam"])})
    if z["n_registry"] != rt.S.shape[0]:
        raise ValueError("Bundle was trained with a different compounds.csv (registry size mismatch).")
    dummy = {"res_sd": torch.ones(3), "y_mu": torch.zeros(3), "y_sd": torch.ones(3)}
    models = []
    for st in z["states"]:
        m = PairwiseMixture(rt.S.shape[1], cfg, dummy, z.get("rank", 32)).to(dev)
        m.load_state_dict(st)
        models.append(m.eval())
    return {"models": models, "cfg": cfg, "tasks": z["tasks"], "conformal": z["conformal"],
            "comp_count": z["comp_count"], "n_max_train": z["n_max_train"], "arch": "pairwise"}


def save_bundle(path, arch, models, cfg, data, core, reg, rt, rank):
    meta = data.meta
    mix = meta["n_comp"].to_numpy()[core] >= 2
    idx = data.idx[core][mix]
    torch.save({"arch": arch, "cfg": asdict(cfg), "tasks": list(cfg.tasks), "conformal": [float("nan")] * 3,
                "comp_count": np.bincount(idx[idx >= 0], minlength=len(reg.ik)),
                "n_max_train": int(meta["n_comp"].to_numpy()[core].max()), "n_registry": len(reg.ik),
                "d_static": int(rt.S.shape[1]), "rank": rank,
                "states": [{k: v.detach().cpu() for k, v in m.state_dict().items()} for m in models]}, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="results/pairwise")
    ap.add_argument("--train-max-comp", type=int, default=2)
    ap.add_argument("--train-splits", nargs="+", default=["S1"])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--rank", type=int, default=32)
    ap.add_argument("--hidden", type=int, default=128)
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    a = ap.parse_args()

    root, out, dev = Path(a.root), Path(a.out), torch.device(a.device)
    out.mkdir(parents=True, exist_ok=True)
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    rt = RegistryTensors(reg, dev)
    data = build_data(root, reg, PurePropertyProvider(reg, root / "data/cache/pure_grid.npz"), rt, dev)
    t = to_device(data, dev)
    meta, lab = data.meta, data.lab
    split, n_comp = meta["split"].to_numpy(), meta["n_comp"].to_numpy()
    tr_splits = [s.upper() for s in a.train_splits]
    cfg = Cfg(mode="residual", seeds=a.seeds, epochs=a.epochs, hidden=a.hidden, embed=a.hidden)

    core = np.where(np.isin(split, tr_splits) & (n_comp <= a.train_max_comp) & lab.any(1))[0]
    sets = {"in_dist": np.where(~np.isin(split, tr_splits) & (split != "NONE") & (n_comp >= 2)
                                & (n_comp <= a.train_max_comp) & lab.any(1))[0],
            "extrap": np.where((n_comp > a.train_max_comp) & lab.any(1))[0]}
    print(f"[extrap] train rows {len(core)} (n_comp<={a.train_max_comp}); in-dist test {len(sets['in_dist'])}; "
          f"extrapolation test {len(sets['extrap'])}")
    if len(sets["extrap"]):
        print(pd.crosstab(n_comp[sets["extrap"]], meta["prop"].to_numpy()[sets["extrap"]]).rename(
            columns=dict(enumerate(SHORT))).to_string())
    norm, tasks = compute_norm(data, core), list(cfg.tasks)
    arch = {"deepsets": DeepSetsMixture, "pairwise": lambda d, c, n: PairwiseMixture(d, c, n, a.rank)}

    store = {}
    for name, fn in arch.items():
        models = []
        for s in range(a.seeds):
            fit_rows, va = group_holdout(meta, core, lab, cfg.val_frac, 1000 + s, tasks)
            models.append(fit(fn, t, rt, fit_rows, va, cfg, norm, s, dev))
        for sn, rows in sets.items():
            if len(rows):
                mu, _ = predict_rows(models, t, rt, torch.as_tensor(rows, device=dev))
                store[(name, sn)] = mu.cpu().numpy()
        if a.save:
            save_bundle(out / f"extrap_{name}_bundle.pt", name, models, cfg, data, core, reg, rt, a.rank)
        print(f"[{name}] trained {a.seeds} members")

    res = []
    for sn, rows in sets.items():
        for p in range(3):
            m = lab[rows, p] if len(rows) else np.zeros(0, bool)
            if m.sum() < 5:
                continue
            r = rows[m]
            y, sid = data.y[r, p], meta["system_id"].to_numpy()[r]
            errs = {"base": data.base_raw[r, p] - y, **{k: store[(k, sn)][m, p] - y for k in arch}}
            d, ci, pb = cluster_bootstrap(errs["pairwise"], errs["deepsets"], sid, p)
            res.append({"set": sn, "prop": SHORT[p], "n": int(m.sum()), "systems": len(np.unique(sid)),
                        **{k: metric(p, v) for k, v in errs.items()}, "pw_minus_ds": d, "ci_lo": ci[0],
                        "ci_hi": ci[1], "P(pw<ds)": pb})
    res = pd.DataFrame(res)
    pd.set_option("display.width", 220, "display.float_format", "{:.4f}".format)
    print("\n" + res.to_string(index=False))
    if {"in_dist", "extrap"} <= set(res["set"]):
        g = res.pivot(index="prop", columns="set", values=["base", "deepsets", "pairwise"])
        print("\nError-growth factor extrap / in-dist (use ~this as --unc-inflate scale in T9):")
        for k in ("base", "deepsets", "pairwise"):
            print(f"  {k:9s}", (g[k]["extrap"] / g[k]["in_dist"]).round(2).to_dict())
    res.to_csv(out / "extrapolation.csv", index=False)
    print("\nNegative pw_minus_ds = pairwise better. Interval containing 0 = not distinguishable.")


if __name__ == "__main__":
    main()
