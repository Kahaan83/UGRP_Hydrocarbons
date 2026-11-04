"""
T8 - Permutation-invariant DeepSets mixture model vs physics + residual baselines
==================================================================================
Install at fuelblend_data/src/fuelblend/mixture/t8_deepsets.py   (next to common.py)

Quick smoke test (minutes):
    python -m fuelblend.mixture.t8_deepsets --root . --seeds 1 --epochs 3 --device cpu
Real run (multi-task, residual over physics, both pooling types):
    python -m fuelblend.mixture.t8_deepsets --root . --pool both --mode residual
Ablation grid (tasks x pooling x residual/direct), writes ablation.csv:
    python -m fuelblend.mixture.t8_deepsets --root . --ablate --seeds 3

Design
  token_i  = [static descriptors, pure ln eta_i(T), pure rho_i(T), w_i, x_i, T-context]
  phi      = per-component MLP;      z = [sum_i w_i phi_i | sum_i x_i phi_i | sum_i w_i attn_i]
  rho      = MLP([z, physics baselines + validity flags, T-context])
  heads    = 3 linear heads (ln eta, rho, flash);  residual mode:  y = baseline + sd * head(z)
  attention is biased by log(w_i), so near-zero-fraction components are ignored smoothly
  (this keeps dense 170-component evaluation in T9 consistent with sparse training mixtures).
Loss: Huber on residual scale per property, masked (each row carries ONE property label),
      quality-weighted, balanced sampling across properties, asymmetric penalty for flash
      OVER-prediction (safety-critical). Ensemble over seeds -> epistemic std. Split-conformal
      absolute-error quantiles per property from a held-out calibration set of systems.
Evaluation: fixed splits from splits/s1s2s3.json, grouped by system_id; baseline vs DeepSets on
      identical rows; mixtures (n_comp>=2) reported separately from pure rows.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from fuelblend.mixture.common import (PROPS, T_HI, T_LO, CompoundRegistry, PurePropertyProvider,
                                      RegistryTensors, assemble_batch, ctx_np, physics_baselines)

SHORT = ("visc", "dens", "flash")
EXCLUDE_SYSTEMS = {"a175950749e8"}      # known 1000x viscosity reporting error (T3 audit); only the gross-outlier rows are dropped
REFERENCE = {                           # from PROJECT_PROGRESS_SUMMARY.md (their row sets may differ slightly)
    "visc S2 log-RMSE": "physics 0.1225 -> HGB residual 0.0781",
    "dens S2 RMSE kg/m3": "ideal-volume 2.48 -> residual 2.46",
    "flash S2 MAE degC": "Le Chatelier 1.354 -> +median shift 0.835 (conformal cov 81.8%, over>5C 0.0%)",
}


# --------------------------------------------------------------------------- config / model
@dataclass
class Cfg:
    pool: str = "both"            # sum | attn | both
    mode: str = "residual"        # residual | direct   (baseline is an input feature in both)
    tasks: tuple = (0, 1, 2)      # 0 visc, 1 dens, 2 flash
    hidden: int = 128
    embed: int = 128
    heads: int = 4
    dropout: float = 0.1
    epochs: int = 80
    bs: int = 512
    lr: float = 2e-3
    wd: float = 1e-4
    patience: int = 12
    balance: float = 0.5          # sampling weight ~ n_prop^-balance (0 = natural frequencies)
    flash_over_weight: float = 3.0
    lam: tuple = (1.0, 1.0, 1.0)
    seeds: int = 5
    val_frac: float = 0.10
    cal_frac: float = 0.15


class MLP(nn.Module):
    def __init__(self, d_in, d_hidden, d_out, depth=2, dropout=0.0):
        super().__init__()
        layers, d = [], d_in
        for _ in range(depth):
            layers += [nn.Linear(d, d_hidden), nn.LayerNorm(d_hidden), nn.SiLU(), nn.Dropout(dropout)]
            d = d_hidden
        layers.append(nn.Linear(d, d_out))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class FractionBiasedAttention(nn.Module):
    """Single self-attention layer over components; keys are biased by log(mass fraction)."""

    def __init__(self, d, heads):
        super().__init__()
        self.h, self.dh = heads, d // heads
        self.qkv = nn.Linear(d, 3 * d)
        self.out = nn.Linear(d, d)

    def forward(self, h, logw, mask):
        B, N, D = h.shape
        q, k, v = self.qkv(h).view(B, N, 3, self.h, self.dh).permute(2, 0, 3, 1, 4)
        s = q @ k.transpose(-1, -2) / math.sqrt(self.dh) + logw[:, None, None, :]
        s = s.masked_fill(~mask[:, None, None, :], -1e9)
        a = torch.softmax(s, -1) @ v
        return self.out(a.transpose(1, 2).reshape(B, N, D))


class DeepSetsMixture(nn.Module):
    def __init__(self, d_static, cfg: Cfg, norm):
        super().__init__()
        self.cfg = cfg
        self.use_sum = cfg.pool in ("sum", "both")
        self.use_attn = cfg.pool in ("attn", "both")
        self.phi = MLP(d_static + 2 + 2 + 2, cfg.hidden, cfg.embed, 2, cfg.dropout)
        if self.use_attn:
            self.attn = FractionBiasedAttention(cfg.embed, cfg.heads)
        n_pool = (2 * cfg.embed if self.use_sum else 0) + (cfg.embed if self.use_attn else 0)
        self.rho = MLP(n_pool + 6 + 2, cfg.hidden, cfg.hidden, 2, cfg.dropout)
        self.heads = nn.ModuleList([nn.Linear(cfg.hidden, 1) for _ in range(3)])
        if cfg.mode == "residual":
            for hd in self.heads:
                nn.init.zeros_(hd.weight)
                nn.init.zeros_(hd.bias)              # starts exactly at the physics baseline
        for k in ("res_sd", "y_mu", "y_sd"):
            self.register_buffer(k, norm[k].detach().clone().float().cpu())

    def forward(self, b):
        w, x, mask = b["w"], b["x"], b["mask"]
        N = w.shape[1]
        tok = torch.cat([b["static"], b["pure"], w[..., None], x[..., None],
                         b["ctx"][:, None, :].expand(-1, N, -1)], -1)
        h = self.phi(tok)
        pooled = []
        if self.use_sum:
            pooled += [(w[..., None] * h).sum(1), (x[..., None] * h).sum(1)]
        if self.use_attn:
            a = self.attn(h, torch.log(w.clamp_min(1e-6)), mask)
            pooled.append((w[..., None] * a).sum(1))
        z = self.rho(torch.cat(pooled + [b["base_in"], b["ctx"]], -1))
        out = torch.cat([hd(z) for hd in self.heads], -1)
        if self.cfg.mode == "residual":
            return b["base_raw"] + self.res_sd * out
        return self.y_mu + self.y_sd * out


# --------------------------------------------------------------------------- data
@dataclass
class MixData:
    meta: pd.DataFrame
    idx: np.ndarray
    w: np.ndarray
    ln_eta_i: np.ndarray
    rho_i: np.ndarray
    ctx: np.ndarray
    base_raw: np.ndarray
    base_valid: np.ndarray
    y: np.ndarray
    lab: np.ndarray
    wt: np.ndarray


def load_split_labels(df, path):
    if "split" in df.columns:
        return df["split"].astype(str).str.upper().to_numpy()
    j = json.load(open(path))
    if isinstance(j, dict) and isinstance(j.get("splits"), dict):
        j = j["splits"]
    sid = df["system_id"].to_numpy()
    if all(isinstance(v, (list, tuple)) for v in j.values()):            # {"S1": [system_id,...], ...}
        first_val = next((v for v in j.values() if len(v) > 0), None)
        if first_val is not None and isinstance(first_val[0], int):
            out = np.full(len(df), "NONE", dtype=object)
            for name, ids in j.items():
                common = set(ids).intersection(df.index)
                if common:
                    out[df.index.isin(common)] = str(name).upper()
            return out
        out = np.full(len(df), "NONE", dtype=object)
        for name, ids in j.items():
            out[np.isin(sid, list(ids))] = str(name).upper()
        return out
    if all(isinstance(v, str) for v in j.values()):                       # {system_id: "S1", ...}
        return np.array([str(j.get(s, "NONE")).upper() for s in sid], dtype=object)
    raise ValueError(f"Unrecognised splits JSON layout (top-level keys: {list(j)[:10]}). "
                     "Edit load_split_labels() to match it.")


def build_data(root, reg, provider, rt, dev, drop_conflict=True) -> MixData:
    df = pd.read_csv(root / "data/processed/model_table.csv")
    splits_path = root / "splits/s1s2s3.json"
    if "split" not in df.columns and splits_path.exists():
        j = json.load(open(splits_path))
        split_col = np.full(len(df), "NONE", dtype=object)
        for s_name in ("S1", "S2", "S3"):
            if s_name in j:
                idx_list = j[s_name]
                if len(idx_list) > 0 and isinstance(idx_list[0], int):
                    split_col[idx_list] = s_name
                elif len(idx_list) > 0 and isinstance(idx_list[0], str):
                    split_col[df["system_id"].isin(idx_list)] = s_name
        df["split"] = split_col
    if "is_hc" in df.columns:
        df = df[df["is_hc"].astype(bool)]
    df = df[df["property"].isin(PROPS)]
    if drop_conflict and "conflict" in df:
        df = df[~df["conflict"].fillna(False).astype(bool)]
    p = df["property"].map({k: i for i, k in enumerate(PROPS)}).to_numpy()
    T, val = df["T_K"].to_numpy(float), df["value"].to_numpy(float)
    ok = np.isfinite(val) & ((p == 2) | (np.isfinite(T) & (T >= T_LO) & (T <= T_HI) & (val > 0)))
    df = df[ok].reset_index(drop=True)

    keep, ids_l, w_l = [], [], []
    for i, (a, b) in enumerate(zip(df["comp_inchikeys"], df["w_mass"])):
        ids, ws = [reg.ik2i.get(k, -1) for k in json.loads(a)], json.loads(b)
        if len(ids) != len(ws) or min(ids) < 0 or not reg.usable[ids].all() or sum(ws) <= 0:
            continue
        keep.append(i)
        ids_l.append(ids)
        w_l.append([v / sum(ws) for v in ws])
    print(f"[data] rows kept {len(keep)}/{len(df)} (unknown / unusable compounds dropped)")
    df = df.iloc[keep].reset_index(drop=True)
    R, N = len(df), max(map(len, ids_l))
    idx, w = np.full((R, N), -1, dtype=np.int64), np.zeros((R, N))
    for r, (a, b) in enumerate(zip(ids_l, w_l)):
        idx[r, :len(a)], w[r, :len(b)] = a, b

    p = df["property"].map({k: i for i, k in enumerate(PROPS)}).to_numpy()
    val = df["value"].to_numpy(float)
    miss = p == 2                                         # flash rows carry no temperature
    T = np.where(miss, 298.15, df["T_K"].to_numpy(float))
    ctx = ctx_np(T, miss)

    provider.ensure(np.unique(idx[idx >= 0]))
    ln_i, rho_i = provider.lookup(idx, T[:, None])

    base, valid = np.zeros((R, 3)), np.zeros((R, 3), bool)
    f32 = lambda a: torch.as_tensor(a, dtype=torch.float32, device=dev)
    with torch.no_grad():
        for s in range(0, R, 8192):
            sl = slice(s, s + 8192)
            b, v, _ = physics_baselines(rt, torch.as_tensor(idx[sl], device=dev), f32(w[sl]),
                                        f32(ln_i[sl]), f32(rho_i[sl]))
            base[sl], valid[sl] = b.cpu().numpy(), v.cpu().numpy()

    r_ = np.arange(R)
    y = np.full((R, 3), np.nan)
    y[r_, p] = np.where(p == 0, np.log(np.clip(val, 1e-12, None)), val)
    lab = np.zeros((R, 3), bool)
    lab[r_, p] = valid[r_, p] & np.isfinite(y[r_, p])
    bad = (df["system_id"].isin(EXCLUDE_SYSTEMS).to_numpy() & (p == 0)
           & (np.abs(y[r_, p] - base[r_, p]) > math.log(100)))
    lab[r_[bad], p[bad]] = False
    print(f"[data] labelled rows per property {dict(zip(SHORT, lab.sum(0)))}; "
          f"dropped for invalid baseline: {int((~valid[r_, p]).sum())}; gross outliers: {int(bad.sum())}")

    meta = pd.DataFrame({"system_id": df["system_id"], "prop": p, "T_K": df["T_K"], "kind": df["kind"],
                         "n_comp": [len(a) for a in ids_l], "value": val,
                         "split": load_split_labels(df, root / "splits/s1s2s3.json")})
    print(pd.crosstab(meta["split"], meta["prop"].map(dict(enumerate(SHORT)))).to_string())
    wt = df["quality_weight"].fillna(1.0).to_numpy(float) if "quality_weight" in df else np.ones(R)
    return MixData(meta, idx, w, ln_i, rho_i, ctx, base, valid, y, lab, wt)


def to_device(d: MixData, dev):
    f = lambda a, dt=torch.float32: torch.as_tensor(np.array(a), dtype=dt, device=dev)
    return dict(idx=f(d.idx, torch.long), w=f(d.w), ln_eta_i=f(np.nan_to_num(d.ln_eta_i)),
                rho_i=f(np.nan_to_num(d.rho_i, nan=750.0)), ctx=f(d.ctx),
                base_raw=f(np.nan_to_num(d.base_raw)), base_valid=f(d.base_valid, torch.bool),
                y=f(np.nan_to_num(d.y)), lab=f(d.lab, torch.bool), wt=f(d.wt))


def make_batch(t, rt, rows):
    b = assemble_batch(rt, t["idx"][rows], t["w"][rows], t["ln_eta_i"][rows], t["rho_i"][rows],
                       t["ctx"][rows], t["base_raw"][rows], t["base_valid"][rows])
    b.update(y=t["y"][rows], lab=t["lab"][rows], wt=t["wt"][rows])
    return b


def group_holdout(meta, rows, lab, frac, seed, tasks):
    """Hold out whole systems (per property, so small properties keep representation)."""
    rng, held = np.random.default_rng(seed), set()
    sids = meta["system_id"].to_numpy()
    for p in tasks:
        u = np.unique(sids[rows[lab[rows, p]]])
        if len(u) >= 4:
            held |= set(rng.choice(u, size=max(1, int(round(frac * len(u)))), replace=False))
    m = np.isin(sids[rows], list(held))
    return rows[~m], rows[m]


# --------------------------------------------------------------------------- training
def task_loss(pred, b, scale, cfg: Cfg):
    total = pred.new_zeros(())
    for p in cfg.tasks:
        m = b["lab"][:, p]
        if not m.any():
            continue
        e = (pred[:, p] - b["y"][:, p]) / scale[p]
        l = F.smooth_l1_loss(e, torch.zeros_like(e), reduction="none", beta=1.0)
        if p == 2 and cfg.flash_over_weight != 1.0:
            l = l * torch.where(e > 0, torch.full_like(e, cfg.flash_over_weight), torch.ones_like(e))
        wm = b["wt"] * m
        total = total + cfg.lam[p] * (l * wm).sum() / wm.sum().clamp_min(1e-9)
    return total


@torch.no_grad()
def predict_rows(models, t, rt, rows, bs=4096):
    for m in models:
        m.eval()
    P = torch.cat([torch.stack([m(make_batch(t, rt, rows[s:s + bs])) for m in models])
                   for s in range(0, len(rows), bs)], 1)
    return P.mean(0), (P.std(0) if P.shape[0] > 1 else torch.zeros_like(P[0]))


def fit_member(t, rt, rows_tr, rows_va, cfg, norm, seed, dev):
    torch.manual_seed(seed)
    model = DeepSetsMixture(rt.S.shape[1], cfg, norm).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.wd)
    tr = torch.as_tensor(rows_tr, device=dev)
    va = torch.as_tensor(rows_va, device=dev)
    steps = math.ceil(len(tr) / cfg.bs)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=cfg.lr, total_steps=cfg.epochs * steps, pct_start=0.1)
    scale = norm["res_sd"].to(dev)
    tasks = list(cfg.tasks)
    lab = t["lab"][tr][:, tasks]
    npk = lab.sum(0).clamp_min(1).float()
    pw = (lab.float() * npk.pow(-cfg.balance)[None]).sum(1)
    best, best_state, bad = float("inf"), None, 0
    for ep in range(cfg.epochs):
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
                best, bad = vl, 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
                if bad >= cfg.patience:
                    break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model.eval()


def compute_norm(data: MixData, rows):
    res, mu, sd = [], [], []
    for p in range(3):
        r = rows[data.lab[rows, p]]
        if len(r) < 5:
            res.append(1.0), mu.append(0.0), sd.append(1.0)
            continue
        y, b = data.y[r, p], data.base_raw[r, p]
        res.append(max(float((y - b).std()), 1e-3)), mu.append(float(y.mean())), sd.append(max(float(y.std()), 1e-3))
    return {k: torch.tensor(v, dtype=torch.float32) for k, v in (("res_sd", res), ("y_mu", mu), ("y_sd", sd))}


# --------------------------------------------------------------------------- evaluation
def _stats(p, e):
    e = np.asarray(e)
    if len(e) == 0:
        return {}
    if p == 0:
        return {"rmse": float(np.sqrt((e ** 2).mean())), "mape_pct": float(100 * np.abs(np.expm1(e)).mean())}
    if p == 1:
        return {"rmse": float(np.sqrt((e ** 2).mean())), "mae": float(np.abs(e).mean())}
    return {"mae": float(np.abs(e).mean()), "rmse": float(np.sqrt((e ** 2).mean())),
            "bias_med": float(np.median(e)), "over5_pct": float(100 * (e > 5).mean())}


def conformal_q(err, alpha=0.1):
    n = len(err)
    if n < 10:
        return float("nan")
    return float(np.quantile(np.abs(err), min(1.0, math.ceil((n + 1) * (1 - alpha)) / n), method="higher"))


def run_experiment(data, t, rt, cfg, train_splits, dev, out_dir, tag, reg):
    meta, lab, tasks = data.meta, data.lab, list(cfg.tasks)
    split = meta["split"].to_numpy()
    usable = lab[:, tasks].any(1)
    tr_all = np.where(np.isin(split, train_splits) & usable)[0]
    core, cal = group_holdout(meta, tr_all, lab, cfg.cal_frac, 12345, tasks)
    norm = compute_norm(data, core)
    t0, models = time.time(), []
    for s in range(cfg.seeds):
        fit, va = group_holdout(meta, core, lab, cfg.val_frac, 1000 + s, tasks)
        models.append(fit_member(t, rt, fit, va, cfg, norm, s, dev))
    print(f"[{tag}] trained {cfg.seeds} members in {time.time() - t0:.0f}s on {len(core)} rows "
          f"(cal {len(cal)}, val ~{cfg.val_frac:.0%} of systems)")

    qs = []
    if len(cal):
        mu_c, _ = predict_rows(models, t, rt, torch.as_tensor(cal, device=dev))
        mu_c = mu_c.cpu().numpy()
        for p in range(3):
            m = lab[cal, p] if p in tasks else np.zeros(len(cal), bool)
            qs.append(conformal_q(mu_c[m, p] - data.y[cal[m], p]) if m.sum() else float("nan"))
    else:
        qs = [float("nan")] * 3

    res, preds = [], []
    for sp in sorted(set(split) - set(train_splits) - {"NONE"}):
        rows = np.where(split == sp)[0]
        if not len(rows):
            continue
        mu, sd = predict_rows(models, t, rt, torch.as_tensor(rows, device=dev))
        mu, sd = mu.cpu().numpy(), sd.cpu().numpy()
        for p in tasks:
            for subset, sel in (("mix", meta["n_comp"].to_numpy()[rows] >= 2), ("all", np.ones(len(rows), bool))):
                m = lab[rows, p] & sel
                if not m.any():
                    continue
                r = rows[m]
                eb, em = data.base_raw[r, p] - data.y[r, p], mu[m, p] - data.y[r, p]
                row = {"tag": tag, "split": sp, "subset": subset, "prop": SHORT[p], "n": int(m.sum()),
                       "systems": int(meta["system_id"].iloc[r].nunique())}
                row.update({f"base_{k}": v for k, v in _stats(p, eb).items()})
                row.update({f"ds_{k}": v for k, v in _stats(p, em).items()})
                if not math.isnan(qs[p]):
                    row["ds_cover90_pct"] = float(100 * (np.abs(em) <= qs[p]).mean())
                row["ds_mean_std"] = float(sd[m, p].mean())
                res.append(row)
        pm = pd.DataFrame({"tag": tag, "split": sp, "system_id": meta["system_id"].iloc[rows].to_numpy(),
                           "prop": [SHORT[i] for i in meta["prop"].to_numpy()[rows]],
                           "n_comp": meta["n_comp"].to_numpy()[rows], "T_K": meta["T_K"].to_numpy()[rows]})
        pidx = meta["prop"].to_numpy()[rows]
        ar = np.arange(len(rows))
        pm["y_true"], pm["baseline"] = data.y[rows, pidx], data.base_raw[rows, pidx]
        pm["deepsets"], pm["deepsets_std"] = mu[ar, pidx], sd[ar, pidx]
        pm["labelled"] = lab[rows, pidx]
        preds.append(pm)

    mixmask = (meta["n_comp"].to_numpy()[core] >= 2)
    comp_count = np.bincount(data.idx[core][mixmask][data.idx[core][mixmask] >= 0], minlength=len(reg.ik))
    bundle = {"cfg": asdict(cfg), "tasks": tasks, "conformal": qs, "comp_count": comp_count,
              "n_max_train": int(meta["n_comp"].to_numpy()[core].max()), "train_splits": list(train_splits),
              "n_registry": len(reg.ik), "d_static": int(rt.S.shape[1]),
              "states": [{k: v.detach().cpu() for k, v in m.state_dict().items()} for m in models]}
    out_dir.mkdir(parents=True, exist_ok=True)
    torch.save(bundle, out_dir / f"{tag}_bundle.pt")
    if preds:
        pd.concat(preds).to_csv(out_dir / f"{tag}_predictions.csv", index=False)
    return pd.DataFrame(res)


def load_bundle(path, rt, dev):
    z = torch.load(path, map_location="cpu", weights_only=False)
    c = z["cfg"]
    cfg = Cfg(**{**c, "tasks": tuple(c["tasks"]), "lam": tuple(c["lam"])})
    if z["n_registry"] != rt.S.shape[0]:
        raise ValueError("Bundle was trained with a different compounds.csv (registry size mismatch).")
    dummy = {"res_sd": torch.ones(3), "y_mu": torch.zeros(3), "y_sd": torch.ones(3)}
    models = []
    for st in z["states"]:
        m = DeepSetsMixture(rt.S.shape[1], cfg, dummy).to(dev)
        m.load_state_dict(st)
        models.append(m.eval())
    return {"models": models, "cfg": cfg, "tasks": z["tasks"], "conformal": z["conformal"],
            "comp_count": z["comp_count"], "n_max_train": z["n_max_train"]}


# --------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="results/t8")
    ap.add_argument("--pool", default="both", choices=["sum", "attn", "both"])
    ap.add_argument("--mode", default="residual", choices=["residual", "direct"])
    ap.add_argument("--tasks", nargs="+", default=["visc", "dens", "flash"], choices=list(SHORT))
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--bs", type=int, default=512)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--hidden", type=int, default=128)
    ap.add_argument("--embed", type=int, default=128)
    ap.add_argument("--dropout", type=float, default=0.1)
    ap.add_argument("--balance", type=float, default=0.5)
    ap.add_argument("--flash-over-weight", type=float, default=3.0)
    ap.add_argument("--train-splits", nargs="+", default=["S1"])
    ap.add_argument("--ablate", action="store_true")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    a = ap.parse_args()

    root, out, dev = Path(a.root), Path(a.out), torch.device(a.device)
    out.mkdir(parents=True, exist_ok=True)
    train_splits = [s.upper() for s in a.train_splits]
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    rt = RegistryTensors(reg, dev)
    provider = PurePropertyProvider(reg, root / "data/cache/pure_grid.npz")
    data = build_data(root, reg, provider, rt, dev)
    t = to_device(data, dev)

    base = Cfg(pool=a.pool, mode=a.mode, tasks=tuple(SHORT.index(s) for s in a.tasks), seeds=a.seeds,
               epochs=a.epochs, bs=a.bs, lr=a.lr, hidden=a.hidden, embed=a.embed, dropout=a.dropout,
               balance=a.balance, flash_over_weight=a.flash_over_weight)
    if a.ablate:
        runs = []
        for pool, mode in itertools.product(["sum", "attn", "both"], ["residual", "direct"]):
            for name, tk in (("multi", (0, 1, 2)), ("visc", (0,)), ("dens", (1,)), ("flash", (2,))):
                runs.append((f"{name}_{pool}_{mode}", Cfg(**{**asdict(base), "pool": pool, "mode": mode,
                                                                 "tasks": tk, "lam": base.lam})))
    else:
        runs = [("deepsets_" + "-".join(a.tasks) + f"_{a.pool}_{a.mode}", base)]

    allres = []
    for tag, cfg in runs:
        allres.append(run_experiment(data, t, rt, cfg, train_splits, dev, out, tag, reg))
    res = pd.concat(allres, ignore_index=True)
    res.to_csv(out / ("ablation.csv" if a.ablate else "metrics.csv"), index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.4f}".format)
    show = res[res["subset"] == "mix"]
    print("\n=== Baseline vs DeepSets, mixtures only (n_comp>=2) ===")
    print(show.drop(columns=["subset"]).to_string(index=False))
    print("\nReference (progress report, may use different row sets):")
    for k, v in REFERENCE.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
