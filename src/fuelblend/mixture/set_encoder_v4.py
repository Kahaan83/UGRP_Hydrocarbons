"""v4 Task 16: permutation-invariant set encoder that outputs RESIDUALS on top of physics baselines.

    prediction_head = baseline_head + scale_head * head(z, T)

* ``component_featurizer`` maps static (temperature-independent) per-component descriptors [B, N, d] -> [B, N, h].
* ``pool``: ``"sum"`` -> z = sum_i w_i phi_i;  ``"attention"`` -> z = softmax_i(s_i + log w_i) . v_i (+ the
  weighted sum).  Both ignore components with w_i == 0 *exactly* (mask / log-weight -inf), and both are invariant to
  the order of the components.
* Each ``Head`` is a small MLP whose final linear layer is ZERO-INITIALISED, so an untrained model returns the
  physics baseline exactly.
* Temperature enters ONLY through the Task 11 structured form: a "structured" head emits (a0, a1) from the mixture
  state z and the correction is  a0 + a1 * (1000/T - T_REF)  -- exactly linear in 1000/T for every z, so
  extrapolation to -40 C cannot develop curvature (same guarantee as ``StructuredResidual``).  A "scalar" head
  (e.g. flash point, no T) emits a0 only.
* Split helpers (``dev_group_folds``, ``assert_no_locked``) live here so that the "locked fuels never enter
  training / calibration" contract is unit-testable without data.
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn

T_REF = 3.7          # 1000/T centring constant (~270 K); a0 is the residual at 1000/T = T_REF


class MLPFeaturizer(nn.Module):
    def __init__(self, d_in: int, hidden: int = 64, depth: int = 2):
        super().__init__()
        layers, d = [], d_in
        for _ in range(depth):
            layers += [nn.Linear(d, hidden), nn.SiLU()]
            d = hidden
        self.net = nn.Sequential(*layers)
        self.d_out = hidden

    def forward(self, x):                      # [B, N, d_in] -> [B, N, hidden]
        return self.net(x)


class Head(nn.Module):
    """Residual head.  structured=True -> (a0, a1) with delta = a0 + a1*(1000/T - T_REF); else delta = a0."""

    def __init__(self, d_in: int, structured: bool, scale: float, hidden: int = 64):
        super().__init__()
        self.structured, self.scale = bool(structured), float(scale)
        self.body = nn.Sequential(nn.Linear(d_in, hidden), nn.SiLU())
        self.out = nn.Linear(hidden, 2 if structured else 1)
        nn.init.zeros_(self.out.weight)
        nn.init.zeros_(self.out.bias)          # zero residual at init -> prediction == physics

    def coefficients(self, z):
        o = self.out(self.body(z))
        return o[:, 0], (o[:, 1] if self.structured else torch.zeros_like(o[:, 0]))

    def forward(self, z, T_K):
        a0, a1 = self.coefficients(z)
        t = 1000.0 / T_K - T_REF
        return self.scale * (a0 + a1 * t)


class SetEncoder(nn.Module):
    def __init__(self, component_featurizer: nn.Module, pool: str = "attention", heads: dict | None = None,
                 d_context: int = 0):
        super().__init__()
        if pool not in ("attention", "sum"):
            raise ValueError("pool must be 'attention' or 'sum'")
        if not heads:
            raise ValueError("at least one head is required")
        self.phi, self.pool = component_featurizer, pool
        h = component_featurizer.d_out
        self.d_context = int(d_context)
        if pool == "attention":
            self.score = nn.Linear(h, 1)
            self.value = nn.Linear(h, h)
        self.heads = nn.ModuleDict(heads)
        for hd in heads.values():
            assert isinstance(hd, Head)
        self.d_z = (2 * h if pool == "attention" else h)
        # heads were built with d_in = d_z + d_context by the caller (see make_encoder)

    def mixture_state(self, comp, w, ctx=None):
        """comp [B,N,d], w [B,N] (>=0, rows need not sum to 1; padded slots have w=0) -> z [B, d_z + d_context]."""
        w = w.clamp_min(0.0)
        phi = self.phi(comp)                                       # [B,N,h]
        wsum = (w.unsqueeze(-1) * phi).sum(1)                      # [B,h]
        if self.pool == "sum":
            z = wsum
        else:
            logit = self.score(phi).squeeze(-1)                    # [B,N]
            mask = w > 0
            logit = torch.where(mask, logit + torch.log(w.clamp_min(1e-30)), torch.full_like(logit, float("-inf")))
            att = torch.softmax(logit, dim=1)
            att = torch.nan_to_num(att, nan=0.0)                   # all-zero-weight row -> zero state
            z = torch.cat([wsum, (att.unsqueeze(-1) * self.value(phi)).sum(1)], dim=1)
        if self.d_context:
            z = torch.cat([z, ctx], dim=1)
        return z

    def residual(self, head: str, comp, w, T_K, ctx=None):
        return self.heads[head](self.mixture_state(comp, w, ctx), T_K)

    def forward(self, head: str, comp, w, T_K, baseline, ctx=None):
        """prediction = baseline + residual (baseline in the head's own units: ln nu for viscosity, degC for flash)."""
        return baseline + self.residual(head, comp, w, T_K, ctx)


def make_encoder(d_comp: int, pool: str = "attention", hidden: int = 64, d_context: int = 0,
                 head_specs: dict | None = None) -> SetEncoder:
    """head_specs: {name: dict(structured=bool, scale=float)}.  Default: viscosity (structured, 0.3 ln-units) and
    flash (scalar, 5 degC)."""
    head_specs = head_specs or {"viscosity": dict(structured=True, scale=0.3), "flash": dict(structured=False, scale=5.0)}
    phi = MLPFeaturizer(d_comp, hidden)
    d_z = (2 * hidden if pool == "attention" else hidden) + d_context
    heads = {k: Head(d_z, hidden=hidden, **v) for k, v in head_specs.items()}
    return SetEncoder(phi, pool, heads, d_context=d_context)


# ----------------------------------------------------------------------------------- split contracts
def dev_group_folds(group_ids, locked, n_folds: int = 9, seed: int = 0):
    """Deterministic grouped K-fold over the NON-locked groups only.  Returns list of (train_pos, test_pos) index
    arrays into ``group_ids``; locked-group rows appear in neither."""
    g = np.asarray(group_ids).astype(str)
    lock = set(map(str, locked))
    dev_groups = np.array(sorted(set(g) - lock))
    k = min(n_folds, len(dev_groups))
    if k < 2:
        raise ValueError("need at least 2 dev groups")
    perm = np.random.default_rng(seed).permutation(len(dev_groups))
    parts = np.array_split(dev_groups[perm], k)
    is_dev = ~np.isin(g, list(lock))
    out = []
    for p in parts:
        te = np.where(np.isin(g, p))[0]
        tr = np.where(is_dev & ~np.isin(g, p))[0]
        out.append((tr, te))
    return out


def assert_no_locked(index_sets: dict, group_ids, locked):
    """Raise if any index set (train / early-stop / calibration / selection) contains a locked group's row."""
    g = np.asarray(group_ids).astype(str)
    lock = set(map(str, locked))
    for name, idx in index_sets.items():
        bad = [x for x in g[np.asarray(idx, dtype=int)] if x in lock]
        if bad:
            raise AssertionError(f"locked group(s) {sorted(set(bad))[:5]} present in '{name}' indices")
    return True
