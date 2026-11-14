"""WS4 neural mixture encoders (plain PyTorch), all residual corrections on a physics value (delta_ws4_mixture_ml).

Models map a padded set of components (features F, mole fractions X, mask M, row-level 1000/T, molecule index, property id) to a scalar
correction per row, in the training-scaled target unit:
  PairDeepSets : sum_{i<j} x_i x_j phi(h_i + h_j, |h_i - h_j|, 1000/T),  h = MLP(F_i)
  MolSetsGNN   : (1 - sum x^2) * head( sum_i a_i z_i ),  a = softmax(w.z_i + ln x_i),  z_i = MLP(graph_embedding_i, F_i, 1000/T)
  SetTransformer: (1 - sum x^2) * head_p( PMA(SAB(SAB(tokens))) ),  token_i = MLP(F_i, x_i, 1000/T), one head per property p
All three are permutation invariant, zero for a pure component and unaffected by masked padding slots.
Research code for a declared comparison; nothing here is used by predict_fuel.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import torch
from torch import nn

N_COMP_FEATS = 14          # 6 static + ln_eta + v_molar + 5 class one-hot + ln_eta_missing
N_ATOM_FEATS = 13


# ------------------------------------------------------------------------------------------------ data containers
@dataclass
class SetBatch:
    F: torch.Tensor      # (N, K, D) float
    X: torch.Tensor      # (N, K) float, 0 at padding
    M: torch.Tensor      # (N, K) bool, True = live component
    invT: torch.Tensor   # (N,) float (standardised 1000/T)
    mol: torch.Tensor    # (N, K) long, index into a GraphBank (0 at padding)
    prop: torch.Tensor   # (N,) long, 0 viscosity / 1 density

    @staticmethod
    def from_numpy(F, X, M, invT, mol, prop) -> "SetBatch":
        return SetBatch(torch.as_tensor(np.asarray(F, np.float32)), torch.as_tensor(np.asarray(X, np.float32)),
                        torch.as_tensor(np.asarray(M, bool)), torch.as_tensor(np.asarray(invT, np.float32)),
                        torch.as_tensor(np.asarray(mol, np.int64)), torch.as_tensor(np.asarray(prop, np.int64)))

    def __len__(self):
        return self.X.shape[0]

    def take(self, idx) -> "SetBatch":
        idx = torch.as_tensor(idx, dtype=torch.long)
        return SetBatch(self.F[idx], self.X[idx], self.M[idx], self.invT[idx], self.mol[idx], self.prop[idx])

    def to(self, dev) -> "SetBatch":
        return SetBatch(*(t.to(dev) for t in (self.F, self.X, self.M, self.invT, self.mol, self.prop)))


class GraphBank:
    """All unique molecules as one disjoint graph: atom features, directed edges (both directions), atom -> molecule index."""

    def __init__(self, atom_feat, edge_index, atom2mol, n_mol):
        self.atom_feat, self.edge_index, self.atom2mol, self.n_mol = atom_feat, edge_index, atom2mol, int(n_mol)

    @staticmethod
    def atom_features(a) -> list:
        deg = min(a.GetDegree(), 4)
        nh = min(a.GetTotalNumHs(), 3)
        return ([1.0 if a.GetAtomicNum() == 6 else 0.0, 0.0 if a.GetAtomicNum() == 6 else 1.0]
                + [1.0 if deg == k else 0.0 for k in range(5)] + [1.0 if nh == k else 0.0 for k in range(4)]
                + [float(a.GetIsAromatic()), float(a.IsInRing())])

    @classmethod
    def from_mols(cls, mols) -> "GraphBank":
        feats, src, dst, a2m, off = [], [], [], [], 0
        for k, m in enumerate(mols):
            if m is None:
                feats.append([0.0] * N_ATOM_FEATS)      # placeholder single atom for an unparsable molecule
                a2m.append(k)
                off += 1
                continue
            for a in m.GetAtoms():
                feats.append(cls.atom_features(a))
                a2m.append(k)
            for b in m.GetBonds():
                i, j = b.GetBeginAtomIdx() + off, b.GetEndAtomIdx() + off
                src += [i, j]
                dst += [j, i]
            off += m.GetNumAtoms()
        return cls(torch.tensor(feats, dtype=torch.float32), torch.tensor([src, dst], dtype=torch.long).reshape(2, -1),
                   torch.tensor(a2m, dtype=torch.long), len(mols))

    @classmethod
    def from_smiles(cls, smiles) -> "GraphBank":
        from rdkit import Chem
        return cls.from_mols([Chem.MolFromSmiles(s) for s in smiles])

    def to(self, dev) -> "GraphBank":
        return GraphBank(self.atom_feat.to(dev), self.edge_index.to(dev), self.atom2mol.to(dev), self.n_mol)


# ------------------------------------------------------------------------------------------------ models
def mlp(i, h, o, depth=2):
    layers, d = [], i
    for _ in range(depth - 1):
        layers += [nn.Linear(d, h), nn.SiLU()]
        d = h
    layers.append(nn.Linear(d, o))
    return nn.Sequential(*layers)


def _nonideal_factor(b: SetBatch) -> torch.Tensor:
    x = b.X * b.M
    return 1.0 - (x ** 2).sum(1)


class PairDeepSets(nn.Module):
    def __init__(self, width=64, d_in=N_COMP_FEATS):
        super().__init__()
        self.enc = mlp(d_in, width, width)
        self.phi = mlp(2 * width + 1, width, 1)

    def forward(self, b: SetBatch) -> torch.Tensor:
        h = self.enc(b.F)
        K = h.shape[1]
        ii, jj = torch.triu_indices(K, K, offset=1)
        hi, hj = h[:, ii], h[:, jj]
        t = b.invT[:, None, None].expand(-1, len(ii), 1)
        g = self.phi(torch.cat([hi + hj, (hi - hj).abs(), t], -1)).squeeze(-1)
        xx = (b.X * b.M)[:, ii] * (b.X * b.M)[:, jj]
        return (xx * g).sum(1)


class MolSetsGNN(nn.Module):
    def __init__(self, bank: GraphBank, width=64, layers=3, d_in=N_COMP_FEATS):
        super().__init__()
        self.bank = bank
        self.inp = nn.Linear(N_ATOM_FEATS, width)
        self.self_l = nn.ModuleList(nn.Linear(width, width) for _ in range(layers))
        self.nb_l = nn.ModuleList(nn.Linear(width, width) for _ in range(layers))
        self.comp = mlp(width + d_in + 1, width, width)
        self.att = nn.Linear(width, 1)
        self.head = mlp(width, width, 1)

    def set_bank(self, bank: GraphBank):
        self.bank = bank

    def mol_embeddings(self) -> torch.Tensor:
        bk = self.bank
        h = torch.relu(self.inp(bk.atom_feat))
        src, dst = bk.edge_index
        deg = torch.zeros(h.shape[0], device=h.device).index_add_(0, dst, torch.ones_like(dst, dtype=h.dtype)).clamp(min=1.0)
        for ls, ln in zip(self.self_l, self.nb_l):
            agg = torch.zeros_like(h).index_add_(0, dst, h[src]) / deg[:, None]
            h = torch.relu(ls(h) + ln(agg))
        return torch.zeros(bk.n_mol, h.shape[1], device=h.device).index_add_(0, bk.atom2mol, h)

    def forward(self, b: SetBatch) -> torch.Tensor:
        e = self.mol_embeddings()[b.mol]
        t = b.invT[:, None, None].expand(-1, b.F.shape[1], 1)
        z = self.comp(torch.cat([e, b.F, t], -1))
        logit = self.att(z).squeeze(-1) + torch.log(b.X.clamp(min=1e-12))
        logit = logit.masked_fill(~b.M, -1e9)
        a = torch.softmax(logit, 1)
        s = self.head((a[..., None] * z).sum(1)).squeeze(-1)
        return _nonideal_factor(b) * s


class _SAB(nn.Module):
    def __init__(self, width, heads):
        super().__init__()
        self.att = nn.MultiheadAttention(width, heads, batch_first=True)
        self.n1, self.n2 = nn.LayerNorm(width), nn.LayerNorm(width)
        self.ff = mlp(width, width, width)

    def forward(self, h, pad):
        a, _ = self.att(h, h, h, key_padding_mask=pad, need_weights=False)
        h = self.n1(h + a)
        return self.n2(h + self.ff(h))


class SetTransformer(nn.Module):
    def __init__(self, width=64, heads=4, blocks=2, n_props=2, d_in=N_COMP_FEATS):
        super().__init__()
        self.tok = mlp(d_in + 2, width, width)
        self.sabs = nn.ModuleList(_SAB(width, heads) for _ in range(blocks))
        self.seed = nn.Parameter(torch.randn(1, 1, width) * 0.1)
        self.pma = nn.MultiheadAttention(width, heads, batch_first=True)
        self.heads = nn.ModuleList(mlp(width, width, 1) for _ in range(n_props))

    def forward(self, b: SetBatch) -> torch.Tensor:
        t = b.invT[:, None, None].expand(-1, b.F.shape[1], 1)
        h = self.tok(torch.cat([b.F, b.X[..., None], t], -1))
        pad = ~b.M
        for s in self.sabs:
            h = s(h, pad)
        q = self.seed.expand(h.shape[0], -1, -1)
        p, _ = self.pma(q, h, h, key_padding_mask=pad, need_weights=False)
        p = p.squeeze(1)
        out = torch.stack([hd(p).squeeze(-1) for hd in self.heads], 1)
        s = out.gather(1, b.prop.clamp(max=len(self.heads) - 1)[:, None]).squeeze(1)
        return _nonideal_factor(b) * s


# ------------------------------------------------------------------------------------------------ training
def system_macro_weights(systems, props=None) -> np.ndarray:
    """Row weight 1/(rows of its system)/(systems of its property): each property's weights sum to 1."""
    systems = np.asarray(systems)
    props = np.zeros(len(systems), int) if props is None else np.asarray(props)
    w = np.zeros(len(systems))
    for p in np.unique(props):
        m = props == p
        s, inv, cnt = np.unique(systems[m], return_inverse=True, return_counts=True)
        w[m] = 1.0 / cnt[inv] / len(s)
    return w


def train_ensemble(make_model, data: SetBatch, y, w, seeds=(0, 1, 2, 3, 4), epochs=60, batch_size=512, lr=2e-3, wd=1e-4,
                   device="cpu", clip=1.0, log=None):
    """Weighted-MSE training (per-batch sum(w e^2)/sum(w)), AdamW + cosine schedule, one model per seed. y is target-scaled."""
    y = torch.as_tensor(np.asarray(y, np.float32)).to(device)
    w = torch.as_tensor(np.asarray(w, np.float32)).to(device)
    data = data.to(device)
    n = len(y)
    models = []
    for s in seeds:
        torch.manual_seed(int(s))
        g = torch.Generator().manual_seed(int(s))
        m = make_model().to(device)
        if hasattr(m, "bank"):
            m.set_bank(m.bank.to(device))
        opt = torch.optim.AdamW(m.parameters(), lr=lr, weight_decay=wd)
        steps = epochs * max(1, math.ceil(n / batch_size))
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps)
        m.train()
        for ep in range(epochs):
            perm = torch.randperm(n, generator=g).to(device)
            tot = 0.0
            for k in range(0, n, batch_size):
                idx = perm[k:k + batch_size]
                b = data.take(idx)
                e = m(b) - y[idx]
                ww = w[idx]
                loss = (ww * e ** 2).sum() / ww.sum().clamp(min=1e-12)
                opt.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(m.parameters(), clip)
                opt.step()
                sch.step()
                tot += float(loss.detach()) * len(idx)
            if log is not None and (ep + 1) % max(1, epochs // 4) == 0:
                log(f"  seed {s} epoch {ep + 1}/{epochs} loss {tot / n:.4g}")
        m.eval()
        models.append(m)
    return models


@torch.no_grad()
def predict_ensemble(models, data: SetBatch, device="cpu", chunk=4096):
    """-> (mean, std) over the ensemble members, numpy, in the training-scaled unit."""
    data = data.to(device)
    out = []
    for m in models:
        parts = [m(data.take(torch.arange(k, min(k + chunk, len(data))))).cpu().numpy() for k in range(0, len(data), chunk)]
        out.append(np.concatenate(parts) if parts else np.zeros(0))
    P = np.vstack(out)
    return P.mean(0), P.std(0)


# ------------------------------------------------------------------------------------------------ evaluation helpers
def gate_pass(boot: dict) -> bool:
    """Gain = MAE_pair - MAE_candidate; pass only if the whole interval is above 0 (a CI below 0 is a significant loss)."""
    return bool(boot["lo"] > 0)


def fold_test_systems(systems, dev_rows, n_splits=5) -> list:
    """Test-system lists of GroupKFold(n_splits).split(dev_rows, groups=systems[dev_rows]) -- the baseline scripts' exact call."""
    from sklearn.model_selection import GroupKFold
    systems = np.asarray(systems)
    return [sorted(set(systems[dev_rows[te]])) for _, te in GroupKFold(n_splits=n_splits).split(dev_rows, groups=systems[dev_rows])]
