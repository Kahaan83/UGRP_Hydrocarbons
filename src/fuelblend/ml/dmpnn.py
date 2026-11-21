"""Directed MPNN (Chemprop-style) in plain torch, multi-task with masked loss. No PyG.

Track B (v6). Targets are standardised per task; missing labels are masked. Output is the standardised prediction
(callers invert with the stored mean/std).
"""
from __future__ import annotations
import numpy as np
import torch
import torch.nn as nn
from rdkit import Chem

ELEM = [6, 7, 8, 9, 15, 16, 17, 35, 53]
HYB = [Chem.HybridizationType.SP, Chem.HybridizationType.SP2, Chem.HybridizationType.SP3]
BT = [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC]


def _oh(x, choices):
    v = [0.0] * (len(choices) + 1)
    v[choices.index(x) if x in choices else len(choices)] = 1.0
    return v


def atom_feat(a):
    return (_oh(a.GetAtomicNum(), ELEM) + _oh(a.GetTotalDegree(), [1, 2, 3, 4]) + _oh(a.GetTotalNumHs(), [0, 1, 2, 3])
            + _oh(a.GetHybridization(), HYB) + [float(a.GetIsAromatic()), float(a.IsInRing()), a.GetMass() / 100.0])


def bond_feat(b):
    return _oh(b.GetBondType(), BT) + [float(b.GetIsConjugated()), float(b.IsInRing())]


AF = len(atom_feat(Chem.MolFromSmiles("C").GetAtomWithIdx(0)))
BF = len(bond_feat(Chem.MolFromSmiles("CC").GetBondWithIdx(0)))


def featurize(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    x = np.array([atom_feat(a) for a in m.GetAtoms()], dtype=np.float32)
    src, dst, bf = [], [], []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        f = bond_feat(b)
        src += [i, j]; dst += [j, i]; bf += [f, f]
    return x, np.array(src, dtype=np.int64), np.array(dst, dtype=np.int64), np.array(bf, dtype=np.float32).reshape(-1, BF)


def collate(feats):
    xs, srcs, dsts, bfs, ab = [], [], [], [], []
    na = ne = 0
    for g, (x, s, d, bf) in enumerate(feats):
        xs.append(x); srcs.append(s + na); dsts.append(d + na); bfs.append(bf)
        ab.append(np.full(len(x), g)); na += len(x); ne += len(s)
    return (torch.from_numpy(np.concatenate(xs)), torch.from_numpy(np.concatenate(srcs)), torch.from_numpy(np.concatenate(dsts)),
            torch.from_numpy(np.concatenate(bfs)), torch.from_numpy(np.concatenate(ab)), len(feats))


class DMPNN(nn.Module):
    def __init__(self, n_tasks: int, hidden: int = 200, depth: int = 4, dropout: float = 0.1, ffn: int = 200):
        super().__init__()
        self.depth = depth
        self.W_i = nn.Linear(AF + BF, hidden, bias=False)
        self.W_h = nn.Linear(hidden, hidden, bias=False)
        self.W_o = nn.Linear(AF + hidden, hidden)
        self.drop = nn.Dropout(dropout)
        self.head = nn.Sequential(nn.Linear(hidden, ffn), nn.ReLU(), nn.Dropout(dropout), nn.Linear(ffn, n_tasks))

    def encode(self, batch):
        x, src, dst, bf, ab, n_g = batch
        n_a, n_e = x.shape[0], src.shape[0]
        if n_e == 0:
            return torch.zeros(n_g, self.W_o.out_features)
        # edge e: src[e]->dst[e]; reverse of e is e^1 (edges were added in pairs)
        h0 = self.W_i(torch.cat([x[src], bf], 1))
        h = torch.relu(h0)
        rev = torch.arange(n_e) ^ 1
        for _ in range(self.depth - 1):
            a_in = torch.zeros(n_a, h.shape[1]).index_add_(0, dst, h)   # sum of incoming messages per atom
            m = a_in[src] - h[rev]                                       # exclude the reverse edge
            h = self.drop(torch.relu(h0 + self.W_h(m)))
        a_in = torch.zeros(n_a, h.shape[1]).index_add_(0, dst, h)
        a_h = self.drop(torch.relu(self.W_o(torch.cat([x, a_in], 1))))
        out = torch.zeros(n_g, a_h.shape[1]).index_add_(0, ab, a_h)
        cnt = torch.zeros(n_g).index_add_(0, ab, torch.ones(n_a)).clamp(min=1).unsqueeze(1)
        return out / cnt

    def forward(self, batch):
        return self.head(self.encode(batch))


def fit_ensemble(feats, Y, M, val_feats, val_Y, val_M, n_models=5, epochs=150, patience=25, lr=1e-3, wd=1e-5, bs=32,
                 hidden=200, depth=4, dropout=0.1, seed=0):
    """Train n_models with different seeds; early-stop each on masked val MSE. Y,M are standardised/ masks (N,T)."""
    models = []
    Yt, Mt = torch.tensor(Y, dtype=torch.float32), torch.tensor(M, dtype=torch.float32)
    vY, vM = torch.tensor(val_Y, dtype=torch.float32), torch.tensor(val_M, dtype=torch.float32)
    vb = collate(val_feats)
    n = len(feats)
    for k in range(n_models):
        torch.manual_seed(seed * 1000 + k)
        rng = np.random.default_rng(seed * 1000 + k)
        net = DMPNN(Y.shape[1], hidden, depth, dropout)
        opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
        best, best_sd, bad = 1e9, None, 0
        for ep in range(epochs):
            net.train()
            perm = rng.permutation(n)
            for i in range(0, n, bs):
                idx = perm[i:i + bs]
                if len(idx) < 2: continue
                out = net(collate([feats[j] for j in idx]))
                loss = (((out - Yt[idx]) ** 2) * Mt[idx]).sum() / Mt[idx].sum().clamp(min=1)
                opt.zero_grad(); loss.backward(); opt.step()
            net.eval()
            with torch.no_grad():
                vo = net(vb)
                vl = float((((vo - vY) ** 2) * vM).sum() / vM.sum().clamp(min=1))
            if vl < best - 1e-5: best, bad, best_sd = vl, 0, {a: b.clone() for a, b in net.state_dict().items()}
            else:
                bad += 1
                if bad >= patience: break
        net.load_state_dict(best_sd); net.eval(); models.append(net)
    return models


def predict(models, feats):
    b = collate(feats)
    with torch.no_grad():
        P = torch.stack([m(b) for m in models]).numpy()
    return P.mean(0), P.std(0)
