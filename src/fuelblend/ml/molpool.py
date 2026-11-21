"""Track C (v6): mole-fraction pooled mixture head on frozen D-MPNN component embeddings (MolPool-style) with a pair term.

Per-component embedding h_i (frozen Track B encoder) -> phi_i = MLP(h_i); mixture vector z = sum_i x_i phi_i;
pair term q = 0.5 * (|u|^2 - sum_i x_i^2 |P phi_i|^2) with u = sum_i x_i P phi_i, i.e. sum_{i<j} x_i x_j (P phi_i).(P phi_j);
output = MLP([z, q, aux]) = residual over a physics baseline. The last layer is zero-initialised so an untrained head returns 0
(the baseline). Plain torch, CPU.
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn

from fuelblend.ml import dmpnn as D


def load_encoder(path):
    ck = torch.load(path, map_location="cpu", weights_only=False)
    c = ck["cfg"]
    models = []
    for sd in ck["states"]:
        m = D.DMPNN(len(ck["tasks"]), c["hidden"], c["depth"], c["dropout"])
        m.load_state_dict(sd)
        m.eval()
        models.append(m)
    return models


def embed(models, smiles) -> np.ndarray:
    """Ensemble-mean frozen embedding, one row per SMILES."""
    b = D.collate([D.featurize(s) for s in smiles])
    with torch.no_grad():
        return torch.stack([m.encode(b) for m in models]).mean(0).numpy()


def pad_sets(emb_rows, x_rows):
    """emb_rows: list of (k_i, E) arrays; x_rows: list of (k_i,) arrays -> padded (B,K,E), (B,K)."""
    K = max(len(x) for x in x_rows)
    E = emb_rows[0].shape[1]
    H = np.zeros((len(x_rows), K, E), np.float32)
    X = np.zeros((len(x_rows), K), np.float32)
    for i, (h, x) in enumerate(zip(emb_rows, x_rows)):
        H[i, :len(x)] = h
        X[i, :len(x)] = x
    return H, X


class MolPool(nn.Module):
    def __init__(self, emb_dim: int, d: int = 16, rank: int = 8, n_aux: int = 0, hidden: int = 32, dropout: float = 0.1):
        super().__init__()
        self.phi = nn.Sequential(nn.Linear(emb_dim, d), nn.ReLU(), nn.Dropout(dropout), nn.Linear(d, d))
        self.P = nn.Linear(d, rank, bias=False)
        self.out = nn.Sequential(nn.Linear(d + 1 + n_aux, hidden), nn.ReLU(), nn.Linear(hidden, 1))
        nn.init.zeros_(self.out[-1].weight)
        nn.init.zeros_(self.out[-1].bias)

    def forward(self, H, X, aux=None):
        phi = self.phi(H)                                     # (B,K,d)
        z = (X.unsqueeze(-1) * phi).sum(1)                    # mole-fraction pooling
        p = self.P(phi)                                       # (B,K,r)
        u = (X.unsqueeze(-1) * p).sum(1)
        q = 0.5 * ((u ** 2).sum(-1) - ((X ** 2).unsqueeze(-1) * p ** 2).sum((1, 2)))
        f = [z, q.unsqueeze(-1)] + ([aux] if aux is not None else [])
        return self.out(torch.cat(f, -1)).squeeze(-1)


def fit_predict(H_tr, X_tr, r_tr, H_te, X_te, aux_tr=None, aux_te=None, d=16, wd=1e-2, epochs=200, lr=3e-3, seeds=(0,),
                rank=8, dropout=0.1):
    """Train on residual targets r_tr (already standardised); returns the mean prediction over seeds for the test sets."""
    Ht, Xt = torch.tensor(H_tr), torch.tensor(X_tr)
    rt = torch.tensor(r_tr, dtype=torch.float32)
    He, Xe = torch.tensor(H_te), torch.tensor(X_te)
    at = None if aux_tr is None else torch.tensor(aux_tr, dtype=torch.float32)
    ae = None if aux_te is None else torch.tensor(aux_te, dtype=torch.float32)
    preds = []
    for s in seeds:
        torch.manual_seed(1000 + s)
        net = MolPool(H_tr.shape[2], d, rank, 0 if aux_tr is None else aux_tr.shape[1], dropout=dropout)
        opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
        n = len(rt)
        g = torch.Generator().manual_seed(s)
        for _ in range(epochs):
            net.train()
            perm = torch.randperm(n, generator=g)
            for i in range(0, n, 64):
                idx = perm[i:i + 64]
                if len(idx) < 2:
                    continue
                loss = ((net(Ht[idx], Xt[idx], None if at is None else at[idx]) - rt[idx]) ** 2).mean()
                opt.zero_grad()
                loss.backward()
                opt.step()
        net.eval()
        with torch.no_grad():
            preds.append(net(He, Xe, ae).numpy())
    return np.mean(preds, 0)
