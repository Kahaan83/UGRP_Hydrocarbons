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
        m.to(D.DEVICE).eval()
        models.append(m)
    return models


def embed(models, smiles) -> np.ndarray:
    """Ensemble-mean frozen embedding, one row per SMILES."""
    b = D.collate([D.featurize(s) for s in smiles])
    with torch.no_grad():
        return torch.stack([m.encode(b) for m in models]).mean(0).cpu().numpy()


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
                rank=8, dropout=0.1, hidden=32, device=None, val_idx=None, patience=None, logger=None, on_epoch=None):
    """Train on residual targets r_tr (already standardised); returns the mean prediction over seeds for the test sets.

    v7 hooks (defaults reproduce v6 exactly): device (None = cpu); val_idx = indices into the training arrays held out for
    early stopping on val MSE (patience, default 50, best state restored; they are not trained on); logger; on_epoch(epoch,
    val_loss or train_loss)."""
    import time
    dev = torch.device("cpu" if device is None else device)
    if dev.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(f"CUDA requested but torch {torch.__version__} has no CUDA")
    Ht, Xt = torch.tensor(H_tr), torch.tensor(X_tr)
    rt = torch.tensor(r_tr, dtype=torch.float32)
    He, Xe = torch.tensor(H_te), torch.tensor(X_te)
    at = None if aux_tr is None else torch.tensor(aux_tr, dtype=torch.float32)
    ae = None if aux_te is None else torch.tensor(aux_te, dtype=torch.float32)
    vi = None
    if val_idx is not None:
        vmask = np.zeros(len(rt), bool)
        vmask[np.asarray(val_idx)] = True
        vi = torch.as_tensor(np.where(vmask)[0])
        ti = torch.as_tensor(np.where(~vmask)[0])
        Hv, Xv, rv, av = Ht[vi], Xt[vi], rt[vi], (None if at is None else at[vi])
        Ht, Xt, rt, at = Ht[ti], Xt[ti], rt[ti], (None if at is None else at[ti])
        Hv, Xv, rv, av = Hv.to(dev), Xv.to(dev), rv.to(dev), (None if av is None else av.to(dev))
    pat = 50 if patience is None else patience
    Ht, Xt, rt, He, Xe = Ht.to(dev), Xt.to(dev), rt.to(dev), He.to(dev), Xe.to(dev)
    at = None if at is None else at.to(dev)
    ae = None if ae is None else ae.to(dev)
    preds = []
    for s in seeds:
        torch.manual_seed(1000 + s)
        net = MolPool(H_tr.shape[2], d, rank, 0 if aux_tr is None else aux_tr.shape[1], hidden=hidden, dropout=dropout).to(dev)
        opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=wd)
        n = len(rt)
        g = torch.Generator().manual_seed(s)
        lg = None if logger is None else logger.child(seed=s)
        best, best_sd, bad = float("inf"), None, 0
        for ep in range(epochs):
            t0 = time.time()
            net.train()
            perm = torch.randperm(n, generator=g).to(dev)
            tl, tn = 0.0, 0
            for i in range(0, n, 64):
                idx = perm[i:i + 64]
                if len(idx) < 2:
                    continue
                loss = ((net(Ht[idx], Xt[idx], None if at is None else at[idx]) - rt[idx]) ** 2).mean()
                opt.zero_grad()
                loss.backward()
                opt.step()
                if lg is not None:
                    tl += float(loss.detach())
                    tn += 1
            vl = None
            if vi is not None:
                net.eval()
                with torch.no_grad():
                    vl = float(((net(Hv, Xv, av) - rv) ** 2).mean())
            if lg is not None:
                lg.log(ep + 1, train_loss=tl / max(tn, 1), **({} if vl is None else {"val_loss": vl}),
                       lr=opt.param_groups[0]["lr"], sec=time.time() - t0)
            if on_epoch is not None:
                on_epoch(ep + 1, vl if vl is not None else tl / max(tn, 1))
            if vl is not None:
                if vl < best - 1e-6:
                    best, bad, best_sd = vl, 0, {a: b.clone() for a, b in net.state_dict().items()}
                else:
                    bad += 1
                    if bad >= pat:
                        break
        if best_sd is not None:
            net.load_state_dict(best_sd)
        net.eval()
        with torch.no_grad():
            preds.append(net(He, Xe, ae).cpu().numpy())
    return np.mean(preds, 0)
