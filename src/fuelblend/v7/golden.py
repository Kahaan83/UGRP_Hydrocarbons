"""Deterministic CPU fixtures whose outputs are frozen in tests/data/v7_golden_*.npy (before the v7 hooks were added)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[3]


def dmpnn_run(**kw) -> np.ndarray:
    from fuelblend.ml import dmpnn as D
    pool = pd.read_csv(ROOT / "results/v6/gnn/pool.csv")
    comp = pool.drop_duplicates("inchikey").sort_values("inchikey").head(30).reset_index(drop=True)
    tasks = ["flash_C", "ysi"]
    Y = np.full((30, 2), np.nan)
    for j, t in enumerate(tasks):
        v = pool[pool.target == t].groupby("inchikey").value.first()
        Y[:, j] = comp.inchikey.map(v).values
    M = ~np.isnan(Y)
    mu = np.nanmean(Y, 0); sd = np.nanstd(Y, 0) + 1e-9
    Z = np.nan_to_num((Y - mu) / sd)
    feats = [D.featurize(s) for s in comp.smiles]
    tr, va = np.arange(22), np.arange(22, 30)
    ms = D.fit_ensemble([feats[i] for i in tr], Z[tr], M[tr], [feats[i] for i in va], Z[va], M[va],
                        n_models=2, epochs=6, patience=5, seed=0, **kw)
    return D.predict(ms, [feats[i] for i in va])[0]


def molpool_run(**kw) -> np.ndarray:
    from fuelblend.ml import molpool as P
    g = np.random.default_rng(0)
    H = g.normal(size=(30, 3, 12)).astype(np.float32)
    X = g.dirichlet(np.ones(3), size=30).astype(np.float32)
    r = g.normal(size=30).astype(np.float32)
    return P.fit_predict(H[:22], X[:22], r[:22], H[22:], X[22:], epochs=5, seeds=(0, 1), **kw)


def set_models_run(**kw) -> np.ndarray:
    from fuelblend.ml import set_models as sm
    g = np.random.default_rng(0)
    n, k = 30, 3
    F = g.normal(size=(n, k, sm.N_COMP_FEATS)).astype(np.float32)
    X = g.dirichlet(np.ones(k), size=n).astype(np.float32)
    b = sm.SetBatch.from_numpy(F, X, np.ones((n, k), bool), g.normal(size=n), g.integers(0, 3, size=(n, k)), np.zeros(n, int))
    y = torch.randn(n, generator=torch.Generator().manual_seed(0))
    w = torch.ones(n) / n
    ens = sm.train_ensemble(lambda: sm.PairDeepSets(width=8), b, y, w, seeds=[0, 1], epochs=4, batch_size=16, **kw)
    return sm.predict_ensemble(ens, b)[0]
