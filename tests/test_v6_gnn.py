import numpy as np, pandas as pd, torch
from pathlib import Path
from fuelblend.ml import dmpnn as D

R = Path(__file__).resolve().parents[1]

def test_featurize_shapes():
    x, s, d, bf = D.featurize("CC(C)O")
    assert x.shape == (4, D.AF) and len(s) == len(d) == 6 and bf.shape == (6, D.BF)

def test_reverse_edge_pairing():
    _, s, d, _ = D.featurize("c1ccccc1C")
    assert all(s[i] == d[i ^ 1] and d[i] == s[i ^ 1] for i in range(len(s)))

def test_permutation_invariance():
    torch.manual_seed(0)
    net = D.DMPNN(2).eval()
    a = net(D.collate([D.featurize("OCC")])); b = net(D.collate([D.featurize("CCO")]))
    assert torch.allclose(a, b, atol=1e-5)

def test_masked_fit_runs():
    f = [D.featurize(s) for s in ["CCCC", "CCCCC", "c1ccccc1", "CC(C)C"]]
    Y = np.random.default_rng(0).normal(size=(4, 2)); M = np.array([[1, 0]] * 2 + [[0, 1]] * 2, float)
    ms = D.fit_ensemble(f, Y, M, f, Y, M, n_models=2, epochs=2)
    assert D.predict(ms, f)[0].shape == (4, 2)

def test_pool_oxygenates_octane_only():
    p = pd.read_csv(R / "results/v6/gnn/pool.csv")
    assert set(p[~p.is_hc].target) <= {"ron", "mon"}
    assert not p.duplicated(["inchikey", "target"]).any()
