"""WS4 set models: structural properties (permutation, pure-component zero, padding, determinism, graph featuriser) and gate helpers."""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from fuelblend.ml import set_models as sm  # noqa: E402


def _batch(n=6, k=3, d=sm.N_COMP_FEATS, seed=0, live=None):
    g = np.random.default_rng(seed)
    F = g.normal(size=(n, k, d)).astype(np.float32)
    X = g.dirichlet(np.ones(k), size=n).astype(np.float32)
    M = np.ones((n, k), bool)
    if live is not None:
        M[:, live:] = False
        X[:, live:] = 0
        X = X / X.sum(1, keepdims=True)
    mol = g.integers(0, 3, size=(n, k))
    return sm.SetBatch.from_numpy(F, X, M, g.normal(size=n).astype(np.float32), mol, np.zeros(n, int))


def _bank():
    return sm.GraphBank.from_smiles(["CCC", "Cc1ccccc1", "CCCCCCCC"])


def _models(bank):
    torch.manual_seed(0)
    return {"pds": sm.PairDeepSets(width=16), "gnn": sm.MolSetsGNN(bank, width=16), "st": sm.SetTransformer(width=16, n_props=2)}


@pytest.mark.parametrize("name", ["pds", "gnn", "st"])
def test_permutation_invariance(name):
    m = _models(_bank())[name].eval()
    b = _batch()
    perm = [2, 0, 1]
    bp = sm.SetBatch(b.F[:, perm], b.X[:, perm], b.M[:, perm], b.invT, b.mol[:, perm], b.prop)
    with torch.no_grad():
        assert torch.allclose(m(b), m(bp), atol=1e-5)


@pytest.mark.parametrize("name", ["pds", "gnn", "st"])
def test_zero_for_pure_component(name):
    m = _models(_bank())[name].eval()
    b = _batch(live=1)                       # one live component with x = 1, other slots masked
    with torch.no_grad():
        assert torch.allclose(m(b), torch.zeros(len(b.X)), atol=1e-6)


@pytest.mark.parametrize("name", ["pds", "gnn", "st"])
def test_padding_invariance(name):
    m = _models(_bank())[name].eval()
    b2 = _batch(k=2, seed=3)
    pad = lambda t, v: torch.cat([t, torch.full_like(t[:, :1], v)], 1)  # noqa: E731
    F3 = torch.cat([b2.F, torch.randn_like(b2.F[:, :1])], 1)
    b3 = sm.SetBatch(F3, pad(b2.X, 0.0), pad(b2.M, False), b2.invT, pad(b2.mol, 0), b2.prop)
    with torch.no_grad():
        assert torch.allclose(m(b2), m(b3), atol=1e-5)


def test_seed_determinism():
    bank = _bank()
    b = _batch(n=40)
    y = torch.randn(40)
    w = torch.ones(40) / 40
    out = []
    for _ in range(2):
        ens = sm.train_ensemble(lambda: sm.PairDeepSets(width=8), b, y, w, seeds=[0], epochs=3, batch_size=16)
        out.append(sm.predict_ensemble(ens, b)[0])
    assert np.allclose(out[0], out[1])
    assert bank.n_mol == 3


def test_graph_featuriser():
    bank = sm.GraphBank.from_smiles(["CCC", "Cc1ccccc1"])
    assert bank.n_mol == 2
    assert int((bank.atom2mol == 0).sum()) == 3 and int((bank.atom2mol == 1).sum()) == 7
    # propane 2 bonds, toluene 7 bonds; edges are stored in both directions
    src = bank.edge_index[0]
    assert int((bank.atom2mol[src] == 0).sum()) == 4 and int((bank.atom2mol[src] == 1).sum()) == 14


def test_gate_direction():
    assert sm.gate_pass(dict(lo=0.001, hi=0.01)) is True
    assert sm.gate_pass(dict(lo=-0.01, hi=0.01)) is False
    assert sm.gate_pass(dict(lo=-0.02, hi=-0.01)) is False      # significant loss is not a pass


def test_fold_lists_reproduce_groupkfold():
    from sklearn.model_selection import GroupKFold
    g = np.repeat(np.array([f"s{i}" for i in range(23)]), 4)
    rows = np.arange(len(g))
    ref = [sorted(set(g[rows[te]])) for _, te in GroupKFold(n_splits=5).split(rows, groups=g)]
    assert sm.fold_test_systems(g, rows) == ref


def test_system_macro_weights_sum_to_one_per_property():
    s = np.array(["a", "a", "b", "c", "c", "c"])
    p = np.array([0, 0, 0, 1, 1, 1])
    w = sm.system_macro_weights(s, p)
    assert np.isclose(w[p == 0].sum(), 1) and np.isclose(w[p == 1].sum(), 1)
    assert np.isclose(w[0], 0.25) and np.isclose(w[2], 0.5)
