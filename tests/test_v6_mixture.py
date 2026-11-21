"""Track C (v6) tests: MolPool maths, lock integrity, baseline sanity."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.ml import molpool as MP  # noqa: E402

OUT = ROOT / "results/v6/mixture"


def test_untrained_head_is_baseline():
    net = MP.MolPool(10, d=4, rank=3)
    net.eval()
    H = torch.randn(5, 3, 10)
    X = torch.softmax(torch.randn(5, 3), -1)
    assert torch.allclose(net(H, X), torch.zeros(5))


def test_pair_term_matches_bruteforce_and_padding_is_inert():
    torch.manual_seed(0)
    net = MP.MolPool(6, d=4, rank=3)
    with torch.no_grad():
        net.out[-1].weight.copy_(torch.eye(1, net.out[-1].weight.shape[1]).fill_(0))
        net.out[-1].weight[0, 4] = 1.0       # read the pair feature q only (index d = 4)
        net.out[0].weight.copy_(torch.eye(32, 5))   # pass z0..z3 and q straight through the hidden layer
        net.out[0].bias.zero_()
    net.eval()
    H = torch.randn(2, 3, 6)
    X = torch.tensor([[0.5, 0.3, 0.2], [0.6, 0.4, 0.0]])
    with torch.no_grad():
        phi = net.phi(H)
        p = net.P(phi)
        want = []
        for b in range(2):
            q = 0.0
            for i in range(3):
                for j in range(i + 1, 3):
                    q += float(X[b, i] * X[b, j] * (p[b, i] * p[b, j]).sum())
            want.append(q)
        got = net(H, X).numpy()
    assert np.allclose(got, want, atol=1e-5)
    # a padded slot (x = 0) must not change the output
    with torch.no_grad():
        H2 = H.clone()
        H2[1, 2] = torch.randn(6)
        assert np.allclose(net(H2, X).numpy(), got, atol=1e-6)


def test_permutation_invariance():
    net = MP.MolPool(6, d=4, rank=3)
    with torch.no_grad():
        net.out[-1].weight.normal_()
    net.eval()
    H = torch.randn(1, 4, 6)
    X = torch.softmax(torch.randn(1, 4), -1)
    perm = torch.tensor([2, 0, 3, 1])
    with torch.no_grad():
        assert torch.allclose(net(H, X), net(H[:, perm], X[:, perm]), atol=1e-6)


def test_pad_sets():
    H, X = MP.pad_sets([np.ones((2, 3)), np.ones((1, 3))], [np.array([.5, .5]), np.array([1.0])])
    assert H.shape == (2, 2, 3) and X[1, 1] == 0 and H[1, 1].sum() == 0


@pytest.mark.skipif(not (OUT / "lock.json").exists(), reason="lock not built")
def test_lock_integrity():
    lock = json.load(open(OUT / "lock.json"))
    assert hashlib.sha256((OUT / "rows.csv").read_bytes()).hexdigest() == lock["rows_csv_sha256"]
    df = pd.read_csv(OUT / "rows.csv")
    for head, g in df.groupby("head"):
        locked = set(lock["heads"][head]["locked"])
        assert set(g[g.role == "locked"].group.astype(str)) == {str(s) for s in locked}
        assert not set(g[g.role == "dev"].group) & set(g[g.role == "locked"].group)
        assert (g.role == "locked").sum() <= lock["row_max"] * len(g) + 1e-9
        assert (g.role == "locked").sum() > 0
    # each (head, prop, group, composition, T) appears once (dedup)
    assert not df.duplicated(["head", "prop", "group", "comps", "x", "T_K"]).any()


@pytest.mark.skipif(not (OUT / "rows.csv").exists(), reason="rows not built")
def test_baseline_sanity():
    df = pd.read_csv(OUT / "rows.csv")
    for s in df.x:
        assert abs(sum(json.loads(s)) - 1.0) < 1e-5
    o = df[df["head"] != "boehm"]
    assert o.base.between(-10, 130).all() and o["cov"].ge(0.9).all()
    b = df[df["head"] == "boehm"]
    assert (b[b.prop == "density"].base.between(500, 1100)).all()
    assert np.isfinite(b.base).all() and b.T_K.notna().all()
    # Refutas of a mixture lies between its pure ln-viscosities is checked through the y/base scale: both finite and < 6 ln units apart
    v = b[b.prop == "ln_visc"]
    assert (v.y - v.base).abs().max() < 1.5
