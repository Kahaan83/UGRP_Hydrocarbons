import torch
from fuelblend.ml import dmpnn as D
from fuelblend.ml import param_heads as P


def _batch(smiles):
    return D.collate([D.featurize(s) for s in smiles])


def test_antoine_head_shape():
    net = P.ParamHeadDMPNN(hidden=32, depth=2, ffn=32).to(D.DEVICE).eval()
    par = net.params(_batch(["CCCCCC", "c1ccccc1", "CC(C)C"]))
    assert par.shape == (3, 5)
    assert (par[:, 1] > 0).all() and (par[:, 2] < 0).all() and (par[:, 2] > -100).all()
    T = torch.full((3,), 150.0, device=par.device)
    assert (T + par[:, 2] > 50).all()


def test_andrade_monotone_in_T():
    net = P.ParamHeadDMPNN(hidden=32, depth=2, ffn=32).to(D.DEVICE).eval()
    par = net.params(_batch(["CCCCCC", "c1ccccc1"]))
    T = torch.tensor([250.0, 300.0, 350.0], device=par.device)
    mu = torch.stack([P.ln_mu(par, torch.full((2,), float(t), device=par.device)) for t in T])
    assert (mu[0] > mu[1]).all() and (mu[1] > mu[2]).all()                    # b > 0: viscosity falls with T
