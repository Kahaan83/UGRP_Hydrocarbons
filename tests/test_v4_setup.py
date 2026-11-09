"""v4 Task 0: environment / determinism guards."""
from pathlib import Path

import torch

from fuelblend.mixture.bilinear_pairwise import BilinearPairwiseInteraction

ROOT = Path(__file__).resolve().parents[1]


def _params(m):
    return {k: v.detach().clone() for k, v in m.state_dict().items()}


def test_bilinear_init_is_deterministic():
    a = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3, seed=0)
    b = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3, seed=0)
    pa, pb = _params(a), _params(b)
    assert pa.keys() == pb.keys()
    for k in pa:
        assert torch.equal(pa[k], pb[k]), f"parameter {k} differs between seeded constructions"


def test_bilinear_init_independent_of_global_rng_and_seed_sensitive():
    torch.manual_seed(1)
    a = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3, seed=0)
    torch.manual_seed(999)
    b = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3, seed=0)
    c = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3, seed=1)
    assert torch.equal(a.U, b.U) and torch.equal(a.lambdas, b.lambdas)
    assert not torch.equal(a.U, c.U)


def test_requirements_lower_bounds():
    text = (ROOT / "requirements.txt").read_text()
    for pkg in ("scikit-learn>=1.3", "torch>=2.1", "joblib>=1.3", "pytest>=7.4", "PyYAML>=6.0"):
        assert pkg in text, pkg
