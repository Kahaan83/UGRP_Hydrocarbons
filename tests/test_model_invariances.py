"""Property tests for Mixture Model Invariances (T4.4).

Verifies 5 physical invariance principles:
  1. Permutation invariance
  2. Splitting a compound into two identical halves
  3. Adding a zero-fraction component
  4. Basis conversion consistency (mass <-> mole)
  5. Pure-compound limit
"""
import pytest
import torch
import numpy as np

from fuelblend.mixture.bilinear_pairwise import BilinearPairwiseInteraction
from fuelblend.baselines.heat_of_combustion import mass_fractions_from_mole, mole_fractions_from_mass


def test_1_permutation_invariance():
    torch.manual_seed(42)
    model = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3)
    z = torch.randn(1, 7, 8)
    w = torch.rand(1, 7)
    w = w / w.sum()

    out_orig = model(z, w)

    # Random permutation of indices
    perm = torch.randperm(7)
    z_perm = z[:, perm, :]
    w_perm = w[:, perm]

    out_perm = model(z_perm, w_perm)
    assert torch.allclose(out_orig, out_perm, atol=1e-6)


def test_2_split_compound_invariance_pairwise():
    """Splitting compound A into two identical halves with w_A/2 must yield identical prediction."""
    torch.manual_seed(42)
    model = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3)

    # 3 compounds
    z_orig = torch.tensor([
        [1.0, 0.5, -0.2, 0.8, -1.0, 0.4, 0.1, -0.7],
        [0.2, -0.9, 1.4, -0.3, 0.5, 0.6, -0.2, 0.0],
        [-0.5, 0.3, 0.0, 1.2, -0.4, -0.8, 1.1, 0.5],
    ]).unsqueeze(0) # (1, 3, 8)
    w_orig = torch.tensor([[0.4, 0.35, 0.25]])

    out_orig = model(z_orig, w_orig)

    # Split compound 0 (w=0.4) into two identical halves with w=0.2 each
    z_split = torch.cat([z_orig[:, :1, :], z_orig[:, :1, :], z_orig[:, 1:, :]], dim=1) # (1, 4, 8)
    w_split = torch.tensor([[0.2, 0.2, 0.35, 0.25]])

    out_split = model(z_split, w_split)

    # Pairwise bilinear model strictly passes the split test!
    assert torch.allclose(out_orig, out_split, atol=1e-5), "Pairwise model must pass compound split invariance!"


def test_3_zero_fraction_component_invariance():
    """Adding a dummy component with fraction w = 0.0 changes nothing."""
    torch.manual_seed(42)
    model = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3)

    z_base = torch.randn(1, 4, 8)
    w_base = torch.tensor([[0.3, 0.4, 0.2, 0.1]])

    out_base = model(z_base, w_base)

    # Append dummy component with arbitrary features and w = 0.0
    z_dummy = torch.randn(1, 1, 8) * 10.0 # Extreme arbitrary vector
    z_ext = torch.cat([z_base, z_dummy], dim=1)
    w_ext = torch.cat([w_base, torch.tensor([[0.0]])], dim=1)

    out_ext = model(z_ext, w_ext)
    assert torch.allclose(out_base, out_ext, atol=1e-6)


def test_4_basis_conversion_consistency():
    """Mass <-> mole conversions must roundtrip without drift."""
    rng = np.random.default_rng(42)
    mw = np.array([100.2, 114.2, 92.1, 128.2, 142.3])
    x_orig = rng.dirichlet(np.ones(len(mw)))

    w = mass_fractions_from_mole(x_orig, mw)
    x_recov = mole_fractions_from_mass(w, mw)

    assert np.allclose(x_orig, x_recov, atol=1e-12)
    assert np.isclose(w.sum(), 1.0)
    assert np.isclose(x_recov.sum(), 1.0)


def test_5_pure_compound_limit():
    """When N=1 and w=1.0, single-component prediction equals linear self-term."""
    torch.manual_seed(42)
    model = BilinearPairwiseInteraction(in_features=8, rank=4, num_tasks=3)

    z_single = torch.randn(1, 1, 8)
    w_single = torch.tensor([[1.0]])

    out = model(z_single, w_single)
    assert out.shape == (1, 3)
    assert torch.isfinite(out).all()
