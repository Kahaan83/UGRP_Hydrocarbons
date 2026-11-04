"""Unit tests for Low-Rank Bilinear Pairwise Interactions (T4.2, T4.3)."""
import torch
import numpy as np

from fuelblend.mixture.bilinear_pairwise import BilinearPairwiseInteraction, ClassSharedInteraction


def test_matching_explicit_O_N2_output():
    torch.manual_seed(42)
    batch_size = 3
    N = 6  # Small component count to test exact O(N^2) match
    d = 12 # Descriptor feature dimension
    rank = 4
    num_tasks = 3

    model = BilinearPairwiseInteraction(in_features=d, rank=rank, num_tasks=num_tasks)
    z = torch.randn(batch_size, N, d)
    w = torch.rand(batch_size, N)
    w = w / w.sum(dim=1, keepdim=True)

    out_linear = model.forward_linear_O_N(z, w)
    out_explicit = model.forward_explicit_O_N2(z, w)

    # Must be mathematically identical within float precision
    assert torch.allclose(out_linear, out_explicit, atol=1e-5), "Linear and explicit O(N^2) outputs must match!"


def test_unseen_pairs_get_nonzero_informed_predictions():
    torch.manual_seed(42)
    d = 8
    model = BilinearPairwiseInteraction(in_features=d, rank=4, num_tasks=2)

    # Two distinct synthetic molecular descriptors (e.g. novel cycloalkane and novel aromatic)
    z1 = torch.tensor([[1.2, -0.5, 0.8, 2.1, 0.0, -1.1, 0.4, 0.9]])
    z2 = torch.tensor([[-0.8, 1.4, -0.2, 0.3, 1.5, 0.2, -0.6, -0.4]])

    z_pair = torch.stack([z1, z2], dim=1) # (1, 2, d)
    w_pair = torch.tensor([[0.5, 0.5]])

    out = model(z_pair, w_pair)
    # Output must be nonzero, continuous, and non-trivial
    assert out.shape == (1, 2)
    assert not torch.allclose(out, torch.zeros_like(out), atol=1e-3)
    assert torch.isfinite(out).all()


def test_permutation_invariance():
    torch.manual_seed(42)
    model = BilinearPairwiseInteraction(in_features=6, rank=3, num_tasks=2)

    z = torch.randn(1, 5, 6)
    w = torch.tensor([[0.1, 0.2, 0.3, 0.25, 0.15]])

    out_orig = model(z, w)

    # Permute order of components
    perm = [3, 0, 4, 1, 2]
    z_perm = z[:, perm, :]
    w_perm = w[:, perm]

    out_perm = model(z_perm, w_perm)
    assert torch.allclose(out_orig, out_perm, atol=1e-6), "Model must be permutation invariant!"


def test_class_shared_interaction():
    torch.manual_seed(42)
    model = ClassSharedInteraction(num_classes=5, num_tasks=3)
    class_fracs = torch.tensor([[0.4, 0.3, 0.2, 0.1, 0.0], [0.2, 0.2, 0.2, 0.2, 0.2]])
    out = model(class_fracs)
    assert out.shape == (2, 3)
    assert torch.isfinite(out).all()
