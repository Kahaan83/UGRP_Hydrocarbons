"""Low-Rank Bilinear Feature-Based Pairwise Interaction Model (T4.2, T4.3).

Computes component-component interactions in O(N * d * r) rather than O(N^2),
enabling scale-up to 40+ components with continuous descriptor-based generalization.

Mathematical formulation:
  Explicit O(N^2) formulation:
      E = sum_i sum_j (w_i * w_j * z_i^T * W * z_j)
  With low-rank decomposition W = U * diag(lambda) * U^T (rank r):
      E = sum_k lambda_k * [ sum_i w_i * (u_k^T * z_i) ]^2
  which evaluates in strictly O(N * d * r) time.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import numpy as np


class BilinearPairwiseInteraction(nn.Module):
    def __init__(self, in_features: int, rank: int = 16, num_tasks: int = 3, seed: int = 0):
        super().__init__()
        self.in_features = in_features
        self.rank = rank
        self.num_tasks = num_tasks
        self.seed = seed

        # Local generator: initialization is deterministic given `seed` and independent of
        # (and does not consume) the global torch RNG (v4 Task 0 / CLM-19 fix).
        gen = torch.Generator().manual_seed(int(seed))

        # Projection matrix U: (in_features, rank)
        self.U = nn.Parameter(torch.randn(in_features, rank, generator=gen) * (1.0 / np.sqrt(in_features)))
        # Task-specific eigenvalues lambda: (num_tasks, rank)
        self.lambdas = nn.Parameter(torch.randn(num_tasks, rank, generator=gen) * 0.05)
        # Linear self-term (single-component contribution). nn.Linear's default init draws from
        # the global RNG, so re-draw its weight from `gen` with the same U(-b, b), b = 1/sqrt(fan_in).
        self.linear_self = nn.Linear(in_features, num_tasks, bias=False)
        bound = 1.0 / np.sqrt(in_features)
        with torch.no_grad():
            self.linear_self.weight.copy_(
                (torch.rand(num_tasks, in_features, generator=gen) * 2.0 - 1.0) * bound
            )

    def forward_linear_O_N(self, z: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        """Compute interaction in O(N * d * r) linear time.

        Parameters:
          z: Component feature tensor of shape (batch, N, in_features)
          w: Component fraction tensor of shape (batch, N)

        Returns:
          Interaction energy tensor of shape (batch, num_tasks)
        """
        # 1. Project features: P = z @ U -> shape (batch, N, rank)
        P = torch.matmul(z, self.U)

        # 2. Pool weighted projections: V_k = sum_i w_i * P_{i, k} -> shape (batch, rank)
        w_expanded = w.unsqueeze(-1)  # (batch, N, 1)
        V = torch.sum(w_expanded * P, dim=1)  # (batch, rank)

        # 3. Bilinear interaction: E_task = sum_k lambda_{task, k} * V_k^2
        # V^2 shape: (batch, rank), lambdas shape: (num_tasks, rank)
        V_sq = V ** 2  # (batch, rank)
        E_pair = torch.matmul(V_sq, self.lambdas.t())  # (batch, num_tasks)

        # 4. Self linear contribution: sum_i w_i * Linear(z_i)
        self_term = torch.sum(w_expanded * self.linear_self(z), dim=1)  # (batch, num_tasks)

        return E_pair + self_term

    def forward_explicit_O_N2(self, z: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        """Explicit O(N^2) double-loop calculation for mathematical verification."""
        batch_size, N, d = z.shape
        W_tasks = torch.matmul(self.U, torch.matmul(torch.diag_embed(self.lambdas), self.U.t()))  # (num_tasks, d, d)

        out = torch.zeros(batch_size, self.num_tasks, device=z.device, dtype=z.dtype)
        for b in range(batch_size):
            for task in range(self.num_tasks):
                W_t = W_tasks[task]  # (d, d)
                task_sum = 0.0
                for i in range(N):
                    for j in range(N):
                        term = w[b, i] * w[b, j] * torch.matmul(z[b, i], torch.matmul(W_t, z[b, j]))
                        task_sum = task_sum + term
                out[b, task] = task_sum

        w_expanded = w.unsqueeze(-1)
        self_term = torch.sum(w_expanded * self.linear_self(z), dim=1)
        return out + self_term

    def forward(self, z: torch.Tensor, w: torch.Tensor) -> torch.Tensor:
        return self.forward_linear_O_N(z, w)


class ClassSharedInteraction(nn.Module):
    """Class-level shared interaction parameters (T4.3 alternative parameterization).

    Classes: n_paraffin, iso_paraffin, cycloalkane, aromatic, alkene_alkyne (5 classes).
    Parameters: Symmetric 5x5 interaction matrix per task.
    """
    def __init__(self, num_classes: int = 5, num_tasks: int = 3):
        super().__init__()
        self.num_classes = num_classes
        self.num_tasks = num_tasks
        # Symmetric interaction matrices
        self.M_raw = nn.Parameter(torch.randn(num_tasks, num_classes, num_classes) * 0.05)

    def forward(self, class_fractions: torch.Tensor) -> torch.Tensor:
        """
        Parameters:
          class_fractions: Mass fractions grouped by chemical class, shape (batch, num_classes)
        Returns:
          Interaction energy, shape (batch, num_tasks)
        """
        # Symmetrize M
        M = 0.5 * (self.M_raw + self.M_raw.transpose(1, 2))  # (tasks, C, C)
        batch_size = class_fractions.shape[0]
        out = torch.zeros(batch_size, self.num_tasks, device=class_fractions.device)
        for t in range(self.num_tasks):
            # c @ M_t @ c^T
            quad = torch.matmul(class_fractions, torch.matmul(M[t], class_fractions.t()))
            out[:, t] = torch.diagonal(quad)
        return out
