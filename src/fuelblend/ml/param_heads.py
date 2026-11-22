"""v7 Task 18: D-MPNN with equation-parameter heads. Outputs Antoine (A, B, C) and Andrade (a, b) per molecule.

log10(P/kPa) = A - B/(T_K + C);   ln(mu/mPa.s) = a + b/T_K.  B > 0 and b > 0 by softplus; C = -100 sigmoid(.) so T_K + C > 50 K for T_K >= 150.
"""
from __future__ import annotations

import torch
import torch.nn.functional as F

from fuelblend.ml.dmpnn import DMPNN

_SP1 = 0.5413248546129181                       # softplus(_SP1) = 1


class ParamHeadDMPNN(DMPNN):
    def __init__(self, hidden: int = 200, depth: int = 4, dropout: float = 0.1, ffn: int = 200):
        super().__init__(5, hidden, depth, dropout, ffn)

    def params(self, batch) -> torch.Tensor:
        """(N, 5) columns A, B, C, a, b."""
        r = self.forward(batch)
        A = 6.0 + 3.0 * r[:, 0]
        B = 1500.0 * F.softplus(r[:, 1] + _SP1)
        C = -100.0 * torch.sigmoid(r[:, 2])
        a = -4.0 + 2.0 * r[:, 3]
        b = 1000.0 * F.softplus(r[:, 4] + _SP1)
        return torch.stack([A, B, C, a, b], dim=1)


def log10_p(par: torch.Tensor, T_K: torch.Tensor) -> torch.Tensor:
    return par[:, 0] - par[:, 1] / (T_K + par[:, 2])


def ln_mu(par: torch.Tensor, T_K: torch.Tensor) -> torch.Tensor:
    return par[:, 3] + par[:, 4] / T_K
