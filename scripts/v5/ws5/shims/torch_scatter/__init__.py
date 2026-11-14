"""Minimal pure-torch stand-in for torch_scatter (no wheel exists for torch 2.15 CPU on Windows).

Only scatter_add and scatter_mean over dim=0 with a 1-D index, which is all GH-GNN's GHGNN_architecture.py uses.
Semantics follow torch_scatter: optional `out` is the initial buffer; scatter_mean divides the accumulated sum by the
per-index count clamped to >= 1.  Written for the WS5 test; not part of the GH-GNN repo.
"""
import torch


def _expand(index, src, dim):
    if dim < 0:
        dim += src.dim()
    shape = [1] * src.dim()
    shape[dim] = -1
    return index.view(shape).expand_as(src), dim


def scatter_add(src, index, dim=0, out=None, dim_size=None):
    idx, dim = _expand(index, src, dim)
    if out is None:
        size = list(src.shape)
        size[dim] = int(dim_size) if dim_size is not None else (int(index.max()) + 1 if index.numel() else 0)
        out = torch.zeros(size, dtype=src.dtype, device=src.device)
    return out.scatter_add_(dim, idx, src)


def scatter_mean(src, index, dim=0, out=None, dim_size=None):
    summed = scatter_add(src, index, dim, out, dim_size)
    d = dim if dim >= 0 else dim + summed.dim()
    ones = torch.ones(index.shape[0], dtype=src.dtype, device=src.device)
    cnt = torch.zeros(summed.shape[d], dtype=src.dtype, device=src.device).scatter_add_(0, index, ones).clamp(min=1)
    shape = [1] * summed.dim()
    shape[d] = -1
    return summed / cnt.view(shape)
