"""Isomer-averaged tile constants (Hall et al. 2022 idea: a GC x GC tile is a probability-weighted mean over its isomers, not one
stand-in compound).  Each tile representative (``SubReps.by[(kind, n)][0]``, the lowest-InChIKey stand-in) gets the mean Tm and dHfus of
every registry compound of the same (sub-kind, carbon number) that has finite constants in the given tmap (equal weights).
SubReps and v5 behaviour are untouched; this builds a modified tmap copy."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]


def get_registry_and_reps():
    """Registry plus the same R3 representatives (SubReps) the v5 G1 mapping uses."""
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    return reg, SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))


def tile_pool(reg, tmap: dict) -> dict:
    """(sub-kind, carbon number) -> list of (tm_K, dh_fus_J_mol) over all registry compounds with finite constants."""
    from fuelblend.data.reconstruction_r3 import subkind
    pool = defaultdict(list)
    for ik, s in zip(reg.ik, reg.smiles):
        t = tmap.get(ik)
        if not t:
            continue
        tm, dh = t.get("tm_K", np.nan), t.get("dh_fus_J_mol", np.nan)
        if not (np.isfinite(tm) and np.isfinite(dh)):
            continue
        m = Chem.MolFromSmiles(s)
        k = subkind(s) if m is not None else None
        if k is None:
            continue
        pool[(k, sum(a.GetSymbol() == "C" for a in m.GetAtoms()))].append((float(tm), float(dh)))
    return pool


def tile_mean_tmap(reg, reps, tmap: dict, min_isomers: int = 2) -> dict:
    """tmap copy where each tile representative carries the isomer-mean constants.  Tiles with fewer than ``min_isomers`` are unchanged."""
    out = dict(tmap)
    pool = tile_pool(reg, tmap)
    for key, members in reps.by.items():
        vals = pool.get(key, [])
        if len(vals) < min_isomers:
            continue
        rep = reps.ik[members[0]]
        base = dict(tmap.get(rep) or {})
        base.update(tm_K=float(np.mean([v[0] for v in vals])), dh_fus_J_mol=float(np.mean([v[1] for v in vals])), n_isomers=len(vals))
        out[rep] = base
    return out
