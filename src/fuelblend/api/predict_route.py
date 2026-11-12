"""v5 API route: route context (registry + frozen provider grid + SubReps) and per-property component mapping.

Mirrors the BASE scoring route (``fuelblend.eval.composition_route``): a component missing own data for a property is
replaced, for that property only, by ``SubReps.get``; non-hydrocarbons and components with no representative are unmapped.
The provider is read-only: ``ensure()`` never runs and the npz is never written.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Optional

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi

from fuelblend.baselines.flash_unifac import _lookup
from fuelblend.data.reconstruction_r3 import SubReps, subkind
from fuelblend.eval.composition_route import NEAREST_ANY, RouteCtx, build_tmap
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
GRID_NPZ = ROOT / "data" / "cache" / "pure_grid.npz"
COMPOUNDS_CSV = ROOT / "data" / "unified" / "compounds.csv"
_IK_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


class ReadOnlyProvider(PurePropertyProvider):
    """Frozen pure-property grid: loads the npz, never computes or writes."""

    def __init__(self, reg, cache_path=GRID_NPZ):
        if not Path(cache_path).exists():
            raise FileNotFoundError(f"pure-property grid not found: {cache_path}")
        super().__init__(reg, cache_path, verbose=False)

    def ensure(self, comp_idx):                       # noqa: D401 - deliberate no-op
        return None


@lru_cache(maxsize=1)
def get_ctx() -> RouteCtx:
    reg = CompoundRegistry(COMPOUNDS_CSV)
    prov = ReadOnlyProvider(reg, GRID_NPZ)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    ctx = RouteCtx(reg, prov)
    ctx.reps = reps
    return ctx


@dataclass
class Component:
    ik: str
    smiles: str
    w: float
    reg_idx: Optional[int]
    subkind: Optional[str]
    n_carbon: int
    is_hc: bool


def _mol_info(smi: str):
    m = Chem.MolFromSmiles(smi)
    can = Chem.MolToSmiles(m)
    nc = sum(a.GetSymbol() == "C" for a in m.GetAtoms())
    hc = nc > 0 and all(a.GetSymbol() in ("C", "H") for a in m.GetAtoms())
    return m, can, nc, hc


def resolve(composition, mole_fractions: Optional[bool] = None) -> list[Component]:
    """composition: {id: weight} or [(id, weight), ...]; id is an InChIKey or SMILES.  Weights are mass fractions unless
    ``mole_fractions`` is truthy.  Duplicates are merged by InChIKey, zero weights dropped, weights normalised."""
    ctx = get_ctx()
    reg = ctx.reg
    items = list(composition.items()) if hasattr(composition, "items") else list(composition)
    acc: dict[str, list] = {}
    for key, wt in items:
        wt = float(wt)
        if not np.isfinite(wt) or wt < 0:
            raise ValueError(f"invalid weight {wt!r} for {key!r}")
        key = str(key).strip()
        if _IK_RE.match(key):
            i = reg.ik2i.get(key)
            if i is None:
                raise ValueError(f"InChIKey not in registry: {key}")
            smi = reg.smiles[i]
        else:
            m = Chem.MolFromSmiles(key)
            if m is None:
                raise ValueError(f"unparseable SMILES: {key!r}")
            smi = Chem.MolToSmiles(m)
            key = inchi.MolToInchiKey(m)
            i = reg.ik2i.get(key)
            if i is not None:
                smi = reg.smiles[i]
        if key in acc:
            acc[key][1] += wt
        else:
            acc[key] = [smi, wt, i]
    out = []
    for ik, (smi, wt, i) in acc.items():
        if mole_fractions:
            wt = wt * ctx.mw(smi)
        if wt > 0:
            out.append([ik, smi, wt, i])
    tot = sum(o[2] for o in out)
    if not out or tot <= 0:
        raise ValueError("composition has no positive weight")
    comps = []
    for ik, smi, wt, i in out:
        _, can, nc, hc = _mol_info(smi)
        comps.append(Component(ik=ik, smiles=smi, w=wt / tot, reg_idx=i, subkind=subkind(smi) if hc else None,
                               n_carbon=nc, is_hc=hc))
    return comps


@dataclass
class PropMap:
    idx: np.ndarray
    w: np.ndarray
    own_mass: float
    substituted_mass: float
    unmapped_mass: float
    substitutions: list


def map_for(comps: list, has_own: Callable[[int], bool], reps) -> PropMap:
    """Per-property mapping.  Own data wins; else a SubReps representative (gap <= 3, 99 for NEAREST_ANY kinds); else unmapped.
    Weights are renormalised over the mapped mass; representatives are aggregated as in ``map_bins``."""
    agg: dict[int, float] = {}
    own = sub = lost = 0.0
    subs = []
    for c in comps:
        if c.reg_idx is not None and has_own(c.reg_idx):
            agg[c.reg_idx] = agg.get(c.reg_idx, 0.0) + c.w
            own += c.w
            continue
        rep = None
        if c.is_hc and c.subkind is not None:
            rep = reps.get(c.subkind, c.n_carbon, max_gap=99 if c.subkind in NEAREST_ANY else 3)
        if rep is None:
            lost += c.w
            continue
        agg[rep] = agg.get(rep, 0.0) + c.w
        sub += c.w
        rep_nc = _mol_info(reps.smiles[rep])[2]
        subs.append((c.ik, reps.ik[rep], abs(rep_nc - c.n_carbon)))
    idx = np.array(sorted(agg), dtype=np.int64)
    w = np.array([agg[i] for i in idx], float)
    w = w / w.sum() if w.size and w.sum() > 0 else w
    return PropMap(idx=idx, w=w, own_mass=own, substituted_mass=sub, unmapped_mass=lost, substitutions=subs)


@lru_cache(maxsize=1)
def _tmap() -> dict:
    ctx = get_ctx()
    ok = [i for i in range(len(ctx.reg.ik)) if ctx.reg.usable[i]]
    return build_tmap(ctx.reg, ok)


@lru_cache(maxsize=1)
def _ysi_smiles() -> frozenset:
    from fuelblend.baselines.sooting import ysi_table      # data/processed/v4/sooting/ysi_hc.csv
    return frozenset(ysi_table())


def predicates(ctx: Optional[RouteCtx] = None) -> dict:
    """Own-data predicates over registry indices: HAS_GRID, HAS_FLASH, HAS_SLE, HAS_YSI."""
    ctx = ctx or get_ctx()
    reg, prov = ctx.reg, ctx.prov

    def has_sle(i):
        t = _tmap().get(reg.ik[i])
        return bool(t) and bool(np.isfinite(t.get("tm_K", np.nan)) and np.isfinite(t.get("dh_fus_J_mol", np.nan)))

    def has_ysi(i):
        m = Chem.MolFromSmiles(reg.smiles[i])
        return m is not None and Chem.MolToSmiles(m) in _ysi_smiles()

    return {"HAS_GRID": lambda i: bool(prov.valid[i]),
            "HAS_FLASH": lambda i: _lookup(reg.smiles[i], None) is not None,
            "HAS_SLE": has_sle, "HAS_YSI": has_ysi}
