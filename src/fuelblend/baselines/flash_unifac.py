"""UNIFAC-corrected Le Chatelier flash point (v4 Task 9).

Model (same solver, units and conventions as ``baselines/flash.py``):
    H(T) = sum_i x_i * gamma_i(T, x) * Psat_i(T) / (P * LFL_i(T)) = 1       (Psat kPa, T K, flash degC)

* Psat: the repo's anchored Antoine constants (``data/processed/v2/compounds_anchored_vp.csv``,
  log10(P/kPa) = A - B/(T/K + C), anchored on pure-liquid vapour pressure at 273-373 K; no flash data used).
* LFL: the same table (measured where available, else the Jones 0.55*Cst estimate), Zabetakis T-correction.
* gamma: ORIGINAL (Hansen 1991) UNIFAC from ``thermo.unifac.UNIFAC`` with ``UFSG``/``UFIP``, version=0.
  NOT UNIFAC-Dortmund.  Groups come from the SMARTS-free RDKit fragmenter ``unifac_subgroups`` below, which
  covers hydrocarbons only: CH3/CH2/CH/C (1-4), alkenes 5-8 and 70, ACH/AC (9,10), ACCH3/ACCH2/ACCH (11-13).

Documented fallback: a molecule the fragmenter cannot assign (alkyne, allene, heteroatom, carbon bonded to two
aromatic carbons, ...) gets gamma = 1 *for that molecule*; UNIFAC is evaluated on the remaining (assignable)
molecules renormalised to sum 1.  This is an approximation, never silent: a ``UnifacFallbackWarning`` is emitted
and ``flash_point_unifac_detail`` reports ``fallback_components``.
"""
from __future__ import annotations

import warnings
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem

from fuelblend.baselines.flash import P_TOTAL_KPA, predict_flash_le_chatelier

REGISTRY_CSV = Path(__file__).resolve().parents[3] / "data" / "processed" / "v2" / "compounds_anchored_vp.csv"
UNIFAC_VERSION = "original UNIFAC (Hansen 1991; thermo.unifac UFSG/UFIP, version=0)"


class UnifacFallbackWarning(UserWarning):
    """A molecule has no UNIFAC group assignment; gamma = 1 was used for it."""


# --------------------------------------------------------------------------- group fragmentation
_ALKENE = {(2, 1): 5, (1, 1): 6, (2, 0): 7, (1, 0): 8, (0, 0): 70}   # (H on a, H on b) sorted high->low -> subgroup


@lru_cache(maxsize=4096)
def _subgroups_cached(canon_smiles: str) -> Optional[Tuple[Tuple[int, int], ...]]:
    mol = Chem.MolFromSmiles(canon_smiles)
    if mol is None:
        return None
    # isotope-labelled hydrogens ([2H]) are explicit atoms in the graph: treat them as ordinary H for grouping
    for a in mol.GetAtoms():
        a.SetIsotope(0)
    mol = Chem.RemoveHs(mol)
    if any(a.GetSymbol() not in ("C", "H") or a.GetFormalCharge() != 0 for a in mol.GetAtoms()):
        return None
    counts: Dict[int, int] = {}

    def add(k: int, n: int = 1) -> None:
        counts[k] = counts.get(k, 0) + n

    consumed = set()                       # aromatic carbons absorbed into an ACCHn group
    done_alkene = set()
    for a in mol.GetAtoms():
        if a.GetIsAromatic():
            continue
        idx, h = a.GetIdx(), a.GetTotalNumHs()
        bonds = list(a.GetBonds())
        orders = [b.GetBondTypeAsDouble() for b in bonds]
        if any(o == 3.0 for o in orders):
            return None                    # alkyne
        n_double = sum(o == 2.0 for o in orders)
        if n_double > 1:
            return None                    # allene / cumulene
        if n_double == 1:
            if idx in done_alkene:
                continue
            b2 = next(b for b in bonds if b.GetBondTypeAsDouble() == 2.0)
            o = b2.GetOtherAtom(a)
            if o.GetIsAromatic() or sum(bb.GetBondTypeAsDouble() == 2.0 for bb in o.GetBonds()) != 1:
                return None
            key = tuple(sorted((h, o.GetTotalNumHs()), reverse=True))
            if key not in _ALKENE:
                return None
            add(_ALKENE[key])
            done_alkene.update((idx, o.GetIdx()))
            continue
        # sp3 carbon
        arom_nb = [n for n in a.GetNeighbors() if n.GetIsAromatic()]
        if len(arom_nb) > 1:
            return None
        if len(arom_nb) == 1:
            if h not in (1, 2, 3):
                return None
            add({3: 11, 2: 12, 1: 13}[h])
            consumed.add(arom_nb[0].GetIdx())
        else:
            add({3: 1, 2: 2, 1: 3, 0: 4}[h])
    for a in mol.GetAtoms():
        if not a.GetIsAromatic() or a.GetIdx() in consumed:
            continue
        add(9 if a.GetTotalNumHs() == 1 else 10)
    if not counts:
        return None
    return tuple(sorted(counts.items()))


def unifac_subgroups(smiles: str) -> Optional[Dict[int, int]]:
    """Original-UNIFAC subgroup counts {thermo UFSG id: count} for a hydrocarbon SMILES, or None if unassignable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    res = _subgroups_cached(Chem.MolToSmiles(mol))
    return None if res is None else dict(res)


# --------------------------------------------------------------------------- registry
@lru_cache(maxsize=1)
def _default_registry() -> pd.DataFrame:
    reg = pd.read_csv(REGISTRY_CSV)
    reg["can"] = reg["smiles"].map(lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)))
    return reg.drop_duplicates("can").set_index("can")


def _lookup(smiles: str, registry: Optional[pd.DataFrame]):
    reg = _default_registry() if registry is None else registry
    can = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    if can not in reg.index:
        return None
    r = reg.loc[can]
    if pd.isna(r["antoine_A"]) or pd.isna(r["lfl_vol_pct"]):
        return None
    return (float(r["antoine_A"]), float(r["antoine_B"]), float(r["antoine_C"])), float(r["lfl_vol_pct"])


# --------------------------------------------------------------------------- gamma
def make_unifac_gamma_fn(smiles: Sequence[str]) -> Tuple[Callable, List[int]]:
    """Return (gamma_fn(T_K, x) -> np.ndarray, fallback_component_indices)."""
    from thermo.unifac import UFIP, UFSG, UNIFAC

    groups = [unifac_subgroups(s) for s in smiles]
    ok = [i for i, g in enumerate(groups) if g is not None]
    fallback = [i for i, g in enumerate(groups) if g is None]
    n = len(smiles)

    def gamma_fn(T_K: float, x: Sequence[float]) -> np.ndarray:
        gam = np.ones(n)
        if len(ok) == 0:
            return gam
        # components with x = 0 are dropped from the UNIFAC call (thermo divides by x) and keep gamma = 1; they add
        # nothing to the hazard sum.  Their absence does not change the other gammas (x_i = 0 is a null species).
        act = [i for i in ok if x[i] > 0.0]
        if len(act) < 2:
            return gam                      # single assignable species: gamma = 1 (pure / infinite-dilution neighbours unknown)
        xs = np.array([x[i] for i in act], float)
        xs = xs / xs.sum()
        g = UNIFAC.from_subgroups(T=float(T_K), xs=list(xs), chemgroups=[groups[i] for i in act],
                                  subgroups=UFSG, interaction_data=UFIP, version=0).gammas()
        for i, gi in zip(act, g):
            gam[i] = gi
        return gam

    return gamma_fn, fallback


# --------------------------------------------------------------------------- public API
def _prepare(smiles, x, registry):
    x = np.asarray(x, float)
    if len(smiles) != len(x):
        raise ValueError("smiles and x must have equal length")
    if x.sum() <= 0:
        raise ValueError("x must have a positive sum")
    x = x / x.sum()
    ant, lfl, missing = [], [], []
    for i, s in enumerate(smiles):
        r = _lookup(s, registry)
        if r is None:
            missing.append(i)
            ant.append((np.nan,) * 3)
            lfl.append(np.nan)
        else:
            ant.append(r[0])
            lfl.append(r[1])
    return x, ant, lfl, missing


def flash_point_unifac_detail(smiles: Sequence[str], x: Sequence[float], P_kPa: float = P_TOTAL_KPA,
                              registry: Optional[pd.DataFrame] = None, use_unifac: bool = True) -> dict:
    """Flash point with diagnostics: flash_C, status, gamma_at_flash, fallback_components, missing_inputs."""
    x, ant, lfl, missing = _prepare(smiles, x, registry)
    out = dict(flash_C=np.nan, status="ok", gamma_at_flash=np.full(len(x), np.nan), fallback_components=[],
               missing_inputs=missing, model="unifac" if use_unifac else "ideal")
    if missing:
        out["status"] = "missing_antoine_or_lfl"
        return out
    gamma_fn, fallback = make_unifac_gamma_fn(smiles) if use_unifac else (None, [])
    out["fallback_components"] = fallback
    if use_unifac and fallback:
        warnings.warn(f"no UNIFAC groups for component(s) {[smiles[i] for i in fallback]}: gamma=1 used for them",
                      UnifacFallbackWarning, stacklevel=3)
    # P_kPa: total pressure in the LFL/hazard ratio. predict_flash_le_chatelier uses P_TOTAL_KPA, so rescale
    # LFL-equivalent via Psat/P: H is linear in 1/P, equivalent to scaling LFL by P/P_TOTAL.
    lfl_eff = [v * (P_kPa / P_TOTAL_KPA) for v in lfl]
    T, status = predict_flash_le_chatelier(list(x), ant, lfl_eff, lfl_temp=True, gamma_fn=gamma_fn)
    out["status"] = status
    if status == "ok":
        out["flash_C"] = T
        out["gamma_at_flash"] = (np.ones(len(x)) if gamma_fn is None else gamma_fn(T + 273.15, x))
    return out


def flash_point_unifac(smiles: Sequence[str], x: np.ndarray, P_kPa: float = P_TOTAL_KPA) -> float:
    """Flash point (degC) of a liquid hydrocarbon blend: Le Chatelier + modified Raoult with original-UNIFAC gamma.

    ``x`` are LIQUID MOLE fractions (normalised internally).  Returns NaN if a component lacks Antoine/LFL
    data or the solver fails (see ``flash_point_unifac_detail`` for the reason)."""
    return float(flash_point_unifac_detail(smiles, x, P_kPa)["flash_C"])


def flash_point_ideal(smiles: Sequence[str], x: np.ndarray, P_kPa: float = P_TOTAL_KPA) -> float:
    """Same model with gamma = 1 (ideal Raoult + Le Chatelier), same inputs and solver."""
    return float(flash_point_unifac_detail(smiles, x, P_kPa, use_unifac=False)["flash_C"])
