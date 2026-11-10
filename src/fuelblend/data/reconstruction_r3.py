"""r3 reconstruction of CRC 647 compositions from DESCRIPTORS ONLY (Table G.1 hydrocarbon types + D86 + olefins).

r2 split the saturates with one fixed n / iso / cycloalkane prior for every fuel and had no two-ring saturate or
indan / tetralin class.  r3 takes each fuel's own measured split from CRC 647 Table G.1: ASTM D6379 totals
(aromatics, saturates, vol%) and ASTM D2425 hydrocarbon-type shares.  D2425 does not separate normal from branched
paraffins, so that one split is still a prior.  Classes with no realistic fuel molecule go to named unexplained
buckets.  Nothing in this module reads a target property.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
from rdkit import Chem

from fuelblend.data.reconstruction_r2 import is_realistic, nominal_sat_split_detail
from fuelblend.mixture.fuel_hindcast import RHO_CLASS, carbon_distribution

log = logging.getLogger(__name__)

SUBKINDS = ("n_paraffin", "iso_paraffin", "cyclo_mono", "cyclo_di", "alkylbenzene", "indan_tetralin",
            "aromatic_di", "alkene_alkyne")
SAT_TYPES = ("paraffins", "cycloparaffins", "dicycloparaffins", "tricycloparaffins")
ARO_TYPES = ("alkylbenzenes", "indans_tetralins", "indenes", "naphthalene", "naphthalenes", "acenaphthenes",
             "acenaphthylenes", "tricyclic_aromatics")
UNEXPLAINED = ("unexplained_tricycloparaffins", "unexplained_minor_aromatics")
SPREAD_DEFAULT = ("n_paraffin", "iso_paraffin", "cyclo_mono", "alkylbenzene", "alkene_alkyne")
SPREAD_NEAREST = ("cyclo_di", "indan_tetralin", "aromatic_di")      # few registry molecules: nearest carbon number


def subkind(smiles: str) -> str | None:
    """Fine class of a realistic fuel molecule (None when it is not one of SUBKINDS)."""
    m = Chem.MolFromSmiles(smiles)              # keep the Mol alive: RingInfo of a temporary reports 0 rings
    if m is None:
        return None
    for k in ("n_paraffin", "iso_paraffin", "alkene_alkyne", "aromatic_di"):
        if is_realistic(smiles, k):
            return k
    ri = m.GetRingInfo()
    n_rings = ri.NumRings()
    if is_realistic(smiles, "cycloalkane"):
        return {1: "cyclo_mono", 2: "cyclo_di"}.get(n_rings)
    if is_realistic(smiles, "aromatic_mono"):
        if n_rings == 1:
            ok = all(b.GetIsAromatic() or b.GetBondType() == Chem.BondType.SINGLE for b in m.GetBonds())
            return "alkylbenzene" if ok else None
        if n_rings == 2:
            r0, r1 = (set(r) for r in ri.BondRings())
            return "indan_tetralin" if r0 & r1 else None
    return None


class SubReps:
    """One realistic registry compound per (subkind, carbon number) bin (lowest InChIKey first).

    Same candidate requirements as r2: usable in the registry, valid pure eta(T)/rho(T) in the provider, and
    ``vp_lookup(smiles)`` is not None.
    """

    def __init__(self, reg, provider, vp_lookup, nc=None):
        self.ik, self.smiles = list(reg.ik), list(reg.smiles)
        has_extra = np.array([vp_lookup(s) is not None for s in self.smiles])
        if nc is None:
            nc = []
            for s in self.smiles:
                m = Chem.MolFromSmiles(s)
                nc.append(sum(a.GetSymbol() == "C" for a in m.GetAtoms()) if m is not None else np.nan)
            nc = np.array(nc, dtype=float)
        cand = [i for i in range(len(self.ik)) if reg.usable[i] and has_extra[i] and np.isfinite(nc[i])]
        provider.ensure(cand)
        self.by: dict[tuple[str, int], list[int]] = {}
        for i in (i for i in cand if provider.valid[i]):
            k = subkind(self.smiles[i])
            if k is not None:
                self.by.setdefault((k, int(nc[i])), []).append(i)
        for k in self.by:
            self.by[k].sort(key=lambda i: self.ik[i])
        self.ns = {k: sorted(n for (kk, n) in self.by if kk == k) for k in SUBKINDS}

    def get(self, kind: str, n_carbon: int, max_gap: int = 3):
        ns = self.ns.get(kind, [])
        if not ns:
            return None
        best = min(ns, key=lambda k: (abs(k - n_carbon), k))
        return self.by[(kind, best)][0] if abs(best - n_carbon) <= max_gap else None


def n_share_prior(jet_class_csv: Path) -> float:
    """n / (n + iso) of the GC jet fuels (a prior: Table G.1 does not separate normal from branched paraffins)."""
    (n, iso, _), _, _ = nominal_sat_split_detail(jet_class_csv)
    return float(n / (n + iso))


def class_masses(g1: Mapping, olefins_vol_pct: float, n_share: float) -> dict | None:
    """Mass fraction of the fuel per class (SUBKINDS plus UNEXPLAINED), summing to 1; None if the inputs are unusable.

    ``g1`` is one fuel's {property: value} from Table G.1 with below-detection cells already set to 0.
    """
    A = float(g1.get("aromatics_vol_pct", np.nan))
    if not np.isfinite(A):
        log.warning("class_masses: no aromatics value")
        return None
    O = float(olefins_vol_pct) if olefins_vol_pct is not None and np.isfinite(olefins_vol_pct) else 0.0
    S = max(100.0 - A - O, 0.0)
    t = {}
    for k in SAT_TYPES + ARO_TYPES:
        v = float(g1.get(f"hc_type_{k}_pct", np.nan))
        if not np.isfinite(v):
            log.warning("class_masses: missing D2425 cell %s", k)
            return None
        t[k] = v
    sat_sum, aro_sum = sum(t[k] for k in SAT_TYPES), sum(t[k] for k in ARO_TYPES)
    if sat_sum <= 0 or (A > 0 and aro_sum <= 0):
        log.warning("class_masses: D2425 type sum is not positive (saturates %.3g, aromatics %.3g)", sat_sum, aro_sum)
        return None
    m_aro, m_ole, m_sat = A * RHO_CLASS["aromatic"], O * RHO_CLASS["olefin"], S * RHO_CLASS["sat"]
    tot = m_aro + m_ole + m_sat
    if tot <= 0:
        return None
    m_aro, m_ole, m_sat = m_aro / tot, m_ole / tot, m_sat / tot
    sa = {k: m_sat * t[k] / sat_sum for k in SAT_TYPES}
    ar = {k: (m_aro * t[k] / aro_sum if aro_sum > 0 else 0.0) for k in ARO_TYPES}
    return {
        "n_paraffin": sa["paraffins"] * n_share,
        "iso_paraffin": sa["paraffins"] * (1.0 - n_share),
        "cyclo_mono": sa["cycloparaffins"],
        "cyclo_di": sa["dicycloparaffins"],
        "alkylbenzene": ar["alkylbenzenes"],
        "indan_tetralin": ar["indans_tetralins"],
        "aromatic_di": ar["naphthalene"] + ar["naphthalenes"],
        "alkene_alkyne": m_ole,
        "unexplained_tricycloparaffins": sa["tricycloparaffins"],
        "unexplained_minor_aromatics": (ar["indenes"] + ar["acenaphthenes"] + ar["acenaphthylenes"]
                                        + ar["tricyclic_aromatics"]),
    }


def reconstruct_r3(row: Mapping, g1: Mapping, n_share: float, reps):
    """-> (DataFrame[inchikey, smiles, fraction, kind], fallback mass, {bucket: fallback mass})."""
    cols = ["inchikey", "smiles", "fraction", "kind"]
    empty = pd.DataFrame(columns=cols)
    cd = carbon_distribution(row)
    olefins = row.get("olefins_vol_pct", np.nan)
    cm = class_masses(g1, olefins, n_share) if cd is not None else None
    if cd is None or cm is None:
        log.warning("reconstruct_r3: skipped (%s)", "fewer than 4 distillation points" if cd is None else "unusable Table G.1 row")
        return empty, 1.0, {"skipped": 1.0}

    out: dict[tuple[int, str], float] = {}
    detail: dict[str, float] = {}
    for kind in SPREAD_DEFAULT + SPREAD_NEAREST:
        mass = cm[kind]
        if mass <= 0:
            continue
        for n, f in cd.items():
            if f <= 0:
                continue
            i = reps.get(kind, n) if kind in SPREAD_DEFAULT else reps.get(kind, n, max_gap=99)
            if i is None:
                detail[f"no_rep_{kind}"] = detail.get(f"no_rep_{kind}", 0.0) + mass * f
            else:
                out[(i, kind)] = out.get((i, kind), 0.0) + mass * f
    for b in UNEXPLAINED:
        if cm[b] > 0:
            detail[b] = detail.get(b, 0.0) + cm[b]
    fallback = float(sum(detail.values()))
    rows = [dict(inchikey=reps.ik[i], smiles=reps.smiles[i], fraction=v, kind=k) for (i, k), v in out.items()]
    return (pd.DataFrame(rows, columns=cols) if rows else empty), fallback, detail
