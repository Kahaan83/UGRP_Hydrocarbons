"""v5 representatives: SubReps plus audit-proposal-8 candidates that the registry already supports with verified constants.

A candidate is admitted only with: Antoine from refit_measured or webbook; LFL that is not a Jones-rule estimate (measured, or tabulated in the
v5 overlay); a measured pure viscosity fit; a valid provider density.  An admitted candidate replaces its own exact (subkind, n_carbon) bin;
all other bins return exactly what SubReps returns.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from fuelblend.data.reconstruction_r3 import SubReps, subkind

ROOT = Path(__file__).resolve().parents[3]
CANDIDATES = {                                     # name -> SMILES, in order of preference within a bin
    "1,2,4-trimethylbenzene": "Cc1ccc(C)c(C)c1",
    "n-propylbenzene": "CCCc1ccccc1",
    "1-methylnaphthalene": "Cc1cccc2ccccc12",
    "2,6-dimethylnaphthalene": "Cc1ccc2cc(C)ccc2c1",
    "perhydrophenanthrene": "C1CCC2C(C1)CCC1CCCCC21",
    "perhydrofluorene": "C1CCC2C(C1)CC1CCCCC12",
    "1-methyltetralin": "CC1CCCc2ccccc21",
    "n-butyltetralin": "CCCCC1CCCc2ccccc21",
    "2-methyldecalin": "CC1CCC2CCCCC2C1",
}
REPLACED_BINS: set[tuple[str, int]] = {("alkylbenzene", 9)}     # bins whose admitted candidate replaces the SubReps choice
GOOD_ANTOINE = ("refit_measured", "webbook")
TABULATED_LFL = ("iec60079_20_1_table", "cooptima_compiled", "zabetakis627_measured")


def _can(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))


def candidate_table(reg, prov, compounds_csv: Path | None = None, lfl_csv: Path | None = None) -> pd.DataFrame:
    from fuelblend.mixture import common
    c = pd.read_csv(compounds_csv or ROOT / "data/unified/compounds.csv")
    c["can"] = [_can(s) for s in c.smiles]
    lfl = pd.read_csv(lfl_csv or ROOT / "data/processed/v5/lfl_v5.csv").set_index("inchikey")
    fits = common._get_measured_fits()[1]
    ik_index = {k: i for i, k in enumerate(reg.ik)}
    rows = []
    for name, smi in CANDIDATES.items():
        can = _can(smi)
        m = Chem.MolFromSmiles(can)
        rec = dict(name=name, smiles=can, kind=subkind(can), n_carbon=int(sum(a.GetSymbol() == "C" for a in m.GetAtoms())))
        hit = c[c.can == can]
        if hit.empty:
            rec.update(has_antoine=False, has_lfl=False, has_viscosity=False, has_density=False, has_fusion=False, admitted=False, reason="not in registry")
            rows.append(rec)
            continue
        h = hit.iloc[0]
        ik = h.inchikey
        abc = np.array([h.antoine_A, h.antoine_B, h.antoine_C], float)
        has_ant = bool(np.isfinite(abc).all() and str(h.antoine_source) in GOOD_ANTOINE)
        has_lfl = bool(ik in lfl.index and str(lfl.lfl_source_v5[ik]) in TABULATED_LFL)
        i = ik_index.get(ik)
        if i is not None:
            prov.ensure([i])
        has_rho = bool(i is not None and prov.valid[i])
        has_visc = bool(ik in fits) and has_rho
        rec.update(has_antoine=has_ant, has_lfl=has_lfl, has_viscosity=has_visc, has_density=has_rho,
                   has_fusion=False)                                      # fusion is not in compounds.csv; not required for admission
        miss = [k for k, ok in (("antoine", has_ant), ("lfl", has_lfl), ("viscosity", has_visc), ("density", has_rho)) if not ok]
        if rec["kind"] is None:
            miss.append("no SubReps subkind for this structure")
        rec["admitted"] = not miss
        rec["reason"] = "admitted" if not miss else "missing: " + ", ".join(miss)
        rows.append(rec)
    return pd.DataFrame(rows, columns=["name", "smiles", "kind", "n_carbon", "has_antoine", "has_lfl", "has_viscosity", "has_density",
                                       "has_fusion", "admitted", "reason"])


class SubRepsV5(SubReps):
    """SubReps with admitted candidates in REPLACED_BINS.  A replaced bin without an admitted candidate falls back to SubReps."""

    def __init__(self, reg, provider, vp_lookup, nc=None, table: pd.DataFrame | None = None):
        super().__init__(reg, provider, vp_lookup, nc)
        t = table if table is not None else candidate_table(reg, provider)
        can_index = {_can(s): i for i, s in enumerate(self.smiles)}
        self.replaced: dict[tuple[str, int], int] = {}
        for r in t[t.admitted.astype(bool)].itertuples():
            b = (r.kind, int(r.n_carbon))
            if b in REPLACED_BINS and b not in self.replaced and r.smiles in can_index:
                self.replaced[b] = can_index[r.smiles]

    def get(self, kind: str, n_carbon: int, max_gap: int = 3):
        b = (kind, int(n_carbon))
        if b in self.replaced:
            return self.replaced[b]
        return super().get(kind, n_carbon, max_gap)
