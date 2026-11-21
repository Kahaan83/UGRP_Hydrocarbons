"""Shared helpers for v6 Track A (freezing).  Reuses v5 code paths; nothing from v5 is reimplemented or modified."""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results/v6/freezing"
V5, V6 = ROOT / "data/real_fuels/v5", ROOT / "data/real_fuels/v6"


def load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def dataset() -> pd.DataFrame:
    return pd.read_csv(OUT / "dataset.csv")


def entities(ds: pd.DataFrame) -> dict[str, dict]:
    """entity -> dict(bins[cls,n_carbon,mass], mapped[ik,smiles,frac], unmapped).  v5 entities come from s5.build_entities();
    new AJFTD fuels (comp_source == 'v6') are mapped with the same G1 map_fuel and R3 representatives."""
    s5 = load("s5_v6", "scripts/v5/s5_freezing_nalkane.py")
    g1 = load("g1_v6", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    ents = s5.build_entities()
    new = ds[ds.comp_source == "v6"]
    if len(new):
        from fuelblend.baselines.flash_unifac import _lookup
        from fuelblend.data.reconstruction_r3 import SubReps
        from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
        reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
        prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
        reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
        comp = pd.read_csv(V6 / "ajftd_gcxgc_composition.csv", dtype={"fuel_id": str})
        for r in new.itertuples():
            b = comp[comp.fuel_id == r.fuel_id]
            agg, unm, _ = g1.map_fuel(b, reps)
            tot = sum(agg.values())
            idx = sorted(agg)
            ents[r.entity] = dict(
                source="AJFTD", fid=r.fuel_id, label=r.fuel_id, unmapped=float(unm),
                bins=pd.DataFrame(dict(cls=[s5.cls3(c) for c in b.hc_class], n_carbon=b.n_carbon.astype(int), mass=b.mass_pct.astype(float))),
                mapped=pd.DataFrame(dict(ik=[reps.ik[i] for i in idx], smiles=[reps.smiles[i] for i in idx], frac=[agg[i] / tot for i in idx])),
                raw_bins=b.rename(columns={"hc_class": "hc_class"}))
    return {e: ents[e] for e in ds.entity}


def unifac_designated(ent: dict, s5, g1, CR, fu, tdes) -> float:
    """Baseline (b): v5 UNIFAC-SvL, designated constants (delta_freeze_t2_unifac.py 'unifac_desig' path, unchanged)."""
    d5 = load("d5_v6", "scripts/v5/base/delta_freeze_t2_unifac.py")
    x, iks, smi, _ = d5.s5_inputs(ent, g1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gam = fu.make_unifac_gamma_fn_vec(smi)
        v, _ = CR.freezing_C(x, iks, smi, tdes, model="svl_unifac", gamma=gam)
    return float(v) if np.isfinite(v) else np.nan
