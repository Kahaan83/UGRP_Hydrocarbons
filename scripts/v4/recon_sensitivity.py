"""v4 review fix 2: sensitivity of CRC 647 property predictions to the nominal saturates split (n, iso, cyclo).

DEV fuels only (locked fuels never enter).  32 Dirichlet draws, rng(20261111).dirichlet(nominal * 50), plus the nominal
split itself (draw = -1).  For every (fuel, draw) the r2 reconstruction is rebuilt with that split and four zero-parameter
physics predictions are scored: density at 15.56 C (ideal volume), ideal-SLE freezing (mole fractions), raw Refutas at
-20 C (kinematic, cSt), ideal Le Chatelier flash.  Output: results/v4/R2/sat_split_sensitivity.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "v4"))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / "v4" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


bulk = _load("eval_bulk_props_v4")
flash = _load("eval_flash_v4")
visc = _load("build_viscosity_leaderboard_v4")
rebuild = _load("rebuild_crc_reconstruction")

from fuelblend.baselines.flash_unifac import _lookup  # noqa: E402
from fuelblend.data.reconstruction_r2 import RealisticReps, nominal_sat_split  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.mixture.cold_flow import compound_meta  # noqa: E402

N_DRAWS, SEED, CONC = 32, 20261111, 50
OUT = ROOT / "results" / "v4" / "R2" / "sat_split_sensitivity.csv"
R2 = ROOT / "data" / "real_fuels" / "v4" / "compositions_crc647_r2.csv"


def comp_frame(out_rows) -> pd.DataFrame:
    """r2 rows (as build_compositions emits them) -> the frame layout the bulk evaluators expect (lump removed)."""
    c = pd.DataFrame(out_rows)
    c = c[~c.is_lumped & (c.inchikey != "LUMP_NO_REALISTIC_REP")].copy()
    c["fuel_id"] = c.system_id.str.replace("fuel:", "", regex=False)
    ex = []
    for smi in c.smiles:
        m, deut = bulk.norm_mol(smi)
        ex.append(dict(smiles_norm=Chem.MolToSmiles(m), inchikey_norm=Chem.MolToInchiKey(m), deuterated=deut,
                       cls=bulk.hc_class(m), n_carbon=sum(a.GetSymbol() == "C" for a in m.GetAtoms()),
                       mw=Descriptors.MolWt(m)))
    return pd.concat([c.reset_index(drop=True), pd.DataFrame(ex)], axis=1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--n-draws", type=int, default=N_DRAWS)
    a = ap.parse_args(argv)

    bulk.COMP = R2
    fuels, role, comp0, props, _ = bulk.load_inputs()
    dev = [f for f in fuels if role.get(f) != "locked"]
    props = props[props.fuel_id.isin(dev)]
    meas = {"density_15.56C": props[(props.property == "density_kg_m3") & np.isclose(props.T_C, 15.56)].set_index("fuel_id").value,
            "freezing": props[props.property == "freezing_point_C"].set_index("fuel_id").value}
    allp = pd.read_csv(rebuild.PROPS)
    kv = allp[(allp.property == "viscosity_kin_cSt") & (allp.provenance == "raw_file") & np.isclose(allp.T_C, -20.0)
              & allp.fuel_id.isin(dev)].groupby("fuel_id").value.mean()
    fl = allp[(allp.property == "flash_point_C") & (allp.provenance == "raw_file") & allp.fuel_id.isin(dev)
              & allp.fuel_id.str.startswith("CRC647")].groupby("fuel_id").value.mean()
    meas["viscosity_-20C_cSt"], meas["flash"] = kv, fl

    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = RealisticReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    _, cls, nc = compound_meta(ROOT / "data/unified/compounds.csv")
    mp = visc.Mapper(reg, cls, nc, prov)
    nominal = nominal_sat_split(rebuild.JETCLASS)
    rows_by_fuel = {f: r for f, r in rebuild.descriptor_rows(pd.read_csv(rebuild.PROPS)).items() if f in dev}

    rng = np.random.default_rng(SEED)
    splits = [(-1, nominal)] + [(i, tuple(rng.dirichlet(np.array(nominal) * CONC))) for i in range(a.n_draws)]
    pd_model = bulk.PureDensity(comp0)
    tmap = bulk.thermal_map_for(comp0)          # same compounds in every draw (reps do not depend on the split)
    recs = []
    for draw, sp in splits:
        out, *_ = rebuild.build_compositions(rows_by_fuel, sp, reps)
        comp = comp_frame(out)
        miss = set(comp.inchikey_norm) - set(tmap)
        if miss:
            tmap.update(bulk.thermal_map_for(comp[comp.inchikey_norm.isin(miss)]))
        raw = pd.DataFrame(out)
        for f in dev:
            g = comp[comp.fuel_id == f]
            if g.empty:
                continue
            w = bulk.mass_frac(g)
            base = dict(fuel_id=f, split_role=role.get(f, "dev"), draw=draw, sat_n=sp[0], sat_iso=sp[1], sat_cyclo=sp[2])
            if f in meas["density_15.56C"].index:
                rho = [pd_model.rho(ik, c, n, 288.71)[0] for ik, c, n in zip(g.inchikey_norm, g.cls, g.n_carbon)]
                recs.append(dict(base, property="density_15.56C", T_C=15.56, predicted=float(1.0 / np.sum(w / np.array(rho))),
                                 measured=float(meas["density_15.56C"][f])))
            if f in meas["freezing"].index:
                recs.append(dict(base, property="freezing", T_C=np.nan,
                                 predicted=bulk.ideal_sle_mole_C(w, g.inchikey_norm.tolist(), g.mw.to_numpy(float), tmap),
                                 measured=float(meas["freezing"][f])))
            if f in meas["flash"].index:
                smi = list(g.smiles)
                x = flash.mass_to_mole(smi, g.fraction.to_numpy(float))
                recs.append(dict(base, property="flash", T_C=np.nan, predicted=flash.run_detail(smi, x, False)["flash_C"],
                                 measured=float(meas["flash"][f])))
            if f in meas["viscosity_-20C_cSt"].index:
                idx, wv, _, _ = visc.fuel_composition(raw[raw.system_id == f"fuel:{f}"], mp, direct=True)
                recs.append(dict(base, property="viscosity_-20C_cSt", T_C=-20.0,
                                 predicted=visc.refutas_nu(prov, idx, wv, 253.15)[0], measured=float(meas["viscosity_-20C_cSt"][f])))
        print(f"draw {draw} split={tuple(round(v, 3) for v in sp)} done", flush=True)
    df = pd.DataFrame(recs, columns=["fuel_id", "split_role", "draw", "sat_n", "sat_iso", "sat_cyclo", "property", "T_C",
                                     "predicted", "measured"])
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(df.groupby(["property", df.draw == -1]).size())


if __name__ == "__main__":
    main()
