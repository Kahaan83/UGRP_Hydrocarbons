"""Loader: ThermoML extractor output for the nine unlisted J. Chem. Eng. Data 2020 files -> measurements_dedup schema.

Input : results/v5/thermoml_2020/thermoml_hc_2020.csv (written by fuelblend.data.thermoml_extract; native ThermoML units)
Output: data/real_fuels/v5/thermoml2020_lab.csv, data/real_fuels/v5/thermoml2020_build_log.json

Rules
* Only mole-fraction and pure rows are read; any other composition scale is counted and dropped.
* A component with fraction 0 is removed.  A "binary" row left with one component is the authors' pure end point: kind = pure.
* Units: Pa*s -> mPa*s, m2/s -> cSt, kg/m3 unchanged, kPa unchanged.
* Kinematic viscosity is also written as dynamic viscosity where the same file has a density at the same composition and
  temperature (dyn = kin x density); otherwise only the kinematic row is kept.
* Repeats at one (property, system, composition, T) collapse to the median; conflict = spread above 5 %.

    PYTHONPATH=src python -m fuelblend.data.loaders_thermoml2020
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.data.loaders_boehm2022 import COLS, system_id

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "results" / "v5" / "thermoml_2020" / "thermoml_hc_2020.csv"
OUT = ROOT / "data" / "real_fuels" / "v5" / "thermoml2020_lab.csv"
LOG = ROOT / "data" / "real_fuels" / "v5" / "thermoml2020_build_log.json"
SOURCE = "thermoml_jced_2020_unlisted"
PROPS = {"Viscosity, Pa*s": ("viscosity_dyn_mPas", 1e3), "Kinematic viscosity, m2/s": ("viscosity_kin_cSt", 1e6),
         "Mass density, kg/m3": ("density_kg_m3", 1.0), "Vapor or sublimation pressure, kPa": ("vapor_pressure_kPa", 1.0)}


def _ident(inchi: str, cache: dict) -> tuple[str, float]:
    if inchi not in cache:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        m = Chem.MolFromInchi(inchi)
        cache[inchi] = (Chem.MolToInchiKey(m), float(Descriptors.MolWt(m)))
    return cache[inchi]


def build(path: Path = RAW, log: dict | None = None) -> pd.DataFrame:
    log = {} if log is None else log
    raw = pd.read_csv(path)
    cache: dict = {}
    recs, dropped = [], {"unsupported composition scale": 0, "unsupported property": 0, "bad composition": 0}
    for r in raw.itertuples():
        if r.prop not in PROPS:
            dropped["unsupported property"] += 1
            continue
        if r.comp_type not in ("mole fraction", "pure"):
            dropped["unsupported composition scale"] += 1
            continue
        x = np.array(json.loads(r.comp_values), float)
        ids = [_ident(i, cache) for i in json.loads(r.inchis)]
        if len(x) != len(ids) or not np.all(np.isfinite(x)) or (x < 0).any() or abs(x.sum() - 1) > 5e-3:
            dropped["bad composition"] += 1
            continue
        keep = x > 0
        x, ids = x[keep] / x[keep].sum(), [i for i, k in zip(ids, keep) if k]
        order = np.argsort([i[0] for i in ids])
        iks = [ids[j][0] for j in order]
        x = x[order]
        w = x * np.array([ids[j][1] for j in order])
        w = w / w.sum()
        prop, f = PROPS[r.prop]
        recs.append(dict(property=prop, value=float(r.value) * f, T_K=float(r.T_K), P_kPa=r.P_kPa, ref=r.file.removesuffix(".xml"),
                         system_id=system_id(iks), kind={1: "pure", 2: "binary"}.get(len(iks), "multi"), n_comp=len(iks),
                         comp_inchikeys=json.dumps(iks), xs=json.dumps([round(float(v), 6) for v in x]),
                         ws=json.dumps([round(float(v), 6) for v in w]), x_key=",".join(f"{v:.3f}" for v in x)))
    d = pd.DataFrame(recs)
    d["T_key"] = (d["T_K"] * 2).round() / 2
    # dynamic viscosity from kinematic x density of the same file, composition and temperature
    key = ["ref", "system_id", "x_key", "T_key"]
    rho = d[d.property == "density_kg_m3"].groupby(key).value.median().rename("rho")
    kin = d[d.property == "viscosity_kin_cSt"].join(rho, on=key)
    made = kin[kin.rho.notna()].assign(property="viscosity_dyn_mPas", value=lambda t: t.value * t.rho / 1000.0).drop(columns="rho")
    log.update(rows_read=int(len(raw)), dropped=dropped, kinematic_rows=int(len(kin)), dynamic_from_kinematic=int(len(made)))
    d = pd.concat([d, made], ignore_index=True)
    rows = []
    for (prop, sid, xk, tk), g in d.groupby(["property", "system_id", "x_key", "T_key"], sort=True):
        v = g["value"].to_numpy()
        med = float(np.median(v))
        q25, q75 = np.percentile(v, [25, 75])
        conflict = bool(len(v) >= 2 and (np.ptp(v) if len(v) == 2 else q75 - q25) > 0.05 * abs(med))
        f = g.iloc[0]
        rows.append(dict(property=prop, system_id=sid, x_key=xk, T_key=tk, P_key=101.0, T_K=f["T_K"], P_kPa=f["P_kPa"],
                         kind=f["kind"], n_comp=f["n_comp"], is_hc=True, comp_inchikeys=f["comp_inchikeys"], x_mole=f["xs"],
                         w_mass=f["ws"], fuel_id=np.nan, ref=";".join(sorted(g["ref"].unique())), composition_level="L3",
                         comp_coverage=1.0, comp_basis="mole", value=med, vmin=float(v.min()), vmax=float(v.max()),
                         q25=float(q25), q75=float(q75), n_obs=len(v), n_refs=g["ref"].nunique(), n_sources=1, sources=SOURCE,
                         n_outlier=0, conflict=conflict))
    out = pd.DataFrame(rows, columns=COLS)
    log.update(rows_out=int(len(out)), by_property_kind={f"{p}|{k}": int(n) for (p, k), n in out.groupby(["property", "kind"]).size().items()},
               mixture_viscosity_systems=int(out[(out.kind != "pure") & (out.property == "viscosity_dyn_mPas")].system_id.nunique()),
               conflicts=int(out.conflict.sum()))
    return out


def main() -> None:
    log: dict = {}
    out = build(log=log)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    LOG.write_text(json.dumps(log, indent=1), encoding="utf-8")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
