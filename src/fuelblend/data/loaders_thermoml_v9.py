"""Loader: ThermoML post-2020 v9 extractor output -> measurements_dedup schema.

Input : data/raw/thermoml_v9/extracted.csv (written by scripts/v9/fetch_thermoml.py)
Output: pd.DataFrame conforming to the project measurements_dedup schema (COLS).

Rules:
* Only mole-fraction and pure rows are read; any other composition scale is counted and dropped.
* A component with fraction 0 is removed. A "binary" row left with one component is kind = pure.
* Units: Pa*s -> mPa*s (x1e3), m2/s -> cSt (x1e6), kg/m3 unchanged, kPa unchanged.
* Kinematic viscosity is also converted to dynamic viscosity where the same file has a density
  at the same composition and temperature (dyn = kin * density / 1000.0).
* Canonical system_id: md5[:12] of sorted InChIKeys.
* Deduplicate on (property, system_id, x_key, T_key, P_key): median value; conflict spread > 5%.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.data.loaders_boehm2022 import COLS, system_id

ROOT = Path(__file__).resolve().parents[3]
RAW_PATH = ROOT / "data" / "raw" / "thermoml_v9" / "extracted.csv"
SOURCE_NAME = "thermoml_v9"

PROPS = {
    "Viscosity, Pa*s": ("viscosity_dyn_mPas", 1e3),
    "Kinematic viscosity, m2/s": ("viscosity_kin_cSt", 1e6),
    "Mass density, kg/m3": ("density_kg_m3", 1.0),
    "Vapor or sublimation pressure, kPa": ("vapor_pressure_kPa", 1.0),
}


def _ident(inchi: str, cache: dict) -> tuple[str, float]:
    if inchi not in cache:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        m = Chem.MolFromInchi(inchi)
        if m is None:
            cache[inchi] = ("", 0.0)
        else:
            cache[inchi] = (Chem.MolToInchiKey(m), float(Descriptors.MolWt(m)))
    return cache[inchi]


def load_thermoml_v9(path: Path = RAW_PATH, log: dict | None = None) -> pd.DataFrame:
    """Load and standardize post-2020 ThermoML extracted hydrocarbon measurements."""
    log = {} if log is None else log
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    raw = pd.read_csv(path)
    cache: dict = {}
    recs = []
    dropped = {"unsupported composition scale": 0, "unsupported property": 0, "bad composition": 0}

    for r in raw.itertuples():
        if r.prop not in PROPS:
            dropped["unsupported property"] += 1
            continue
        if r.comp_type not in ("mole fraction", "pure"):
            dropped["unsupported composition scale"] += 1
            continue
        
        # parse mole fractions and components
        x_vals = np.array(json.loads(r.comp_values), float)
        inchis = json.loads(r.inchis)
        ids = [_ident(i, cache) for i in inchis]
        
        # validate composition
        if len(x_vals) != len(ids) or not np.all(np.isfinite(x_vals)) or (x_vals < 0).any() or abs(x_vals.sum() - 1.0) > 5e-3:
            dropped["bad composition"] += 1
            continue
            
        keep = x_vals > 0
        x_vals = x_vals[keep] / x_vals[keep].sum()
        ids = [i for i, k in zip(ids, keep) if k]
        
        if not ids or any(not i[0] for i in ids):
            dropped["bad composition"] += 1
            continue

        order = np.argsort([i[0] for i in ids])
        iks = [ids[j][0] for j in order]
        x_vals = x_vals[order]
        w_vals = x_vals * np.array([ids[j][1] for j in order])
        w_vals = w_vals / w_vals.sum()
        
        prop, f = PROPS[r.prop]
        ref_name = str(r.file).removesuffix(".xml")
        p_val = float(r.P_kPa) if pd.notna(r.P_kPa) else 101.325

        recs.append({
            "property": prop,
            "value": float(r.value) * f,
            "T_K": float(r.T_K),
            "P_kPa": p_val,
            "ref": ref_name,
            "system_id": system_id(iks),
            "kind": {1: "pure", 2: "binary"}.get(len(iks), "multi"),
            "n_comp": len(iks),
            "comp_inchikeys": json.dumps(iks),
            "xs": json.dumps([round(float(v), 6) for v in x_vals]),
            "ws": json.dumps([round(float(v), 6) for v in w_vals]),
            "x_key": ",".join(f"{v:.3f}" for v in x_vals),
        })

    d = pd.DataFrame(recs)
    if d.empty:
        return pd.DataFrame(columns=COLS)

    d["T_key"] = (d["T_K"] * 2).round() / 2
    d["P_key"] = d["P_kPa"].apply(lambda p: 101.0 if (pd.isna(p) or p <= 150.0) else round(p / 100.0) * 100.0)

    # Dynamic viscosity from kinematic * density of the same file, composition and temperature
    key = ["ref", "system_id", "x_key", "T_key"]
    rho = d[d["property"] == "density_kg_m3"].groupby(key)["value"].median().rename("rho")
    kin = d[d["property"] == "viscosity_kin_cSt"].join(rho, on=key)
    made = kin[kin["rho"].notna()].assign(
        property="viscosity_dyn_mPas",
        value=lambda t: t["value"] * t["rho"] / 1000.0
    ).drop(columns="rho")

    log.update(
        rows_read=int(len(raw)),
        dropped=dropped,
        kinematic_rows=int(len(kin)),
        dynamic_from_kinematic=int(len(made)),
    )
    d = pd.concat([d, made], ignore_index=True)

    rows = []
    for (prop, sid, xk, tk, pk), g in d.groupby(["property", "system_id", "x_key", "T_key", "P_key"], sort=True):
        v = g["value"].to_numpy()
        med = float(np.median(v))
        q25, q75 = np.percentile(v, [25, 75])
        conflict = bool(len(v) >= 2 and (np.ptp(v) if len(v) == 2 else q75 - q25) > 0.05 * abs(med))
        f = g.iloc[0]
        rows.append({
            "property": prop,
            "system_id": sid,
            "x_key": xk,
            "T_key": tk,
            "P_key": pk,
            "T_K": f["T_K"],
            "P_kPa": f["P_kPa"],
            "kind": f["kind"],
            "n_comp": f["n_comp"],
            "is_hc": True,
            "comp_inchikeys": f["comp_inchikeys"],
            "x_mole": f["xs"],
            "w_mass": f["ws"],
            "fuel_id": np.nan,
            "ref": ";".join(sorted(g["ref"].unique())),
            "composition_level": "L3",
            "comp_coverage": 1.0,
            "comp_basis": "mole",
            "value": med,
            "vmin": float(v.min()),
            "vmax": float(v.max()),
            "q25": float(q25),
            "q75": float(q75),
            "n_obs": len(v),
            "n_refs": g["ref"].nunique(),
            "n_sources": 1,
            "sources": SOURCE_NAME,
            "n_outlier": 0,
            "conflict": conflict,
        })

    out = pd.DataFrame(rows, columns=COLS)
    log.update(
        rows_out=int(len(out)),
        by_property_kind={f"{p}|{k}": int(n) for (p, k), n in out.groupby(["property", "kind"]).size().items()},
        mixture_viscosity_systems=int(out[(out["kind"] != "pure") & (out["property"] == "viscosity_dyn_mPas")]["system_id"].nunique()),
        conflicts=int(out["conflict"].sum()),
    )
    return out
