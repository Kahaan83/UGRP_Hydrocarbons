"""v5 BASE Task 3: one target table for every T2 fuel (measured composition + measured properties).

Entities are the 155 fuels of results/v5/S3/corrected_2026-11-12_fuels.csv (AJFTD GCxGC 113, AV-33-22 38, FuelLib 4).
One row per (entity, property, T_K) holding the MEDIAN of the non-censored values; censored values (qualifier not null) are bounds
and are written to t2_censored_excluded.csv instead.  Nothing is fitted.

    .venv/Scripts/python.exe scripts/v5/base/build_t2_targets.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import _lookup  # noqa: E402
from fuelblend.data.reconstruction_r3 import SubReps  # noqa: E402
from fuelblend.eval.composition_route import kind_bins, map_bins  # noqa: E402
from fuelblend.eval.group_stats import CALC, MEAS, classify_method, flash_family  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402

V5 = ROOT / "data" / "real_fuels" / "v5"
S3 = ROOT / "results" / "v5" / "S3"
OUT = ROOT / "results" / "v5" / "BASE"
PRIORITY = {"av3322": 0, "fuellib": 1, "ajftd": 2}                    # tie-break for the fuel kept per double-count group
SRC_KEY = {"AV-33-22": "av3322", "FuelLib": "fuellib", "AJFTD": "ajftd"}
PROPS = ["flash_point_C", "freezing_point_C", "viscosity_kin_cSt", "density_kg_m3", "heat_of_combustion_net_MJ_kg",
         "smoke_point_mm", "vapor_pressure_kPa", "distillation_T10_C", "distillation_T50_C", "distillation_T90_C"]
T_DEPENDENT = {"viscosity_kin_cSt", "density_kg_m3", "vapor_pressure_kPa"}
FUELLIB_RENAME = {"net_heat_of_combustion_MJ_kg": "heat_of_combustion_net_MJ_kg"}


def _long(fuels: pd.DataFrame, log: dict) -> pd.DataFrame:
    parts = []
    # ---- AV-33-22
    a = pd.read_csv(V5 / "crc_av3322_properties.csv")
    a = a[a.property.isin(PROPS) & a.value.notna()]
    a = a.assign(entity="AV-33-22:" + a.fuel_id, source="av3322")[["entity", "source", "property", "value", "T_K", "method", "qualifier"]]
    parts.append(a)
    # ---- AJFTD
    ids = set(fuels[fuels.source == "AJFTD"].fuel_id.astype(int))
    j = pd.read_csv(V5 / "ajftd_fuels_long.csv", low_memory=False)
    j = j[j.fuel_id.isin(ids) & j.property.isin(PROPS) & j.value.notna()]
    drop_visc = (j.property == "viscosity_kin_cSt") & (j.T_source.astype(str) != "reported")
    drop_noT = j.property.isin(["viscosity_kin_cSt", "density_kg_m3"]) & j.T_K.isna()
    log["ajftd_viscosity_T_not_reported_skipped"] = int(drop_visc.sum())
    log["ajftd_T_dependent_without_T_skipped"] = int((drop_noT & ~drop_visc).sum())
    j = j[~(drop_visc | drop_noT)]
    j = j.assign(entity="AJFTD:" + j.fuel_id.astype(int).astype(str), source="ajftd")[["entity", "source", "property", "value", "T_K", "method", "qualifier"]]
    parts.append(j)
    # ---- FuelLib
    ids = set(fuels[fuels.source == "FuelLib"].fuel_id)
    f = pd.read_csv(V5 / "fuellib_properties.csv")
    f["property"] = f.property.replace(FUELLIB_RENAME)
    f = f[f.fuel_id.isin(ids) & f.property.isin(PROPS) & f.value.notna()].copy()
    f["T_K"] = f.T_C + 273.15
    noT = f.property.isin(T_DEPENDENT) & f.T_C.isna()
    log["fuellib_T_dependent_without_T_skipped"] = int(noT.sum())
    f = f[~noT]
    f = f.assign(entity="FuelLib:" + f.fuel_id, source="fuellib", qualifier=np.nan)[["entity", "source", "property", "value", "T_K", "method", "qualifier"]]
    parts.append(f)
    return pd.concat(parts, ignore_index=True)


def _mode(s: pd.Series) -> str:
    s = s.dropna().astype(str)
    return "" if s.empty else sorted(s.value_counts().sort_index().items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def _aggregate(raw: pd.DataFrame, fuels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = raw.copy()
    raw["T_K"] = raw.T_K.round(2)
    cens = raw[raw.qualifier.notna()].copy()
    ok = raw[raw.qualifier.isna()].copy()
    ok["Tk"] = ok.T_K.fillna(-1.0)
    cens["Tk"] = cens.T_K.fillna(-1.0)
    rows = []
    for (ent, prop, tk), g in ok.groupby(["entity", "property", "Tk"]):
        v = g.value.to_numpy(float)
        fam = ""
        val = float(np.median(v))
        if prop == "flash_point_C":
            fam = _mode(g.method.map(flash_family))
        elif prop == "heat_of_combustion_net_MJ_kg":
            cls = g.method.map(classify_method)
            pick = MEAS if (cls == MEAS).any() else CALC if (cls == CALC).any() else None   # measured > calculated > unclassified
            if pick:
                val = float(np.median(v[(cls == pick).to_numpy()]))
            fam = (pick or "UNCLASSIFIED").lower()
        nd = int(((cens.entity == ent) & (cens.property == prop) & (cens.Tk == tk)).sum())
        rows.append(dict(entity=ent, source=g.source.iloc[0], property=prop, T_K=np.nan if tk < 0 else tk, measured=val, n_values=len(v),
                         method="+".join(sorted(set(g.method.dropna().astype(str).str.replace(" ", "")))) or "unstated",
                         method_family=fam if fam else np.nan, qualifier_dropped=nd))
    t = pd.DataFrame(rows)
    fu = fuels.set_index("entity")
    t["fuel_id"] = t.entity.map(fu.fuel_id)
    t["label"] = t.entity.map(fu.label)
    t["group"] = t.entity.map(fu.group)
    t["already_read_locked"] = t.entity.map(fu.locked_ajftd.astype(bool) | fu.group_locked.astype(bool))
    t = t[["entity", "source", "fuel_id", "label", "group", "property", "T_K", "measured", "n_values", "method", "method_family",
           "qualifier_dropped", "already_read_locked"]].sort_values(["entity", "property", "T_K"]).reset_index(drop=True)
    return t, cens.drop(columns="Tk")


def _coverage(fuels: pd.DataFrame, reps) -> pd.DataFrame:
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    rows = []
    for f in fuels.itertuples():
        src = SRC_KEY[f.source]
        c = comps[src]
        fid = int(f.fuel_id) if src == "ajftd" else f.fuel_id
        sub = c[c.fuel_id == fid]
        if sub.empty:
            rows.append(dict(entity=f.entity, unmapped_fraction=1.0, n_components=0))
            continue
        m = map_bins(kind_bins(src, sub), reps)
        rows.append(dict(entity=f.entity, unmapped_fraction=m.unmapped_fraction, n_components=int(len(m.idx))))
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fuels = pd.read_csv(S3 / "corrected_2026-11-12_fuels.csv")
    log: dict = {}
    t, cens = _aggregate(_long(fuels, log), fuels)
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    cov = _coverage(fuels, reps)
    fl = fuels[["entity", "source", "fuel_id", "label", "group"]].merge(cov, on="entity")
    fl["source"] = fl.source.map(SRC_KEY)
    fl["low_coverage"] = fl.unmapped_fraction > 0.10                    # COVERAGE_MIN = 0.90 of mass must map
    fl["n_target_rows"] = fl.entity.map(t.groupby("entity").size()).fillna(0).astype(int)
    order = fl.assign(_a=(fl.n_components > 0).map({True: 0, False: 1}), _b=-fl.n_target_rows, _c=fl.source.map(PRIORITY)) \
              .sort_values(["group", "_a", "_b", "_c", "entity"])
    fl["chosen_for_group"] = fl.entity.isin(order.groupby("group").head(1).entity)
    fl = fl[["entity", "source", "fuel_id", "label", "group", "chosen_for_group", "unmapped_fraction", "low_coverage",
             "n_components", "n_target_rows"]]
    t.to_csv(OUT / "t2_targets.csv", index=False, float_format="%.10g")
    fl.to_csv(OUT / "t2_fuels.csv", index=False, float_format="%.10g")
    cens.to_csv(OUT / "t2_censored_excluded.csv", index=False)
    log.update(n_target_rows=len(t), n_entities=int(fuels.shape[0]), n_groups=int(fuels.group.nunique()),
               n_censored_excluded=int(len(cens)), n_low_coverage=int(fl.low_coverage.sum()),
               rows_by_property=t.groupby("property").size().to_dict())
    json.dump(log, open(OUT / "t2_build_log.json", "w"), indent=1, default=int)
    print(json.dumps(log, indent=1, default=int))


if __name__ == "__main__":
    main()
