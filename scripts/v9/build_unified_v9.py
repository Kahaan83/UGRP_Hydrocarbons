"""Build data/unified_v9/: v5 unified rows plus verified v9 sources, re-deduplicated.  data/unified/ is read only.

New sources entering v9: ThermoML 2020 files not in v5 (loaders_thermoml_v9) and Boehm 2022 neat real fuels
(loaders_cold_visc_v9).  Edwards 2017 and pure_v9 drafts stay out (unverifiable, results/v9/verification.md).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.loaders_cold_visc_v9 import load_cold_visc_v9  # noqa: E402
from fuelblend.data.loaders_thermoml_v9 import load_thermoml_v9  # noqa: E402

V5 = ROOT / "data" / "unified"
OUT = ROOT / "data" / "unified_v9"
KEY = ["property", "system_id", "x_key", "T_key", "P_key"]
ABS_TOL = {"flash_point_C": 2.0, "pour_point_C": 3.0, "cloud_point_C": 3.0, "freezing_point_C": 1.0}
ID_MAP_COLS = ["source_id", "source_sample", "canonical_system_id", "match_key"]
_LABELS = None


def _digits(label):
    m = re.search(r"(\d{4,5})", str(label))
    return m.group(1) if m else None


def _gcxgc_ids():
    f = pd.read_csv(ROOT / "data/real_fuels/v5/ajftd_gcxgc_fuels.csv", usecols=["fuel_id"])
    return {f"fuel:{i}" for i in f.fuel_id}


def _label_table():
    """POSF number -> canonical AJFTD fuel id (prefer a fuel with GCxGC composition, else lowest id)."""
    det = pd.read_csv(ROOT / "data/raw/ajftd/ajftd_fuel_details.csv", usecols=["fuel_id", "fuel_label"],
                      low_memory=False).drop_duplicates()
    comp = {int(i.split(":")[1]) for i in _gcxgc_ids()}
    out = {}
    for fid, lab in zip(det.fuel_id, det.fuel_label):
        d = _digits(lab)
        if d:
            out.setdefault(d, []).append(int(fid))
    return {d: f"fuel:{min(ids, key=lambda i: (i not in comp, i))}" for d, ids in out.items()}


def canonical_id(source_id: str, sample: str) -> str:
    """Canonical system_id for a real-fuel sample: POSF number matched against AJFTD fuel labels."""
    global _LABELS
    if _LABELS is None:
        _LABELS = _label_table()
    d = _digits(sample)
    if d and d in _LABELS and "POSF" in str(sample).upper():
        return _LABELS[d]
    return f"fuel:{source_id}:{sample}"


def _stats(vals):
    return dict(value=float(np.median(vals)), vmin=float(vals.min()), vmax=float(vals.max()),
                q25=float(np.percentile(vals, 25)), q75=float(np.percentile(vals, 75)))


def _conflict(prop, vals):
    spread = float(np.percentile(vals, 75) - np.percentile(vals, 25)) if len(vals) > 2 else float(vals.max() - vals.min())
    if prop in ABS_TOL:
        return spread > ABS_TOL[prop]
    med = abs(float(np.median(vals)))
    return spread > 0.05 * med if med > 0 else spread > 0


def _real_fuel_keys(df):
    """v5 real-fuel key convention: x_key = numeric fuel id, T rounded to 1 K, P_key 101.0."""
    df = df.copy()
    rf = df.kind == "real_fuel"
    num = df.system_id.str.extract(r"^fuel:(\d+)$")[0]
    m = rf & num.notna()
    df.loc[m, "x_key"] = num[m]
    df.loc[m, "fuel_id"] = num[m]
    df.loc[rf, "T_key"] = df.loc[rf, "T_K"].round(0).astype(float).astype(str)
    df.loc[rf, "P_key"] = "101.0"
    return df


def build():
    v5 = pd.read_csv(V5 / "measurements_dedup.csv", low_memory=False)
    raw5 = pd.read_csv(V5 / "measurements.csv", low_memory=False)
    for c in ("T_key", "P_key", "x_key"):
        v5[c] = v5[c].astype(str)

    cold = load_cold_visc_v9()
    idmap = []
    for s in sorted(cold.sample_id.unique()):
        cid = canonical_id("boehm2022", s)
        idmap.append(dict(source_id="boehm2022", source_sample=s, canonical_system_id=cid,
                          match_key="none" if cid.startswith("fuel:boehm2022") else "posf"))
    idmap = pd.DataFrame(idmap, columns=ID_MAP_COLS)
    cold["system_id"] = cold.sample_id.map(dict(zip(idmap.source_sample, idmap.canonical_system_id)))
    cold = _real_fuel_keys(cold)
    cold["source"] = "boehm2022_v9"
    cold["no_composition"] = ~cold.system_id.isin(_gcxgc_ids())

    tml = load_thermoml_v9()
    tml["source"] = "thermoml_v9"
    new = pd.concat([cold, tml], ignore_index=True)
    for c in ("T_key", "P_key", "x_key"):
        new[c] = new[c].astype(str)

    g_rows = []
    for _, grp in new.groupby(KEY, sort=False):
        vals = grp.value.to_numpy(float)
        base = grp.iloc[0].to_dict()
        base.update(_stats(vals))
        base.update(n_obs=len(vals), n_refs=grp.ref.nunique(), n_sources=1, sources=base["source"], n_outlier=0,
                    conflict=_conflict(base["property"], vals))
        g_rows.append(base)
    newg = pd.DataFrame(g_rows)

    idx = {k: i for i, k in zip(v5.index, v5[KEY].itertuples(index=False, name=None))}
    merged, exact_repeats, keep_new = 0, 0, []
    v5["v9_merged"] = False
    for r in newg.to_dict("records"):
        k = tuple(r[c] for c in KEY)
        if k not in idx:
            keep_new.append(r)
            continue
        i = idx[k]
        if "thermoml" in str(v5.at[i, "sources"]) and r["source"] == "thermoml_v9"                 and abs(r["value"] - v5.at[i, "value"]) <= 1e-9 * max(1.0, abs(r["value"])):
            exact_repeats += 1  # same source, same value: exact repeat, dropped (AGENTS.md dedup rule)
            continue
        pooled = np.array([v5.at[i, "value"], r["value"]], float)
        v5.at[i, "value"] = float(np.median(pooled))
        v5.at[i, "vmin"] = min(v5.at[i, "vmin"], r["vmin"])
        v5.at[i, "vmax"] = max(v5.at[i, "vmax"], r["vmax"])
        v5.at[i, "n_obs"] = int(v5.at[i, "n_obs"]) + int(r["n_obs"])
        v5.at[i, "n_refs"] = int(v5.at[i, "n_refs"]) + 1
        v5.at[i, "n_sources"] = int(v5.at[i, "n_sources"]) + 1
        v5.at[i, "sources"] = f"{v5.at[i, 'sources']}+{r['source']}"
        v5.at[i, "conflict"] = bool(v5.at[i, "conflict"]) or _conflict(r["property"], pooled)
        v5.at[i, "v9_merged"] = True
        merged += 1
    v5["v9_new"] = False
    add = pd.DataFrame(keep_new)
    add["v9_new"] = True
    add["v9_merged"] = False
    cols = list(v5.columns) + ["no_composition"]
    dd = pd.concat([v5, add.reindex(columns=cols)], ignore_index=True)
    dd["no_composition"] = dd["no_composition"].fillna(False).astype(bool)

    OUT.mkdir(parents=True, exist_ok=True)
    dd.to_csv(OUT / "measurements_dedup.csv", index=False)
    keep = [c for c in new.columns if c in raw5.columns]
    pd.concat([raw5, new[keep]], ignore_index=True).to_csv(OUT / "measurements.csv", index=False)
    idmap.to_csv(OUT / "fuel_identity_map.csv", index=False)
    _compounds_added(add)
    comp = pd.read_csv(ROOT / "data/real_fuels/v5/ajftd_gcxgc_composition.csv")
    comp["no_carbon_split"] = False
    comp.to_csv(OUT / "real_fuel_composition_v9.csv", index=False)
    log = dict(v5_dedup_rows=int(len(v5)), new_groups=int(len(add)), merged_into_v5=int(merged), exact_repeats_dropped=int(exact_repeats),
               new_by_source={k: int(v) for k, v in add.groupby("source").size().items()},
               new_real_fuel_groups_no_composition=int((add.kind.eq("real_fuel") & add.no_composition.fillna(False)).sum()),
               total_dedup_rows=int(len(dd)))
    (OUT / "build_log.json").write_text(json.dumps(log, indent=2))
    return log


def _compounds_added(add):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    known = set(pd.read_csv(V5 / "compounds.csv", usecols=["inchikey"]).inchikey)
    raw = pd.read_csv(ROOT / "data/raw/thermoml_v9/extracted.csv", usecols=["inchis"]).inchis.drop_duplicates()
    rows = {}
    for lst in raw:
        for inchi in json.loads(lst):
            m = Chem.MolFromInchi(inchi)
            if m is None:
                continue
            ik = Chem.MolToInchiKey(m)
            if ik in known or ik in rows:
                continue
            rows[ik] = dict(inchikey=ik, smiles=Chem.MolToSmiles(m), inchi=inchi,
                            formula=rdMolDescriptors.CalcMolFormula(m), mw=round(Descriptors.MolWt(m), 3),
                            is_hc=all(a.GetSymbol() in ("C", "H") for a in m.GetAtoms()))
    pd.DataFrame(list(rows.values())).to_csv(OUT / "compounds_v9_added.csv", index=False, quoting=1)


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
