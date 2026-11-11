"""Boehm et al. 2022 (Front. Energy Res., doi 10.3389/fenrg.2022.1074699, CC BY) cold-temperature viscosity loader.

Reads ``data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx`` (sheet ``This Work`` ONLY; the ``Literature``,
``Relative Difference`` and ``Extrapolated_*`` sheets are other people's / derived numbers and are ignored) with
openpyxl and writes ``data/real_fuels/v5/boehm2022_viscosity.csv`` in the ``measurements_dedup.csv`` schema
(one row per property/system/composition/T), so it can be concatenated with that table.

Measured at 15, 0, -10, -20, -40 C: kinematic viscosity (cSt) and density (g/cm3).  Properties written:
``viscosity_kin_cSt`` (as measured), ``density_kg_m3`` (as measured), ``viscosity_dyn_mPas`` = kin x measured
density of the SAME sample at the SAME T (cSt * g/cm3 = mPa s).  ``scripts/v4/eval_lab_known_composition.py``
consumes ``viscosity_dyn_mPas`` rows with ``kind in (binary, multi)``, ``is_hc``, ``conflict == False``.

Composition: the sheet gives VOLUME percent at preparation (pipetted volumes, ambient).  Volume -> mass/mole uses the
lab's own measured 15 C pure-component densities (one composition per blend, not one per measurement T).  Nested
blends ("butylcyclohexane/heptylbenzene 50/50") are expanded to pure leaves assuming ideal volume additivity of the
sub-blend.  Any blend containing a real fuel / blend stock (A-1, A-3, C-1, C-8, C-9, HRJ, POSF*) is dropped.
"""
from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "data" / "raw" / "boehm2022" / "HEAT_LAB_ViscosityData_2022.xlsx"
OUT = REPO / "data" / "real_fuels" / "v5" / "boehm2022_viscosity.csv"
SHEET = "This Work"
REF = "boehm2022_this_work"
T_REF_DENSITY_C = 15

SMILES = {
    "n-pentane": "CCCCC", "n-hexane": "CCCCCC", "n-octane": "CCCCCCCC",
    "o-xylene": "Cc1ccccc1C", "m-xylene": "Cc1cccc(C)c1", "ethylbenzene": "CCc1ccccc1",
    "methylcyclohexane": "CC1CCCCC1", "ethylcyclohexane": "CCC1CCCCC1", "butylcyclohexane": "CCCCC1CCCCC1",
    "2,3-dimethylpentane": "CCC(C)C(C)C", "2,4-dimethylpentane": "CC(C)CC(C)C", "2-methylpentane": "CCCC(C)C",
    "trans-decalin": "C1CC[C@@H]2CCCC[C@H]2C1",
    "2,2,4-trimethylpentane": "CC(C)CC(C)(C)C", "2,2,4,6,6-pentamethylheptane": "CC(C)(C)CC(C)CC(C)(C)C",
    "2,2,4,4,6,8,8-heptamethylnonane": "CC(C)(C)CC(C)(C)CC(C)CC(C)(C)C",
    "hexylbenzene": "CCCCCCc1ccccc1", "heptylbenzene": "CCCCCCCc1ccccc1",
}
FUEL_RE = re.compile(r"POSF|^A-\d|^C-\d|HRJ", re.I)
COLS = ["property", "system_id", "x_key", "T_key", "P_key", "T_K", "P_kPa", "kind", "n_comp", "is_hc", "comp_inchikeys",
        "x_mole", "w_mass", "fuel_id", "ref", "composition_level", "comp_coverage", "comp_basis", "value", "vmin", "vmax",
        "q25", "q75", "n_obs", "n_refs", "n_sources", "sources", "n_outlier", "conflict"]


def system_id(iks) -> str:
    return hashlib.md5("|".join(sorted(iks)).encode()).hexdigest()[:12]


def _chem():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    out = {}
    for name, smi in SMILES.items():
        m = Chem.MolFromSmiles(smi)
        out[name] = (Chem.MolToInchiKey(m), float(Descriptors.MolWt(m)), smi)
    return out


def norm_name(s: str) -> str:
    s = str(s).strip().replace("\n", "")
    s = re.sub(r"\s*\(.*\)$", "", s)
    return {"iso-octane": "2,2,4-trimethylpentane"}.get(s, s)


def parse_component(label) -> list[tuple[str, float]]:
    """Header label -> [(leaf name, volume fraction within this component)].  'a/b 50/50' -> two leaves."""
    s = str(label).strip().replace("\n", "")
    if "/" in s and not re.search(r"POSF", s):
        body = re.sub(r"\s*\d+/\d+\s*$", "", s)
        parts = [norm_name(p) for p in body.split("/")]
        return [(p, 1.0 / len(parts)) for p in parts]
    return [(norm_name(s), 1.0)]


def read_blocks(path: Path = RAW):
    """-> list of (label_tuple, [data rows]); data row = (fractions[0..2], {T_C: (kin, rho)})."""
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[SHEET]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    temps = None
    blocks, cur = [], None
    for r in rows:
        r = list(r) + [None] * (16 - len(r))
        if isinstance(r[2], str) and isinstance(r[5], str) and r[5].startswith("Density at"):
            temps = [int(re.search(r"(-?\d+)C", r[k]).group(1)) for k in (5, 7, 9, 11, 13)]
            cur = {"labels": [r[k] for k in (2, 3, 4) if r[k] is not None], "rows": []}
            blocks.append(cur)
        elif isinstance(r[2], str) and cur is not None and all(c is None for c in r[5:15]):
            cur = {"labels": [r[k] for k in (2, 3, 4) if r[k] is not None], "rows": []}   # sub-header w/o T columns
            blocks.append(cur)
        elif cur is not None and isinstance(r[2], (int, float)):
            n = len(cur["labels"])
            fr = [float(x) for x in r[2:2 + n]]
            meas = {}
            for j, T in enumerate(temps):
                rho, kin = r[5 + 2 * j], r[6 + 2 * j]
                rho = float(rho) if isinstance(rho, (int, float)) else None
                kin = float(kin) if isinstance(kin, (int, float)) else None
                meas[T] = (kin, rho)
            cur["rows"].append((fr, meas))
    return blocks


def build(path: Path = RAW, log: dict | None = None) -> pd.DataFrame:
    log = {} if log is None else log
    chem = _chem()
    blocks = read_blocks(path)
    # 15 C pure densities (measured here), median over every block where a leaf appears as a 100 % single compound
    dens15 = defaultdict(list)
    for b in blocks:
        for fr, meas in b["rows"]:
            for k, f in enumerate(fr):
                if abs(f - 100) < 1e-9 and sum(fr) == 100:
                    leaves = parse_component(b["labels"][k])
                    rho = meas[T_REF_DENSITY_C][1]
                    if len(leaves) == 1 and rho and leaves[0][0] in chem:
                        dens15[leaves[0][0]].append(rho)
    rho15 = {k: float(np.median(v)) for k, v in dens15.items()}
    log["pure_density_15C_g_cm3"] = rho15

    dropped_fuel, dropped_unknown, dropped_no_density = set(), set(), set()
    n_rows_in = n_rows_kept = 0
    recs = []
    for b in blocks:
        comps = [parse_component(l) for l in b["labels"]]
        is_fuel = [bool(FUEL_RE.search(str(l).strip())) for l in b["labels"]]
        for fr, meas in b["rows"]:
            n_rows_in += 1
            if any(f > 0 and fu for f, fu in zip(fr, is_fuel)):
                dropped_fuel |= {str(l).strip() for l, f, fu in zip(b["labels"], fr, is_fuel) if fu and f > 0}
                continue
            names = {n for f, c in zip(fr, comps) if f > 0 for n, _ in c}
            unk = {n for n in names if n not in chem}
            if unk:
                dropped_unknown |= unk
                continue
            if any(n not in rho15 for n in names):
                dropped_no_density |= {n for n in names if n not in rho15}
                continue
            vol = defaultdict(float)
            for f, c in zip(fr, comps):
                for n, sub in c:
                    vol[n] += f / 100.0 * sub
            vol = {n: v for n, v in vol.items() if v > 1e-12}
            if abs(sum(vol.values()) - 1) > 1e-6:
                log.setdefault("bad_fraction_sum_rows", 0)
                log["bad_fraction_sum_rows"] += 1
                continue
            n_rows_kept += 1
            mass = {n: v * rho15[n] for n, v in vol.items()}
            tm = sum(mass.values())
            w = {n: m / tm for n, m in mass.items()}
            mol = {n: w[n] / chem[n][1] for n in w}
            tmol = sum(mol.values())
            x = {n: m / tmol for n, m in mol.items()}
            order = sorted(vol, key=lambda n: chem[n][0])
            iks = [chem[n][0] for n in order]
            base = dict(system_id=system_id(iks), kind={1: "pure", 2: "binary"}.get(len(iks), "multi"), n_comp=len(iks),
                        comp_inchikeys=json.dumps(iks), x=[round(x[n], 6) for n in order], w=[round(w[n], 6) for n in order])
            for T_C, (kin, rho) in meas.items():
                T_K = T_C + 273.15
                if kin is not None:
                    recs.append(dict(base, property="viscosity_kin_cSt", T_K=T_K, value=kin))
                    if rho is not None:
                        recs.append(dict(base, property="viscosity_dyn_mPas", T_K=T_K, value=kin * rho))
                    else:
                        log.setdefault("dyn_dropped_no_sample_density", 0)
                        log["dyn_dropped_no_sample_density"] += 1
                if rho is not None:
                    recs.append(dict(base, property="density_kg_m3", T_K=T_K, value=rho * 1000.0))
    log.update(dropped_fuel_or_blend_stock_components=sorted(dropped_fuel), dropped_unresolved_components=sorted(dropped_unknown),
               dropped_components_without_15C_density=sorted(dropped_no_density), rows_read=n_rows_in, rows_kept=n_rows_kept)
    d = pd.DataFrame(recs)
    d["T_key"] = (d["T_K"] * 2).round() / 2
    d["x_key"] = d["x"].map(lambda v: ",".join(f"{a:.3f}" for a in v))
    d["xs"] = d["x"].map(json.dumps)
    d["ws"] = d["w"].map(json.dumps)
    rows = []
    for (prop, sid, xk, tk), g in d.groupby(["property", "system_id", "x_key", "T_key"], sort=True):
        v = g["value"].to_numpy()
        med = float(np.median(v))
        q25, q75 = np.percentile(v, [25, 75])
        conflict = bool(len(v) >= 2 and (np.ptp(v) if len(v) == 2 else q75 - q25) > 0.05 * abs(med))
        f = g.iloc[0]
        rows.append(dict(property=prop, system_id=sid, x_key=xk, T_key=tk, P_key=101.0, T_K=f["T_K"], P_kPa=np.nan,
                         kind=f["kind"], n_comp=f["n_comp"], is_hc=True, comp_inchikeys=f["comp_inchikeys"], x_mole=f["xs"],
                         w_mass=f["ws"], fuel_id=np.nan, ref=REF, composition_level="L3", comp_coverage=1.0,
                         comp_basis="volume", value=med, vmin=float(v.min()), vmax=float(v.max()), q25=float(q25),
                         q75=float(q75), n_obs=len(v), n_refs=1, n_sources=1, sources="boehm2022", n_outlier=0,
                         conflict=conflict))
    return pd.DataFrame(rows, columns=COLS)


def t_bin(T_K: float) -> str:
    return "T_lt_263K" if T_K < 263.15 else ("T_263_283K" if T_K < 283.15 else "T_ge_283K")


def main():
    log: dict = {}
    out = build(RAW, log)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False, float_format="%.10g")
    v = out[(out.property == "viscosity_dyn_mPas") & out.kind.isin(["binary", "multi"])]
    summ = {b: dict(points=int((v["T_K"].map(t_bin) == b).sum()),
                    systems=int(v[v["T_K"].map(t_bin) == b].system_id.nunique())) for b in ("T_lt_263K", "T_263_283K", "T_ge_283K")}
    log["dyn_viscosity_binary_multi_by_bin"] = summ
    (OUT.parent / "boehm2022_build_log.json").write_text(json.dumps(log, indent=1, default=str))
    print(json.dumps(log, indent=1, default=str))
    print(f"wrote {OUT} ({len(out)} rows)")


if __name__ == "__main__":
    main()
