"""v6 flash arm: Leenhouts benchmark flash-point mixtures, ALL chemistries (hydrocarbon, oxygenate, other).

Source: Leenhouts, Larsson, Verhelst, Vermeire, Fuel 381 (2025) 133218, Zenodo benchmark
``datasets/flashpoint/processed_data/flashpoint_mix_exp.csv`` (read from ``datasets.zip`` in memory, never extracted).
Native schema: ``value`` (unit string ``Celcius``, sic), ``fraction_type`` (all ``molar`` in the real file), up to four
``component_inchi_k`` / ``fraction_k`` pairs, no temperature or pressure column.

Output is the canonical ``measurements.csv`` schema (AGENTS.md) plus ``chem_class`` (``hc`` | ``oxy`` | ``other``).
Staging only: this loader never writes to ``data/unified/``. The unified tables keep ``is_hc`` rows only; the oxygenate arm
reads this staged table instead (see ``results/v6/flash_oxy/PREDECLARATION.md``).

Conventions: flash value in degC (canonical unit), P = 101.325 kPa (closed-cup, the benchmark gives none), T_K = NaN (the
flash temperature is the label, not a condition). ``system_id`` = md5[:12] of the sorted InChIKey set, imported from the
production helper. Dedup on (system_id, x rounded to 3 decimals): median, ``n_obs``, flag (never drop) ``conflict`` when
the range exceeds 2 degC. The real file has no repeated compositions, so dedup is a no-op there; it is tested on a
constructed repeat.
"""
from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

from fuelblend.data.thermoml_sle import system_id

ZIP_MEMBER = "datasets/flashpoint/processed_data/flashpoint_mix_exp.csv"
MAX_COMP = 4
FLASH_CONFLICT_C = 2.0
P_KPA = 101.325
COLS = ["source", "property", "value", "T_K", "P_kPa", "kind", "n_comp", "is_hc", "chem_class", "system_id",
        "comp_inchikeys", "x_mole", "w_mass", "fuel_id", "ref", "raw_value", "raw_unit",
        "vmin", "vmax", "n_obs", "conflict"]


def _read_raw(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as z:
            return pd.read_csv(io.BytesIO(z.read(ZIP_MEMBER)))
    return pd.read_csv(path)


def _chem(inchi: str, cache: dict) -> tuple | None:
    if inchi not in cache:
        m = Chem.MolFromInchi(inchi)
        if m is None:
            cache[inchi] = None
        else:
            syms = {a.GetSymbol() for a in m.GetAtoms()}
            cls = "hc" if syms <= {"C", "H"} else ("oxy" if syms <= {"C", "H", "O"} else "other")
            cache[inchi] = (Chem.MolToInchiKey(m), float(MolWt(m)), cls)
    return cache[inchi]


def build(path: Path, log: dict | None = None) -> pd.DataFrame:
    """Canonical flash-mixture rows (deduplicated) from the Leenhouts zip or a csv with the same columns."""
    raw = _read_raw(path)
    log = {} if log is None else log
    assert set(raw["unit"].unique()) <= {"Celcius"}, raw["unit"].unique()  # native unit string, sic
    assert set(raw["fraction_type"].unique()) <= {"molar"}, raw["fraction_type"].unique()
    cache: dict = {}
    recs, dropped = [], 0
    for idx, r in raw.iterrows():
        items = []
        for k in range(MAX_COMP):
            inchi, x = r.get(f"component_inchi_{k}"), r.get(f"fraction_{k}")
            if isinstance(inchi, str) and pd.notna(x) and float(x) > 0:
                c = _chem(inchi.strip(), cache)
                if c is None:
                    items = None
                    break
                items.append((c[0], c[1], c[2], float(x)))
        if not items or pd.isna(r["value"]):
            dropped += 1
            continue
        items.sort(key=lambda t: t[0])
        tot = sum(t[3] for t in items)
        x_mole = [t[3] / tot for t in items]
        m = [xi * t[1] for xi, t in zip(x_mole, items)]
        w_mass = [mi / sum(m) for mi in m]
        classes = {t[2] for t in items}
        cc = "hc" if classes == {"hc"} else ("oxy" if classes <= {"hc", "oxy"} else "other")
        iks = [t[0] for t in items]
        recs.append(dict(source="leenhouts_flash_mix", property="flash_point_C", value=float(r["value"]), T_K=np.nan,
                         P_kPa=P_KPA, kind="binary" if len(items) == 2 else "multi", n_comp=len(items),
                         is_hc=(cc == "hc"), chem_class=cc, system_id=system_id(iks), comp_inchikeys=json.dumps(iks),
                         x_mole=json.dumps([round(v, 6) for v in x_mole]),
                         w_mass=json.dumps([round(v, 6) for v in w_mass]), fuel_id=f"leenhouts_flash_mix_{idx}",
                         ref="Leenhouts2025:flashpoint_mix_exp", raw_value=float(r["value"]), raw_unit="Celcius",
                         _xkey=tuple(round(v, 3) for v in x_mole)))
    df = pd.DataFrame(recs)
    n_in = len(df)
    out = []
    for _, g in df.groupby(["system_id", "_xkey"], sort=False):
        row = g.iloc[0].copy()
        row["value"] = float(g["value"].median())
        row["vmin"], row["vmax"], row["n_obs"] = float(g["value"].min()), float(g["value"].max()), len(g)
        row["conflict"] = bool(row["vmax"] - row["vmin"] > FLASH_CONFLICT_C)
        out.append(row)
    res = pd.DataFrame(out).drop(columns=["_xkey"]).reset_index(drop=True)[COLS]
    log.update(rows_raw=int(len(raw)), rows_dropped=int(dropped), rows_pre_dedup=int(n_in), rows_out=int(len(res)),
               systems=int(res.system_id.nunique()),
               by_class={k: dict(rows=int(len(v)), systems=int(v.system_id.nunique())) for k, v in res.groupby("chem_class")},
               n_conflict=int(res.conflict.sum()))
    return res
