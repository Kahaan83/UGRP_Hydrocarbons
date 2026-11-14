"""Class-mix gate for the real-fuel smoke margin (v5 delta_smoke_cert 1.4; reused by Task 1.5).

The real-fuel margin only applies to recipes whose cycloalkane, n-alkane and iso-alkane mass fractions lie
inside the observed [min, max] over the 113 mapped AJFTD GCxGC fuels (fixed before any design is scored).
Fractions per fuel = class mass_pct / sum of the fuel's mass_pct (normalised to 1)."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

REASON = "class_mix_outside_real_fuel_range"
_ROOT = Path(__file__).resolve().parents[3]
COMP = _ROOT / "data" / "real_fuels" / "v5" / "ajftd_gcxgc_composition.csv"
MAPPED = _ROOT / "data" / "real_fuels" / "v5" / "ajftd_gcxgc_mapped.csv"
_EPS = 1e-9
_CYCLO = ("monocycloparaffin", "dicycloparaffin", "tricycloparaffin")


def class_fractions(comp_path=COMP, mapped_path=MAPPED) -> pd.DataFrame:
    """Per-fuel mass fractions (cyclo, nparaffin, isoparaffin) for the mapped GCxGC fuels."""
    comp = pd.read_csv(comp_path)
    ids = set(pd.read_csv(mapped_path, usecols=["fuel_id"])["fuel_id"])
    comp = comp[comp["fuel_id"].isin(ids)]
    tot = comp.groupby("fuel_id")["mass_pct"].sum()
    piv = comp.pivot_table(index="fuel_id", columns="hc_class", values="mass_pct", aggfunc="sum", fill_value=0.0)
    out = pd.DataFrame({
        "cyclo": piv[[c for c in _CYCLO if c in piv]].sum(axis=1),
        "nparaffin": piv.get("n_paraffin", 0.0),
        "isoparaffin": piv.get("iso_paraffin", 0.0)})
    return out.div(tot.reindex(out.index), axis=0)


@lru_cache(maxsize=1)
def real_fuel_class_ranges() -> dict:
    f = class_fractions()
    return {k: (float(f[k].min()), float(f[k].max())) for k in f.columns}


def class_mix_in_real_fuel_range(cyclo_mass_frac: float, nparaffin_mass_frac: float,
                                 isoparaffin_mass_frac: float) -> tuple:
    """(True, None) if all three mass fractions are inside the real-fuel range, else (False, REASON)."""
    r = real_fuel_class_ranges()
    vals = {"cyclo": cyclo_mass_frac, "nparaffin": nparaffin_mass_frac, "isoparaffin": isoparaffin_mass_frac}
    ok = all(r[k][0] - _EPS <= v <= r[k][1] + _EPS for k, v in vals.items())    # _EPS: float round-off on the range edge
    return (True, None) if ok else (False, REASON)


# ---------------------------------------------------------------------------------------------- recipe side (Task 1.5)
# Library class names -> the GCxGC basis used above.  cycloalkane = mono+di+tri-cycloparaffin (decalins, butylcyclohexane,
# bicyclohexyl are all `cycloalkane` in the library); every other class (aromatic, oxygenate, alkene, unknown) only counts in the
# mass total, exactly as alkylbenzene/cycloaromatic/diaromatic/diene/alkene do in the per-fuel normalisation.
RECIPE_CLASS = {"cycloalkane": "cyclo", "n_paraffin": "nparaffin", "iso_paraffin": "isoparaffin"}


def recipe_class_fractions(pool, w, classes: dict) -> tuple:
    """(cyclo, nparaffin, isoparaffin) mass fractions of a recipe; ``w`` is normalised to a total of 1 (all members count)."""
    import numpy as np
    w = np.asarray(w, float)
    tot = float(w.sum())
    out = {"cyclo": 0.0, "nparaffin": 0.0, "isoparaffin": 0.0}
    if tot > 0:
        for n, wi in zip(pool, w):
            k = RECIPE_CLASS.get(classes.get(n, "unknown"))
            if k is not None:
                out[k] += float(wi) / tot
    return out["cyclo"], out["nparaffin"], out["isoparaffin"]


def class_mix_ok_for_recipe(pool, w, classes: dict) -> tuple:
    """(True, None) | (False, REASON): the recipe's class mix lies inside the real-fuel GCxGC range."""
    return class_mix_in_real_fuel_range(*recipe_class_fractions(pool, w, classes))
