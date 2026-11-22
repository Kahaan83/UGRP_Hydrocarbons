"""Exploratory blend properties from pure-compound predictions (linear mole-fraction rule).

Evidence: results/v7/final/blend_from_pure (PREDECLARATION.md, READING.md). Added to `predict_fuel` as extra keys only; no
existing key, value or margin changes. Every entry is uncertified (no conformal margin, no lo/hi, not in_domain).
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SUMMARY = ROOT / "results/v7/final/blend_from_pure/summary.csv"
# output key -> (pure property, task whose evidence applies, unit); None task = tested and rejected
SPEC = {
    "dcn_blend_exploratory": ("dcn", "dcn_leenhouts", "-"),
    "ron_blend_exploratory": ("ron", "ron_kaust", "-"),
    "mon_blend_exploratory": ("mon", None, "-"),
}
_NO_SKILL = {"mon_blend_exploratory": "MON: usable on KAUST surrogates (5 systems) but no skill on Leenhouts mixtures (32 systems); "
                                      "evidence disagrees across sources, status: no skill shown"}


@lru_cache(maxsize=1)
def _evidence() -> dict:
    import pandas as pd
    return pd.read_csv(SUMMARY).set_index("task").to_dict("index") if SUMMARY.exists() else {}


def _entry(key: str, status: str, reason: str, value=None, task=None, extra=None) -> dict:
    ev = _evidence().get(task) if task else None
    d = {"value": value, "std": None, "unit": SPEC[key][2], "method": "linear_mole_fraction_of_best_pure_models" if value is not None else "none",
         "status": status, "scope": "exploratory", "uncertified": True, "lo": math.nan, "hi": math.nan, "margin": None, "in_domain": False,
         "warnings": [reason, "exploratory: dev-row evidence only, no conformal margin, no lab validation"]}
    if ev:
        d["evidence"] = {"task": task, "n_systems": int(ev["n_systems"]), "n_rows": int(ev["n_rows"]), "mae": float(ev["mae_model"]),
                         "mae_constant": float(ev["mae_const"]), "source": "results/v7/final/blend_from_pure/summary.csv"}
    if extra:
        d.update(extra)
    return d


def exploratory_props(comps, ctx, composition_kind: str = "known") -> dict:
    """{key: entry} for the exploratory blend properties. Never raises: a missing model gives status 'model_missing'."""
    from fuelblend.api import predict_pure as PP
    out = {}
    mw = np.array([ctx.mw(c.smiles) for c in comps])
    x = np.array([c.w for c in comps]) / mw
    x = x / x.sum()
    for key, (prop, task, _) in SPEC.items():
        if task is None:
            out[key] = _entry(key, "no_skill_shown", _NO_SKILL[key])
            continue
        res = PP.predict_pure([c.smiles for c in comps], props=[prop])
        vals = [res[c.smiles][prop] for c in comps]
        if any(v["status"] != "ok" for v in vals):
            why = next(v.get("reason") or v["status"] for v in vals if v["status"] != "ok")
            out[key] = _entry(key, "model_missing", f"pure model unavailable: {why}", task=task)
            continue
        e = _entry(key, "exploratory", "linear mole-fraction rule over best pure-compound predictions; tested on Leenhouts "
                   "and KAUST mixtures only" if task else "", value=float(np.dot(x, [v["value"] for v in vals])), task=task,
                   extra={"components": [{"inchikey": c.ik, "x_mole": float(xi), "pure_value": v["value"]}
                                         for c, xi, v in zip(comps, x, vals)]})
        if composition_kind != "known":
            e["warnings"].append(f"composition_kind={composition_kind}: evidence is for known mixtures")
        if any(not c.is_hc for c in comps):
            e["warnings"].append("non_hc_component: mixtures tested were largely hydrocarbon")
        out[key] = e
    return out
