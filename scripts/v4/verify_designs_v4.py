"""v4 Task 17: INDEPENDENT verifier of results/v4/D/designs_v4.csv.

Recomputes, with the ``thermo`` / ``chemicals`` libraries only (their own Psat, Tm, Hfus, rho, LFL correlations),
  * Psat(T) of every component and the Le Chatelier flash point (ideal solution, gamma = 1, LFL by the library or
    the Jones 0.55 * Cst rule),
  * the ideal solid-liquid-equilibrium liquidus (freezing) temperature,
  * the 15 C density by ideal volume mixing,
and reports every design/property where the engine value differs from the recomputation by MORE than the applied
margin the engine used (margin 0 / missing -> the tolerance in TOL_FALLBACK is reported, never silently skipped).

This file deliberately imports NOTHING from ``fuelblend`` (tests/test_v4_verify_imports.py enforces it).  It reads the
CSV outputs and data/processed/v4/pure_props_v4.csv (only for SMILES, to give the libraries an unambiguous identity).

Usage: python scripts/v4/verify_designs_v4.py  ->  results/v4/D/verification_v4.csv (+ printed summary)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "v4" / "D"
PURE = ROOT / "data" / "processed" / "v4" / "pure_props_v4.csv"

R_GAS = 8.314462618
P_ATM = 101325.0
TOL_FALLBACK = {"flash_point_C": 3.0, "freezing_point_C": 5.0, "density_15C_kg_m3": 10.0}   # used when the engine had no margin
SMILES_EXTRA = {"ethanol": "CCO", "isooctane": "CC(C)CC(C)(C)C", "toluene": "Cc1ccccc1"}


def _chem(name, smiles):
    from thermo import Chemical
    try:
        return Chemical(name)                        # registry names are IUPAC/common names the libraries know
    except Exception:                                # noqa: BLE001
        return Chemical(smiles)


def _lfl_vol_frac(ch, formula_counts):
    """LFL (mole fraction) from the library when it has one, else Jones' rule 0.55 * C_st."""
    try:
        from chemicals import LFL as lfl_fn          # chemicals.safety.LFL (DIPPR/IEC tables)
        v = lfl_fn(CASRN=ch.CAS)
        if v:
            return float(v), "chemicals.LFL"
    except Exception:
        pass
    c, h = formula_counts.get("C", 0), formula_counts.get("H", 0)
    o = formula_counts.get("O", 0)
    n_o2 = c + h / 4.0 - o / 2.0
    cst = 1.0 / (1.0 + 4.76 * n_o2)                  # stoichiometric mole fraction in air
    return 0.55 * cst, "Jones 0.55*Cst"


def _flash_le_chatelier(chs, x, lfls):
    def f(T):
        s = 0.0
        for ch, xi, L in zip(chs, x, lfls):
            if xi <= 0:
                continue
            s += xi * ch.VaporPressure(T) / P_ATM / L
        return s - 1.0
    try:
        return brentq(f, 200.0, 600.0, xtol=1e-4) - 273.15
    except ValueError:
        return float("nan")


def _liquidus_ideal(chs, x):
    """Highest temperature at which any component saturates in the ideal liquid (x_i = exp(-Hfus/R (1/T - 1/Tm)))."""
    best = -math.inf
    for ch, xi in zip(chs, x):
        if xi <= 0:
            continue
        if not (ch.Tm and ch.Hfus):
            return float("nan")                       # cannot check without Tm/Hfus: reported as skipped
        hf = ch.Hfus * ch.MW / 1000.0             # thermo Chemical.Hfus is J/kg -> J/mol
        T = 1.0 / (1.0 / ch.Tm - R_GAS * math.log(xi) / hf)
        best = max(best, T)
    return best - 273.15


def main(out=None):
    global OUT
    if out is not None:
        OUT = Path(out) if Path(out).is_absolute() else ROOT / out
    designs = pd.read_csv(OUT / "designs_v4.csv")
    cons = pd.read_csv(OUT / "constraints_v4.csv")
    smi = pd.read_csv(PURE).drop_duplicates("compound").set_index("compound")["smiles"].to_dict()
    smi.update({k: v for k, v in SMILES_EXTRA.items() if k not in smi})
    rows = []
    for _, d in designs.iterrows():
        w = json.loads(d.weights)
        names = [n for n, v in w.items() if v > 0]
        try:
            chs = [_chem(n, smi.get(n)) for n in names]
        except Exception as e:                        # noqa: BLE001
            rows.append(dict(design_id=d.design_id, property="all", status="skipped", note=f"library lookup failed: {e}"))
            continue
        m = np.array([c.MW for c in chs])
        wt = np.array([w[n] for n in names])
        wt = wt / wt.sum()
        x = (wt / m) / (wt / m).sum()
        c = cons[cons.design_id == d.design_id].set_index("constraint")

        def report(prop, mine, theirs_key, unit):
            if theirs_key not in c.index:
                return
            eng = c.loc[theirs_key, "prediction"]
            mg = c.loc[theirs_key, "applied_margin"]
            tol = TOL_FALLBACK[theirs_key] if not (isinstance(mg, float) and math.isfinite(mg) and mg > 0) else float(mg)
            if not (math.isfinite(mine) and isinstance(eng, float) and math.isfinite(eng)):
                rows.append(dict(design_id=d.design_id, property=prop, engine=eng, independent=mine, status="skipped",
                                 note="engine or independent value unavailable"))
                return
            diff = eng - mine
            rows.append(dict(design_id=d.design_id, property=prop, engine=eng, independent=mine, diff=diff,
                             tolerance=tol, tolerance_is_engine_margin=tol != TOL_FALLBACK[theirs_key],
                             status="DISAGREES" if abs(diff) > tol else "agrees", note=unit))

        # flash
        lfl = []
        srcs = []
        from chemicals import simple_formula_parser
        for ch in chs:
            L, s = _lfl_vol_frac(ch, simple_formula_parser(ch.formula))
            lfl.append(L)
            srcs.append(s)
        report("flash_point_C", _flash_le_chatelier(chs, x, lfl), "flash_point_C", "Le Chatelier, gamma=1, LFL: " + "/".join(sorted(set(srcs))))
        # freezing (ideal SvL)
        report("freezing_point_C", _liquidus_ideal(chs, x), "freezing_point_C", "ideal SvL liquidus from library Tm/Hfus")
        # density 15 C by ideal volume mixing
        try:
            rho = np.array([ch.MW / ch.VolumeLiquid.T_dependent_property(288.15) / 1000.0 for ch in chs])
            report("density_15C_kg_m3", float(1.0 / (wt / rho).sum()), "density_15C_kg_m3", "ideal volume mixing, thermo VolumeLiquid")
        except Exception as e:                        # noqa: BLE001
            rows.append(dict(design_id=d.design_id, property="density_15C_kg_m3", status="skipped", note=str(e)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "verification_v4.csv", index=False)
    print(out.status.value_counts().to_string())
    bad = out[out.status == "DISAGREES"]
    if len(bad):
        print("\nDisagreements larger than the margin / tolerance:")
        print(bad[["design_id", "property", "engine", "independent", "diff", "tolerance"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--out", default=None, help="folder holding designs_v4.csv / constraints_v4.csv (default results/v4/D)")
    sys.exit(main(ap.parse_args().out))
