"""v4 Task 17: re-run the jet / diesel / gasoline (ASTM + BIS) designs on the v4 engine.

Specs: configs/specs_astm.yaml, configs/specs_bis.yaml (only source).  Margins: configs/uncertainty_v4.yaml via
applied_margin (only source).  Physics: src/fuelblend/design/models_v4.py (Task 4-15 library code).

Each (standard family, pool) is designed twice-over by name: ``v3pool`` = the compounds the v3 designs used (so the
comparison is like-for-like) and ``v4pool`` = a pool chosen so more inputs are verified (R16).  The v3 recipe is also
RE-SCORED with the v4 engine (same weights, v4 models and margins) to separate "the recipe was never good" from "the
optimiser changed".

Outputs (results/v4/D/): designs_v4.csv, constraints_v4.csv, robustness_v4.csv, comparison_v3_v4.csv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.design import engine_v4 as E  # noqa: E402
from fuelblend.design import models_v4 as M  # noqa: E402

OUT = ROOT / "results" / "v4" / "D"
V3 = ROOT / "results" / "v3" / "U6" / "designs_v3.csv"
V3_ALIAS = {"2,2,4-trimethylpentane": "isooctane"}      # same molecule, the name that carries the v4 constants

# (v3 design id, spec name, pool id -> compounds)
DESIGNS = [
    ("jet_astm_d1655", "astm:jet_a", {
        "v3pool": ["n-decane", "n-dodecane", "decalin_trans", "tetralin"],
        "v4pool": ["n-decane", "n-dodecane", "n-tetradecane", "butylcyclohexane", "1,2,4-trimethylbenzene"]}),
    ("jet_bis_is1571", "bis:jet_a1", {
        "v3pool": ["n-decane", "n-dodecane", "decalin_trans", "tetralin"],
        "v4pool": ["n-decane", "n-dodecane", "n-tetradecane", "butylcyclohexane", "1,2,4-trimethylbenzene"]}),
    ("diesel_astm_d975", "astm:diesel_no2", {
        "v3pool": ["n-hexadecane", "n-tetradecane", "n-dodecane", "decalin_cis", "tetralin"],
        "v4pool": ["n-hexadecane", "n-tridecane", "mesitylene", "1,2,4-trimethylbenzene"]}),
    ("diesel_bis_is1460", "bis:diesel_bs6", {
        "v3pool": ["n-hexadecane", "n-tetradecane", "n-dodecane", "decalin_cis", "tetralin"],
        "v4pool": ["n-hexadecane", "n-tridecane", "mesitylene", "1,2,4-trimethylbenzene"]}),
    ("gasoline_astm_d4814", "astm:gasoline", {
        "v3pool": ["isooctane", "toluene", "isopentane", "n-heptane", "1,2,4-trimethylbenzene"],
        "v4pool": ["isooctane", "toluene", "n-heptane", "n-pentane", "cyclohexane"],
        "v4pool_dvpe": ["n-pentane", "n-heptane", "cyclohexane"]}),
    ("gasoline_bis_is2796_e10", "bis:gasoline_e10", {
        "v3pool": ["isooctane", "toluene", "isopentane", "n-heptane", "ethanol"],
        "v4pool": ["isooctane", "toluene", "n-heptane", "n-pentane", "ethanol"]}),
]


def v3_weights(design_id: str):
    d = pd.read_csv(V3).set_index("design_id").loc[design_id]
    return {V3_ALIAS.get(k, k): v for k, v in json.loads(d.components_json).items()}, d


def summarise(cons):
    st = [c.status for c in cons]
    return {s: st.count(s) for s in (E.STATUS_PASS, E.STATUS_CONSTRUCTION, E.STATUS_SCOPE, E.STATUS_FAIL, E.STATUS_UNCERTIFIED)}


def point_estimate_view(cons):
    """(n constraints whose point estimate violates a limit, n with no point estimate).  Margins ignored."""
    viol = nopred = 0
    for c in cons:
        p = c.prediction
        if p is None or not np.isfinite(p):
            nopred += 1
            continue
        if (c.lo is not None and p < c.lo) or (c.hi is not None and p > c.hi):
            viol += 1
    return viol, nopred


def category(res, viol, nopred):
    if res.passed:
        return "C_certified"
    if viol:
        return "A_fails_at_point_estimate"
    if nopred:
        return "D_undetermined_some_properties_have_no_point_estimate"
    return "B_meets_at_point_estimate_not_certifiable"


def main(out=None):
    global OUT
    if out is not None:
        OUT = Path(out) if Path(out).is_absolute() else ROOT / out
    OUT.mkdir(parents=True, exist_ok=True)
    L = M.library()
    models = M.build_models(L)
    designs, long_rows, rob_rows, cmp_rows = [], [], [], []
    for v3id, spec_name, pools in DESIGNS:
        v3w, v3row = v3_weights(v3id)
        spec = E.load_spec(spec_name)
        margins = E.load_margins()
        for pool_id, pool in pools.items():
            did = f"{v3id}__{pool_id}"
            print(f"[{did}] designing {pool} ...", flush=True)
            res = E.design(spec_name, pool, models, n_starts=30, seed=42)
            cnt = summarise(res.constraints)
            viol, nopred = point_estimate_view(res.constraints)
            designs.append(dict(
                design_id=did, v3_design_id=v3id, spec_name=spec_name, standard=spec.get("standard"),
                grade=spec.get("grade"), pool_id=pool_id, pool=json.dumps(pool),
                weights=json.dumps({k: round(v, 4) for k, v in res.weights.items()}), passed=res.passed,
                meets_at_point_estimate=bool(viol == 0 and nopred == 0), n_point_violation=viol, n_no_point_estimate=nopred,
                category=category(res, viol, nopred), extrapolating=res.extrapolating, **{f"n_{k}": v for k, v in cnt.items()},
                best_min_slack=res.optimiser["best_min_slack"],
                n_excluded_from_optimiser=len(res.optimiser["excluded_from_optimiser"]), n_converged=res.optimiser["n_converged"],
                binding_constraints=json.dumps(res.binding_constraints, default=str),
                excluded_from_optimiser=json.dumps(res.optimiser["excluded_from_optimiser"]),
                absent_class_fraction=json.dumps(res.absent_class_fraction)))
            for c in res.constraints:
                long_rows.append(dict(design_id=did, spec_name=spec_name, pool_id=pool_id, constraint=c.key, model=c.model,
                                      lo=c.lo, hi=c.hi, prediction=c.prediction, applied_margin=c.applied_margin,
                                      margin_keys=json.dumps(list(c.margin_keys)), status=c.status, slack=c.slack,
                                      test_method=c.test_method, spec_verification=c.spec_verification, reason=c.reason))
            for r in E.robustness_table(spec_name, pool, [res.weights[n] for n in pool], models):
                rob_rows.append(dict(design_id=did, **r))
            # --- v3 recipe re-scored by the v4 engine (same weights, v4 models and margins)
            if pool_id == "v3pool":
                w3 = np.array([v3w[n] for n in pool])
                w3 = w3 / w3.sum()
                re = E.evaluate_blend(spec, pool, w3, models, margins)
                rc = summarise(re)
                bad = [c for c in re if c.status not in E.OK_STATUSES]
                cmp_rows.append(dict(
                    v3_design_id=v3id, spec_name=spec_name,
                    v3_all_constraints_passed=bool(v3row.all_constraints_passed),
                    v3_recipe=json.dumps({k: round(v, 4) for k, v in v3w.items()}),
                    v3_recipe_rescored_n_pass=rc["pass"] + rc["by_construction"], v3_recipe_rescored_n_fail=rc["fail"],
                    v3_recipe_rescored_n_uncertified=rc["uncertified"],
                    v3_recipe_rescored_passed=False if bad else True,
                    v3_recipe_rescored_blocking=json.dumps({c.key: f"{c.status}: {c.reason}" for c in bad}),
                    v3_recipe_rescored_predictions=json.dumps({c.key: round(c.prediction, 3) for c in re
                                                               if np.isfinite(c.prediction)}),
                    v4_design_id=did, v4_redesigned_recipe=json.dumps({k: round(v, 4) for k, v in res.weights.items()}),
                    v4_passed=res.passed, v4_n_pass=cnt["pass"] + cnt["by_construction"], v4_n_fail=cnt["fail"],
                    v4_n_uncertified=cnt["uncertified"],
                    v4_blocking=json.dumps({c.key: f"{c.status}: {c.reason}" for c in res.constraints
                                            if c.status not in E.OK_STATUSES})))
    # reproducible check of the ridge extrapolation claim (pure n-decane, -20 C)
    fns = M.make_fns(L)
    env_rows = []
    for pool_, w_ in ((["n-decane"], [1.0]), (["n-hexadecane"], [1.0]), (["n-decane", "1,2,4-trimethylbenzene"], [0.5, 0.5])):
        for T_ in (-20.0, 20.0):
            try:
                guarded = fns["visc_ridge"](pool_, np.array(w_), T_)
                gnote = "ok"
            except E.ModelUnavailable as e:
                guarded, gnote = float("nan"), str(e)[:140]
            env_rows.append(dict(pool=json.dumps(pool_), w=json.dumps(w_), T_C=T_,
                                 ridge_unguarded_cSt=fns["visc_ridge_unguarded"](pool_, np.array(w_), T_),
                                 ridge_guarded_cSt=guarded, ridge_guard_note=gnote,
                                 recalibrated_refutas_cSt=fns["visc_recal"](pool_, np.array(w_), T_),
                                 raw_refutas_cSt=fns["visc_raw"](pool_, np.array(w_), T_)))
    pd.DataFrame(env_rows).to_csv(OUT / "viscosity_envelope_check_v4.csv", index=False)
    pd.DataFrame(designs).to_csv(OUT / "designs_v4.csv", index=False)
    pd.DataFrame(long_rows).to_csv(OUT / "constraints_v4.csv", index=False)
    pd.DataFrame(rob_rows).to_csv(OUT / "robustness_v4.csv", index=False)
    pd.DataFrame(cmp_rows).to_csv(OUT / "comparison_v3_v4.csv", index=False)
    print(f"wrote {len(designs)} designs, {len(long_rows)} constraint rows, {len(rob_rows)} robustness rows")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--out", default=None, help="output folder (default results/v4/D)")
    main(ap.parse_args().out)
