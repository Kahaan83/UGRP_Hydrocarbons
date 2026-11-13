#!/usr/bin/env python
"""v5 delta_api_route Task 7: inverse design on a user pool (names / SMILES / InChIKeys) through the route adapter.

Models: ``fuelblend.design.library_route.build_route_models`` (every constant carries a status; a model that touches a
non-verified constant is uncertified).  Optimiser and statuses: ``engine_v4.design`` through ``design_routed``.
Margins: results/v5/S1/uncertainty_v5.yaml (as D4).  Freezing has no valid margin (8 systems -> at most an 88.9 % level),
so a design is category A (fails a limit at the design OR the API point estimate), B (meets every limit on both, cannot be
certified) or certified.  The forward check judges every constrained property against the ``predict_fuel`` value.  ``--freezing-guard-K 4`` optimises against a freezing limit tightened by 4 K and reports against the real
spec: "indicative guard band, not certified".

Pools: ``v5pool`` (the D4 pool), ``ext`` (v5pool plus the lowest-InChIKey registry compound with own grid, flash and SLE data
in each declared (subkind, nC) bin; an empty bin is recorded in bins.csv and skipped), or a comma list / .txt file of ids.

Outputs in --out: designs.csv, constraints.csv, forward_check.csv, runtime.txt (and bins.csv for ``ext``).
"""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.design import engine_v4 as E  # noqa: E402
from fuelblend.design import library_route as LR  # noqa: E402

_s = importlib.util.spec_from_file_location("run_designs_v4", ROOT / "scripts/v4/run_designs_v4.py")
R4 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(R4)

V5POOL = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
BINS = [("n_paraffin", 13), ("n_paraffin", 14), ("cyclo_mono", 9), ("cyclo_mono", 10), ("indan_tetralin", 10),
        ("aromatic_di", 11)]
D4_DESIGNS = ROOT / "results/v5/D4/designs_v4.csv"
D4_ID = {"astm:jet_a": "jet_astm_d1655__v5pool", "bis:jet_a1": "jet_bis_is1571__v5pool"}
MARGINS = ROOT / "configs/uncertainty_v5.yaml"      # S1 yaml + flash lab margin 2.2069 C (user decision Q2)


def extended_pool():
    """v5pool + one registry compound per declared bin.  Returns (pool ids, bins dataframe)."""
    from rdkit import Chem
    from rdkit.Chem import inchi as rdinchi
    from fuelblend.api import predict_route as pr
    from fuelblend.data.reconstruction_r3 import subkind
    ctx = pr.get_ctx()
    reg, P = ctx.reg, pr.predicates(ctx)
    have = {rdinchi.MolToInchiKey(Chem.MolFromSmiles(s)) for s in
            pd.read_csv(LR.PURE_PROPS_V5).drop_duplicates("compound").set_index("compound").loc[V5POOL, "smiles"]}
    cand: dict = {}
    for i in range(len(reg.ik)):
        if not reg.usable[i]:
            continue
        k = subkind(reg.smiles[i])
        if k is None:
            continue
        m = Chem.MolFromSmiles(reg.smiles[i])
        nc = sum(a.GetSymbol() == "C" for a in m.GetAtoms())
        if (k, nc) not in BINS:
            continue
        flags = dict(grid=P["HAS_GRID"](i), flash=P["HAS_FLASH"](i), sle=P["HAS_SLE"](i))
        cand.setdefault((k, nc), []).append((reg.ik[i], reg.smiles[i], flags))
    pool, rows = list(V5POOL), []
    for b in BINS:
        c = sorted(cand.get(b, []))
        ok = [x for x in c if all(x[2].values())]
        row = dict(subkind=b[0], nC=b[1], n_registry=len(c), **{f"n_{k}": sum(x[2][k] for x in c) for k in ("grid", "flash", "sle")},
                   n_qualifying=len(ok), inchikey="", smiles="", outcome="")
        if not ok:
            row["outcome"] = "empty bin, skipped (no registry compound has own grid, flash and SLE data)"
        else:
            ik, smi, _ = ok[0]
            row.update(inchikey=ik, smiles=smi)
            if ik in have or ik in {r["inchikey"] for r in rows}:
                row["outcome"] = "already in v5pool, adds nothing"
            else:
                row["outcome"] = "added"
                pool.append(ik)
        rows.append(row)
    return pool, pd.DataFrame(rows)


def resolve_pool(arg):
    if arg == "v5pool":
        return list(V5POOL), None
    if arg == "ext":
        return extended_pool()
    p = Path(arg)
    if p.suffix == ".txt" and p.exists():
        return [x.strip() for x in p.read_text(encoding="utf-8").splitlines() if x.strip() and not x.startswith("#")], None
    return [x.strip() for x in arg.split(",") if x.strip()], None


def guard_configs(spec_name, K, tmp):
    """Copy of the spec yaml with the freezing limit tightened by K (optimiser side only)."""
    std, fam = spec_name.split(":")
    doc = yaml.safe_load(open(E.CONFIGS / f"specs_{std}.yaml", encoding="utf-8"))
    c = doc["families"][fam]["constraints"]["freezing_point_C"]
    c["max"] = float(c["max"]) - K
    out = Path(tmp)
    for f in E.CONFIGS.glob("specs_*.yaml"):
        if f.name != f"specs_{std}.yaml":
            (out / f.name).write_text(f.read_text(encoding="utf-8"), encoding="utf-8")
    yaml.safe_dump(doc, open(out / f"specs_{std}.yaml", "w", encoding="utf-8"), sort_keys=False)
    return out


def blocking_gaps(res):
    return {c.key: f"{c.status}: {c.reason}" for c in res.constraints if c.status not in E.OK_STATUSES}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--spec", default="astm:jet_a,bis:jet_a1", help="astm:jet_a|bis:jet_a1 (comma list)")
    ap.add_argument("--pool", default="v5pool", help="v5pool | ext | comma ids | .txt file")
    ap.add_argument("--n-starts", type=int, default=None, help="default 30, or 8 above 20 components")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--freezing-guard-K", type=float, default=0.0, choices=[0.0, 4.0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--synthetic-blend", action="store_true",
                    help="add the D7566 aromatics minimum (8 vol %%); a pool with no aromatic member is blocked")
    ap.add_argument("--alt-limits", action="store_true", help="honour the D1655 alternative smoke limit (use_alt_limits)")
    ap.add_argument("--recipe-margins", action="store_true",
                    help="viscosity and smoke margins depend on the recipe (margin_fn) and are seen by the optimiser (default off)")
    a = ap.parse_args(argv)
    out = Path(a.out) if Path(a.out).is_absolute() else ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)
    pool_ids, bins = resolve_pool(a.pool)
    if bins is not None:
        bins.to_csv(out / "bins.csv", index=False)
        print(bins[["subkind", "nC", "n_qualifying", "outcome", "inchikey"]].to_string(index=False), flush=True)
    t0 = time.time()
    models = LR.build_route_models(pool_ids)
    pool = models["pool"]
    n_starts = a.n_starts or LR.default_n_starts(pool)
    t_models = time.time() - t0
    print(f"pool ({len(pool)}): {pool}\nn_starts={n_starts} seed={a.seed} guard={a.freezing_guard_K}K models {t_models:.1f}s", flush=True)
    d_rows, c_rows, f_rows, rt = [], [], [], [f"pool_size={len(pool)} n_starts={n_starts} seed={a.seed} guard_K={a.freezing_guard_K} "
                                              f"build_models_s={t_models:.1f}"]
    d4 = pd.read_csv(D4_DESIGNS).set_index("design_id") if D4_DESIGNS.exists() else None
    label = a.label or ("indicative guard band, not certified" if a.freezing_guard_K else "")
    for spec_name in [s.strip() for s in a.spec.split(",")]:
        did = f"{spec_name.replace(':', '_')}__{a.pool if len(a.pool) < 12 else 'custom'}"
        t1 = time.time()
        if a.freezing_guard_K:
            with tempfile.TemporaryDirectory() as tmp:
                cdir = guard_configs(spec_name, a.freezing_guard_K, tmp)
                use = LR.with_recipe_margins(models, E.load_margins(MARGINS)) if a.recipe_margins else models
                opt = LR.annotate(E.design(spec_name, pool, use, n_starts=n_starts, seed=a.seed, margins_path=MARGINS,
                                           configs_dir=cdir,
                                           synthetic_blend=a.synthetic_blend, alt_limits=a.alt_limits),
                                  models=models, alt_limits=a.alt_limits, recipe_margins=a.recipe_margins,
                                  margins_path=MARGINS if a.recipe_margins else None)
            res = LR.evaluate_recipe(spec_name, models, opt.weights, margins_path=MARGINS,
                                    synthetic_blend=a.synthetic_blend, alt_limits=a.alt_limits,
                                    recipe_margins=a.recipe_margins)   # judged on the real spec
            res.optimiser = opt.optimiser
        else:
            res = LR.design_routed(spec_name, models, n_starts=n_starts, seed=a.seed, margins_path=MARGINS,
                                   synthetic_blend=a.synthetic_blend, alt_limits=a.alt_limits,
                                   recipe_margins=a.recipe_margins)
        dt = time.time() - t1
        if res.category:                                  # blocked: no optimisation, no forward check
            d_rows.append(dict(design_id=did, spec_name=spec_name, pool_arg=a.pool, pool_size=len(pool), pool=json.dumps(pool),
                               n_starts=0, seed=a.seed, category=res.category, passed=False,
                               blocking=json.dumps({"aromatics_vol_pct": res.binding_constraints[0]["reason"]})))
            print(f"[{did}] {res.category}: {res.binding_constraints[0]['reason']}", flush=True)
            continue
        viol, nopred = R4.point_estimate_view(res.constraints)
        cat = R4.category(res, viol, nopred)
        t2 = time.time()
        fc = LR.design_forward_check(res)
        fc_t = time.time() - t2
        judged = fc[fc["tol"].notna()]
        fwd_pass = bool(judged["pass"].all()) and len(judged) > 0          # like-for-like tolerance rows only
        lim = fc[fc["limit"] != ""]
        both = LR.limits_met_by_both(fc)
        api_failed = lim.loc[~lim["api_meets"].astype(bool), "constraint"].tolist()
        design_failed = lim.loc[~lim["design_meets"].astype(bool), "constraint"].tolist()
        if res.passed and both:
            cat_final = "C_certified"
        elif not both:
            cat_final = "A_fails_a_limit_on_design_or_API_point_estimate"
        else:
            cat_final = "B_meets_both_not_certifiable"
        gaps = blocking_gaps(res)
        w = {k: round(v, 4) for k, v in res.weights.items()}
        d_rows.append(dict(design_id=did, spec_name=spec_name, pool_arg=a.pool, pool_size=len(pool), pool=json.dumps(pool),
                           n_starts=n_starts, seed=a.seed, freezing_guard_K=a.freezing_guard_K, label=label,
                           weights=json.dumps(w), passed=res.passed, meets_at_point_estimate=bool(viol == 0 and nopred == 0),
                           n_point_violation=viol, n_no_point_estimate=nopred, category_design_side=cat, category=cat_final,
                           meets_limits_design_and_api=both, design_failed_limits=json.dumps(design_failed),
                           api_failed_limits=json.dumps(api_failed),
                           best_min_slack=res.optimiser["best_min_slack"], n_converged=res.optimiser["n_converged"],
                           blocking=json.dumps(gaps), forward_check_pass=fwd_pass,
                           forward_failed=json.dumps(judged.loc[~judged["pass"].astype(bool), "constraint"].tolist()),
                           runtime_design_s=round(dt, 2), runtime_forward_s=round(fc_t, 2)))
        for c in res.constraints:
            c_rows.append(dict(design_id=did, constraint=c.key, model=c.model, lo=c.lo, hi=c.hi, prediction=c.prediction,
                               applied_margin=c.applied_margin, status=c.status, slack=c.slack, reason=c.reason,
                               spec_verification=c.spec_verification))
        f_rows.append(fc.assign(design_id=did))
        rt.append(f"{did}: design {dt:.1f}s ({res.optimiser['n_converged']}/{n_starts} converged), forward check {fc_t:.1f}s")
        print(f"[{did}] {cat_final} (design side {cat}) passed={res.passed} fwd_pass={fwd_pass} api_failed={api_failed} {dt:.1f}s\n  weights {w}", flush=True)
        for k, v in gaps.items():
            print(f"  blocked: {k}: {v[:160]}")
        if d4 is not None and a.pool == "v5pool" and not a.freezing_guard_K:
            ref = json.loads(d4.loc[D4_ID[spec_name], "weights"]) if spec_name in D4_ID else {}
            if ref:
                diff = {k: round(w.get(k, 0.0) - ref.get(k, 0.0), 4) for k in sorted(set(ref) | set(w))}
                print(f"  difference to D4 recipe (route minus D4): {diff}  max |d| = {max(abs(x) for x in diff.values()):.4f}")
                d_rows[-1]["max_abs_diff_to_D4"] = max(abs(x) for x in diff.values())
                d_rows[-1]["diff_to_D4"] = json.dumps(diff)
    pd.DataFrame(d_rows).to_csv(out / "designs.csv", index=False)
    pd.DataFrame(c_rows).to_csv(out / "constraints.csv", index=False)
    pd.concat(f_rows).to_csv(out / "forward_check.csv", index=False) if f_rows else None   # all blocked: no rows
    (out / "runtime.txt").write_text("\n".join(rt) + f"\ntotal_s={time.time() - t0:.1f}\n", encoding="utf-8")
    print("wrote", out)


if __name__ == "__main__":
    main()
