"""Generates tests/data/margin_fn_baseline.json: pre-Task-1.5 results (evaluate_recipe on the 8 stored delta_design_ml_visc
recipes, alt off and on; robustness table; a stub design), captured BEFORE the margin_fn change.
Run: PYTHONPATH=src python tests/data/gen_margin_fn_baseline.py"""
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
import test_v5_margin_fn as T  # noqa: E402
from fuelblend.design import engine_v4 as E, library_route as LR  # noqa: E402


def main():
    out = {"evaluate_recipe": {}, "robustness": [], "stub": {}}
    for case in T._cases():
        folder, spec, pool, weights = case.values
        models = T._models_for(pool)
        for alt in (False, True):
            res = LR.evaluate_recipe(spec, models, weights, margins_path=T.MARGINS, alt_limits=alt)
            out["evaluate_recipe"][f"{case.id}|alt={alt}"] = {
                c.key: dict(status=c.status, pred=c.prediction, slack=c.slack, margin=c.applied_margin, reason=c.reason,
                            model=c.model, margin_keys=list(c.margin_keys)) for c in res.constraints}
    pool = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
    models = T._models_for(pool)
    rows = E.robustness_table("astm:jet_a", pool, [0.15, 0.13, 0.03, 0.69, 0.0, 0.0], models, margins_path=T.MARGINS)
    out["robustness"] = [[r["constraint"], r["model"], r["status"], r["prediction"], r["applied_margin"]] for r in rows]
    with tempfile.TemporaryDirectory() as tmp:
        m = T._stub_models(margin_keys=["k"])
        cdir, mp = T._stub_cfg(tmp, arom_min=10.0, with_alt=False)
        mp.write_text(yaml.safe_dump({"margins": {"k": {"applied_margin": 2.0, "unit": "mm", "direction": "over"}}}))
        r = E.design("test:jet_stub", T.POOL3, m, n_starts=5, seed=42, margins_path=mp, configs_dir=cdir)
        out["stub"] = dict(weights=[r.weights[n] for n in T.POOL3], best=r.optimiser["best_min_slack"])
    (ROOT / "tests/data/margin_fn_baseline.json").write_text(json.dumps(out, indent=1, default=str))
    print({k: len(v) for k, v in out.items()})


if __name__ == "__main__":
    main()
