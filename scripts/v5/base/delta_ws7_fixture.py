"""delta_ws7_design_realism R-W7-4: freeze the feature-off behaviour BEFORE the engine / spec change.

Writes tests/data/ws7_feature_off_fixture.json: every spec family as ``engine_v4.load_spec`` returns it (default flags and with
synthetic_blend / alt_limits), and the pool C stored recipes (6dc66f8 r_pool.csv weights, held in memory) judged by
``library_route.evaluate_recipe`` with the smoke-pool flags (statuses, predictions, slacks, margins).  Run once, before the change."""
import json
import math
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")

import yaml  # noqa: E402

from fuelblend.design import engine_v4 as E  # noqa: E402
from fuelblend.design import library_route as LR  # noqa: E402

POOL_C = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene", "isooctane"]
POOL_C_W = {
    "astm:jet_a": {"n-decane": 0.1597, "n-undecane": 0.1357, "n-dodecane": 0.0283, "butylcyclohexane": 0.5568, "decalin_trans": 0.0,
                   "1,2,4-trimethylbenzene": 0.0934, "isooctane": 0.0261},
    "bis:jet_a1": {"n-decane": 0.0942, "n-undecane": 0.0795, "n-dodecane": 0.0144, "butylcyclohexane": 0.6945, "decalin_trans": 0.0,
                   "1,2,4-trimethylbenzene": 0.0917, "isooctane": 0.0257},
}
MARGINS = ROOT / "configs/uncertainty_v5.yaml"
OUT = ROOT / "tests/data/ws7_feature_off_fixture.json"


def _num(x):
    if x is None:
        return None
    x = float(x)
    return None if math.isnan(x) else x


def spec_snapshot() -> dict:
    out = {}
    for std in ("astm", "bis", "diesel", "gasoline", "jet"):
        p = ROOT / "configs" / f"specs_{std}.yaml"
        if not p.exists():
            continue
        fams = yaml.safe_load(open(p, encoding="utf-8")).get("families") or {}
        for fam in fams:
            for sb in (False, True):
                for al in (False, True):
                    s = E.load_spec(f"{std}:{fam}", synthetic_blend=sb, alt_limits=al)
                    out[f"{std}:{fam}|sb={sb}|alt={al}"] = json.loads(json.dumps(s, default=str))
    return out


def pool_c_snapshot() -> dict:
    models = LR.build_route_models(POOL_C)
    out = {}
    for spec, w in POOL_C_W.items():
        r = LR.evaluate_recipe(spec, models, w, margins_path=MARGINS, synthetic_blend=True, alt_limits=True, recipe_margins=True)
        out[spec] = {c.key: dict(status=c.status, prediction=_num(c.prediction), slack=_num(c.slack), margin=_num(c.applied_margin),
                                 reason=c.reason) for c in r.constraints}
    return out


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(dict(specs=spec_snapshot(), pool_c=pool_c_snapshot()), open(OUT, "w", encoding="utf-8"), indent=1, sort_keys=True)
    print("wrote", OUT)
