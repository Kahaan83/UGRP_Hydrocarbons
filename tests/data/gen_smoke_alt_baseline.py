"""Generates tests/data/smoke_alt_baseline.json: engine_v4.evaluate_blend results for the v5pool design recipes
(weights from results/v5/BASE/delta_design_ml_visc/designs/v5pool/designs.csv), captured BEFORE the alt-limit change."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from fuelblend.design import engine_v4 as E, library_route as LR

ROOT = Path(__file__).resolve().parents[2]
MARGINS = ROOT / "configs" / "uncertainty_v5.yaml"
CSV = ROOT / "results/v5/BASE/delta_design_ml_visc/designs/v5pool/designs.csv"


def snapshot():
    d = pd.read_csv(CSV)
    pool = json.loads(d.pool.iloc[0])
    models = LR.build_route_models(pool)
    marg = E.load_margins(MARGINS)
    out = {}
    for r in d.itertuples():
        wd = json.loads(r.weights)
        w = np.array([wd[n] for n in models["pool"]], float); w /= w.sum()
        spec = E.load_spec(r.spec_name)
        res = E.evaluate_blend(spec, models["pool"], w, models, marg)
        out[r.design_id] = {c.key: dict(status=c.status, pred=c.prediction, slack=c.slack, margin=c.applied_margin,
                                        reason=c.reason) for c in res}
    return out


if __name__ == "__main__":
    s = snapshot()
    (ROOT / "tests/data/smoke_alt_baseline.json").write_text(json.dumps(s, indent=1, default=str))
    print({k: len(v) for k, v in s.items()})
