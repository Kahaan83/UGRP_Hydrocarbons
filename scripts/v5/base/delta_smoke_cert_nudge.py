"""delta_smoke_cert review 3 diagnostic: nudge every iso-paraffin pool member of a returned optimum to 1 wt % (others renormalised)
and re-judge with the design flags.  Declared in results/v5/BASE/delta_smoke_cert/README.md before it was run.
usage: python delta_smoke_cert_nudge.py <tag>   (writes results/v5/BASE/delta_smoke_cert/nudge_<tag>.csv)"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.design import library_route as LR  # noqa: E402

MARGINS = ROOT / "configs/uncertainty_v5.yaml"
D = ROOT / "results/v5/BASE/delta_smoke_cert/designs"
RUNS = ["v5pool", "guard4_v5pool", "ext_n30", "ext_n8", "poolA_toluene_isooctane", "poolB_hexylbenzene_isooctane"]
TARGET = 0.01


def judge(spec, models, w):
    r = LR.evaluate_recipe(spec, models, w, margins_path=MARGINS, synthetic_blend=True, alt_limits=True, recipe_margins=True)
    return {c.key: c for c in r.constraints}


def main(tag):
    rows = []
    for run in RUNS:
        dd = pd.read_csv(D / run / "designs.csv")
        for _, d in dd.iterrows():
            pool = json.loads(d["pool"])
            w0 = json.loads(d["weights"])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                models = LR.build_route_models(pool)
                iso = [n for n in models["pool"] if models["classes"].get(n) == "iso_paraffin"]
                base = judge(d["spec_name"], models, w0)
                if not iso:
                    rows.append(dict(run=run, design=d["design_id"], nudged="none: no iso-paraffin member in the pool",
                                     smoke_status_before=base["smoke_point_mm"].status, note="class-mix gate cannot be met in this pool"))
                    continue
                w1 = {n: v for n, v in w0.items()}
                for n in iso:
                    w1[n] = max(w1.get(n, 0.0), TARGET)
                tot = sum(w1.values())
                w1 = {n: v / tot for n, v in w1.items()}
                nud = judge(d["spec_name"], models, w1)
            s0, s1 = base["smoke_point_mm"], nud["smoke_point_mm"]
            row = dict(run=run, design=d["design_id"], nudged=",".join(iso), iso_wt_before=sum(w0.get(n, 0) for n in iso),
                       iso_wt_after=sum(w1[n] for n in iso),
                       smoke_pred_before=s0.prediction, smoke_pred_after=s1.prediction,
                       smoke_status_before=s0.status, smoke_status_after=s1.status,
                       smoke_margin_before=s0.applied_margin, smoke_margin_after=s1.applied_margin,
                       smoke_reason_after=s1.reason)
            for k, c in nud.items():
                if k != "smoke_point_mm":
                    row[f"dslack_{k}"] = c.slack - base[k].slack
                    row[f"status_{k}_after"] = c.status
            rows.append(row)
            print(run, d["design_id"], s0.status, "->", s1.status, "pred", round(s0.prediction, 2), "->", round(s1.prediction, 2),
                  "margin", s1.applied_margin, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / f"results/v5/BASE/delta_smoke_cert/nudge_{tag}.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
