#!/usr/bin/env python
"""v5: jet-only design re-run on the v5 pure-constant overlay with a pool that can physically reach the spec.

Same engine, specs, optimiser settings and margin rules as scripts/v4/run_designs_v4.py; three things differ:
  * constants: data/processed/v5/pure_props_v5.csv (scripts/v5/audit_pure_constants_v5.py)
  * margins:   results/v5/S1/uncertainty_v5.yaml (as results/v5/D3)
  * pool `v5pool`: n-decane, n-undecane, n-dodecane, butylcyclohexane, trans-decalin, 1,2,4-trimethylbenzene.
    trans-decalin supplies density without aromatics; the v4pool had no such component and could not meet density,
    freezing and smoke point together (results/v5/D3).  The pool was chosen after seeing the D3 failure, so it is
    a post-hoc choice, and the v4pool rows are re-run here unchanged for comparison.
Diesel and gasoline are not run (no validated diesel data; no multi-component octane model).
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("run_designs_v4", ROOT / "scripts/v4/run_designs_v4.py")
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)

V4POOL = ["n-decane", "n-dodecane", "n-tetradecane", "butylcyclohexane", "1,2,4-trimethylbenzene"]
V5POOL = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene"]
R.DESIGNS = [("jet_astm_d1655", "astm:jet_a", {"v4pool": V4POOL, "v5pool": V5POOL}),
             ("jet_bis_is1571", "bis:jet_a1", {"v4pool": V4POOL, "v5pool": V5POOL})]

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--out", default="results/v5/D4")
    ap.add_argument("--margins", default="results/v5/S1/uncertainty_v5.yaml")
    ap.add_argument("--pure-props", default="data/processed/v5/pure_props_v5.csv")
    a = ap.parse_args()
    R.main(a.out, a.margins, a.pure_props)
