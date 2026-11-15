"""delta_smoke_branch_opt (declared in results/v5/BASE/delta_smoke_branch_opt/README.md before any run).

R-SB-2: P2/P3 x ASTM/BIS with the optimiser smoke stand-in on (``alt_standin``), WS7 protocol otherwise unchanged.
R-SB-3: P0x-a = pool C + tetralin + n-tridecane (pure, D86 on).  R-SB-4: P0x-b = P0x-a + hexylbenzene (only after its constants are verified).

    python scripts/v5/base/delta_smoke_branch_opt.py design --pool P2 --spec astm:jet_a
    python scripts/v5/base/delta_smoke_branch_opt.py table
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from fuelblend.design import blendstock_v5 as B  # noqa: E402

_s = importlib.util.spec_from_file_location("ws7run", ROOT / "scripts" / "v5" / "base" / "delta_ws7_design_realism.py")
W = importlib.util.module_from_spec(_s)
_s.loader.exec_module(W)

OUT = ROOT / "results" / "v5" / "BASE" / "delta_smoke_branch_opt"
POOL_X = {"P0xa": W.POOL_C + ["tetralin", "n-tridecane"], "P0xb": W.POOL_C + ["tetralin", "n-tridecane", "hexylbenzene"]}


def make_pool(pool: str) -> B.BlendPool:
    if pool in POOL_X:
        return B.BlendPool([B.pure_component(n) for n in POOL_X[pool]])
    return W.make_pool(pool)


def run_one(pool: str, spec: str, n_starts: int = 10):
    t0 = time.time()
    bp = make_pool(pool)
    models = B.add_d86_models(B.build_models(bp))
    t1 = time.time()
    cap = B.SYNTHETIC_CAP if pool in ("P2", "P3") else None
    res, top, info = B.design_pool(spec, bp, models, W.MARGINS, n_starts=n_starts, seed=42, d86=True, lam=B.LAMBDA_UNC,
                                   synthetic_cap=cap, alt_standin=True)
    d, c = W._rows(res, top, pool + "_standin", spec)
    d.update(build_s=round(t1 - t0, 1), design_s=round(time.time() - t1, 1), pool=pool, n_compounds=len(bp.comp_ids), d86=True,
             alt_standin=True)
    out = OUT / "designs" / pool
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([d]).to_csv(out / f"{d['design_id']}__designs.csv", index=False)
    pd.DataFrame(c).to_csv(out / f"{d['design_id']}__constraints.csv", index=False)
    print(json.dumps(d, indent=1, default=str))
    for x in c:
        print(f"  {x['constraint']:30s} {x['status']:16s} pred {x['prediction']:.5g} slack {x['slack']:+.5f} margin {x['applied_margin']}  "
              f"{str(x['reason'])[:150]}")


def cmd_table():
    rows = []
    for f in sorted((OUT / "designs").glob("*/*__designs.csv")):
        d = pd.read_csv(f).iloc[0]
        c = pd.read_csv(str(f).replace("__designs.csv", "__constraints.csv")).set_index("constraint")
        sl = c.slack.astype(float)
        row = dict(d)
        row.update(min_slack=float(sl.min()), binding=",".join(k for k, v in sl.items() if v <= sl.min() + 0.003),
                   n_pass=int((c.status == "pass").sum()), n_fail=int((c.status == "fail").sum()),
                   n_uncertified=int((c.status == "uncertified").sum()),
                   failing=",".join(c.index[c.status == "fail"]))
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "r_designs.csv", index=False)
    pd.set_option("display.width", 250)
    print(out[["design_id", "passed", "min_slack", "binding", "n_pass", "n_fail", "n_uncertified", "failing", "top_weights"]].to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["design", "table"])
    ap.add_argument("--pool", default="P2")
    ap.add_argument("--spec", default="astm:jet_a")
    ap.add_argument("--n-starts", type=int, default=10)
    a = ap.parse_args()
    if a.cmd == "design":
        if a.pool == "P0xb":           # R-SB-4: declared to run only after hexylbenzene's verified constants are recorded
            from fuelblend.eval.composition_route import build_tmap_verified
            if "LTEQMZWBSYACLV-UHFFFAOYSA-N" not in build_tmap_verified():
                raise SystemExit("hexylbenzene is not in tm_dhfus_verified.csv; R-SB-4 cannot run")
        run_one(a.pool, a.spec, a.n_starts)
    else:
        cmd_table()


if __name__ == "__main__":
    main()
