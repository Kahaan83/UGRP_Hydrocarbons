"""Pre-check for decision 6 (``strict_constraints=True`` by default): which committed design pools drop a spec constraint in ``_Sides``?

Builds the route models and ``engine_v4._Sides`` exactly as ``blendstock_v5.design_pool`` does for every pool behind a committed design run
(specs ``astm:jet_a`` and ``bis:jet_a1``), with NO optimisation, and prints ``pool, spec, constraint, reason`` for every ``_Sides.excluded``
entry whose reason does not start with ``"no model"`` (those are the entries ``strict_constraints`` raises on).  Returns their number.

    python scripts/v5/base/check_sides_excluded.py

The registry-additions overlay (``FUELBLEND_FLASH_REGISTRY_ADDITIONS``) is a per-process cached switch, so the overlay pools (P0y, P0y-b) run in a
child process of this script (``--group overlay``) and the exit code of that child is its count.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SPECS = ["astm:jet_a", "bis:jet_a1"]
# (label, group, loader-pool, alt_standin): the settings each committed run used (ws7: no stand-in; smoke / bicyclohexyl runs: stand-in on)
POOLS = [("P0", "plain", "P0", False), ("P1", "plain", "P1", False), ("P2", "plain", "P2", False), ("P3", "plain", "P3", False),
         ("P0x-a", "plain", "P0xa", True), ("P0x-b", "plain", "P0xb", True),
         ("P0y", "overlay", "P0y", True), ("P0y-b", "overlay", "P0yb", True)]


def _count(group: str) -> int:
    if group == "overlay":
        os.environ["FUELBLEND_FLASH_REGISTRY_ADDITIONS"] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(ROOT / "scripts" / "v5" / "base"))
    warnings.simplefilter("ignore")
    from fuelblend.design import blendstock_v5 as B
    from fuelblend.design import engine_v4 as E
    from fuelblend.design import library_route as LR
    import delta_smoke_branch_opt as SB
    import delta_pool_bicyclohexyl as BC
    W = SB.W
    n = 0
    for label, g, key, standin in POOLS:
        if g != group:
            continue
        bp = BC_pool(BC, SB, key)
        models = B.add_d86_models(B.build_models(bp))
        margins = E.load_margins(W.MARGINS)
        use = LR.with_recipe_margins(models, margins)
        for spec_name in SPECS:
            if E.blocked_result(spec_name, list(models["pool"]), models, 42, True) is not None:
                print(f"{label}, {spec_name}, -, blocked pool (not a _Sides case)", flush=True)
                continue
            spec = E.load_spec(spec_name, synthetic_blend=True, alt_limits=True, d86=True)
            sides = E._Sides(spec, list(models["pool"]), use, margins, graded=True, alt_point_standin=bool(standin))
            bad = {k: v for k, v in sides.excluded.items() if not str(v).startswith("no model")}
            for k, v in bad.items():
                print(f"{label}, {spec_name}, {k}, {v}", flush=True)
            n += len(bad)
            print(f"# {label} {spec_name}: excluded={dict(sides.excluded)}", flush=True)
    return n


def BC_pool(BC, SB, key):
    from fuelblend.design import blendstock_v5 as B
    if key in ("P0y", "P0yb"):
        return B.BlendPool([B.pure_component(n) for n in BC.pool_members(key)])
    return SB.make_pool(key)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", choices=["all", "plain", "overlay"], default="all")
    a = ap.parse_args(argv)
    if a.group != "all":
        return _count(a.group)
    n = _count("plain")
    child = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--group", "overlay"], check=False)
    n += child.returncode
    print(f"total excluded entries not starting with 'no model': {n}")
    return n


if __name__ == "__main__":
    sys.exit(main())
