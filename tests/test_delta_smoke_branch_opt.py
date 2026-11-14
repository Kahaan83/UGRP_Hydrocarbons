"""delta_smoke_branch_opt R-SB-1: the optimiser smoke stand-in (``_Sides(alt_point_standin=True)``) is off by default, changes nothing on a pool
whose alternative requirements all have margins, and on the recorded WS7 P2 ASTM recipe makes the optimiser's smoke slack equal the judge's."""
import importlib.util
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.design import blendstock_v5 as B
from fuelblend.design import engine_v4 as E
from fuelblend.design import library_route as LR

ROOT = Path(__file__).resolve().parents[1]
warnings.simplefilter("ignore")
_s = importlib.util.spec_from_file_location("ws7run_sb", ROOT / "scripts/v5/base/delta_ws7_design_realism.py")
W = importlib.util.module_from_spec(_s)
_s.loader.exec_module(W)
SPEC = "astm:jet_a"


def _labels(sd):
    return [f"{it[0]}_{nm}" for it in sd.items for nm, b in (("min", it[1].get("min")), ("max", it[1].get("max"))) if b is not None]


def _sides(models, standin=None):
    sp = E.load_spec(SPEC, synthetic_blend=True, alt_limits=True, d86=True)
    mg = E.load_margins(W.MARGINS)
    use = LR.with_recipe_margins(models, mg)
    kw = {} if standin is None else dict(alt_point_standin=standin)
    return E._Sides(sp, list(models["pool"]), use, mg, graded=True, **kw)


@pytest.fixture(scope="module")
def p2():
    bp = W.make_pool("P2")
    models = B.add_d86_models(B.build_models(bp))
    d = pd.read_csv(ROOT / "results/v5/BASE/delta_ws7_design_realism/designs/P2/astm_jet_a__P2__designs.csv").iloc[0]
    top = json.loads(d.top_weights)
    u = np.array([top[n] for n in bp.names], float)
    u = u / u.sum()
    return bp, models, u


def test_default_is_off_and_matches_explicit_off(p2):
    bp, models, u = p2
    w = bp.expand(u)
    a, b = _sides(models), _sides(models, standin=False)
    assert a.alt_point_standin is False and a.alt.get("smoke_point_mm") is None
    assert np.array_equal(a.vector(w), b.vector(w))


def test_standin_changes_only_smoke_and_matches_judge(p2):
    bp, models, u = p2
    w = bp.expand(u)
    off, on = _sides(models), _sides(models, standin=True)
    assert on.alt.get("smoke_point_mm") is not None
    vo, vn, lab = off.vector(w), on.vector(w), _labels(off)
    changed = [k for k, a, b in zip(lab, vo, vn) if a != b]
    assert changed == ["smoke_point_mm_min"]
    r = B.judge(SPEC, models, w, W.MARGINS, True, True, True, True, bp.represented_mass(u))
    judge_smoke = {c.key: c for c in r.constraints}["smoke_point_mm"]
    assert judge_smoke.status == E.STATUS_UNCERTIFIED                        # the judge is unchanged: never certifies this branch
    assert abs(vn[lab.index("smoke_point_mm_min")] - judge_smoke.slack) < 1e-9


def test_standin_no_effect_on_pure_pool_c():
    bp = W.make_pool("C")
    models = B.build_models(bp)
    w = np.full(len(models["pool"]), 1.0 / len(models["pool"]))
    assert np.array_equal(_sides(models).vector(w), _sides(models, standin=True).vector(w))
