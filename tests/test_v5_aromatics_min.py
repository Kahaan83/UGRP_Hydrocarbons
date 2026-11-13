"""v5 delta_smoke_cert 1.3: aromatics minimum (8 vol %) as a declared design option (synthetic_blend), --alt-limits flag."""
import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from fuelblend.baselines.sooting import OutOfDomainWarning, smoke_point_from_bulk
from fuelblend.design import engine_v4 as E, library_route as LR

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "tests" / "data" / "smoke_alt_baseline.json"
MARGINS = ROOT / "configs" / "uncertainty_v5.yaml"


def _aro_warnings(a):
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        smoke_point_from_bulk(a)
    return [w for w in rec if issubclass(w.category, OutOfDomainWarning) and "aromatics" in str(w.message)]


def test_aromatics_edge_no_warning():
    assert not _aro_warnings(8.0)
    assert not _aro_warnings(8.0 - 1e-9)
    assert not _aro_warnings(23.0 + 1e-9)
    assert _aro_warnings(7.9)


def test_pool_without_aromatics_blocked():
    models = LR.build_route_models(["n-decane", "n-undecane", "n-dodecane"])
    res = LR.design_routed("astm:jet_a", models, n_starts=2, margins_path=MARGINS, synthetic_blend=True)
    assert res.category == "blocked_aromatics_min_infeasible"
    assert res.passed is False
    assert any(b["reason"] == "blocked: aromatics_min infeasible for pool" for b in res.binding_constraints)
    assert res.optimiser.get("n_starts") == 0            # no optimisation run


def test_synthetic_blend_adds_constraint_only_when_set():
    for sp in ("astm:jet_a", "bis:jet_a1"):
        c = E.load_spec(sp, synthetic_blend=True)["constraints"]["aromatics_vol_pct"]
        assert c["min"] == 8.0 and c["max"] == 25.0 and c["verification"] == "secondary_sources_only"
    assert "aromatics_vol_pct" not in E.load_spec("astm:jet_a")["constraints"]


def test_bis_spec_unchanged_when_flag_off():
    import subprocess
    old = subprocess.run(["git", "show", "HEAD:configs/specs_bis.yaml"], capture_output=True, text=True, cwd=ROOT).stdout
    import yaml
    for fam, d in yaml.safe_load(old)["families"].items():
        assert E.load_spec(f"bis:{fam}")["constraints"] == d["constraints"], fam
    old_a = yaml.safe_load(subprocess.run(["git", "show", "HEAD:configs/specs_astm.yaml"], capture_output=True, text=True, cwd=ROOT).stdout)
    for fam, d in old_a["families"].items():
        assert E.load_spec(f"astm:{fam}")["constraints"] == d["constraints"], fam


def test_aromatics_min_enforced_with_synthetic_blend():
    pool = ["n-decane", "n-dodecane", "1,2,4-trimethylbenzene"]
    models = LR.build_route_models(pool)
    spec = E.load_spec("astm:jet_a", synthetic_blend=True)
    m = E.load_margins(MARGINS)
    wd = {"n-decane": 0.5, "n-dodecane": 0.49, "1,2,4-trimethylbenzene": 0.01}
    lo = E.evaluate_blend(spec, models["pool"], np.array([wd[n] for n in models["pool"]]), models, m)
    ar = [c for c in lo if c.key == "aromatics_vol_pct"][0]
    assert ar.status == E.STATUS_FAIL


def test_alt_limits_flag_sets_spec_without_mutating_caller():
    spec = E.load_spec("astm:jet_a")
    assert "use_alt_limits" not in spec or spec["use_alt_limits"] is False
    s2 = E.load_spec("astm:jet_a", alt_limits=True)
    assert s2["use_alt_limits"] is True and not spec.get("use_alt_limits", False)


def test_flag_off_identical():
    """R-SC-4: synthetic_blend=False and alt_limits=False give exactly the results frozen before the change."""
    import importlib.util
    s = importlib.util.spec_from_file_location("gen_base", ROOT / "tests" / "data" / "gen_smoke_alt_baseline.py")
    g = importlib.util.module_from_spec(s)
    s.loader.exec_module(g)
    now = json.loads(json.dumps(g.snapshot(), default=str))
    old = json.loads(BASE.read_text())
    assert now.keys() == old.keys()
    for did in old:
        assert now[did].keys() == old[did].keys(), did
        for k in old[did]:
            a, b = now[did][k], old[did][k]
            assert a["status"] == b["status"] and a["reason"] == b["reason"], (did, k)
            for f in ("pred", "slack", "margin"):
                x, y = a[f], b[f]
                assert (x is None and y is None) or (str(x) == "nan" and str(y) == "nan") or x == y, (did, k, f)
