import json

import numpy as np
import pandas as pd
import pytest

from fuelblend.eval import freeze_multicomp as FM
from fuelblend.eval.composition_route import ideal_sle_C

C16 = (291.34, 53359.0)
C13 = (267.79, 36162.0)
IC8 = (165.78, 9200.0)


def test_class_mix():
    assert FM.class_mix(["CCCCCCCCCCCCCCCC", "CC(C)CC(C)(C)C"]) == "alkane_only"
    assert FM.class_mix(["CCCCCCCCCCCCCCCC", "C1CCC(CC1)C1CCCCC1"]) == "alkane_naphthene"
    assert FM.class_mix(["CCCCCC", "Cc1ccc(C)cc1", "CCCCCCCCCCCCCCCC"]) == "with_aromatic"
    assert FM.class_mix(["Cc1cccc2ccccc12", "C1CCC2CCCCC2C1", "C" * 36]) == "with_aromatic"


def test_svl_ideal_matches_designated_function():
    x = [0.0017, 0.0104, 0.9879]
    cc = [C16, C13, IC8]
    iks = ["a", "b", "c"]
    tmap = {k: {"tm_K": t, "dh_fus_J_mol": h} for k, (t, h) in zip(iks, cc)}
    assert FM.svl_ideal_C(x, cc) == pytest.approx(ideal_sle_C(x, iks, tmap), abs=1e-9)


def test_activity_reduces_to_ideal_when_gamma_is_one():
    x = np.array([0.0017, 0.0104, 0.9879])
    cc = [C16, C13, IC8]
    one = lambda T, xx: np.ones(len(xx))                                   # noqa: E731
    assert FM.svl_activity_C(x, cc, one) == pytest.approx(FM.svl_ideal_C(x, cc), abs=1e-3)


def test_solid_solution_not_below_pure_solid_svl():
    # an ideal solid solution of n-alkanes always crystallises at or above the highest pure-solid SvL temperature of those alkanes
    x = [0.3, 0.3, 0.4]
    cc = [C16, C13, IC8]
    smi = ["C" * 16, "C" * 13, "CC(C)CC(C)(C)C"]
    assert FM.solid_solution_nalkanes_C(x, smi, cc) >= FM.svl_ideal_C(x, cc) - 1e-6


def test_missing_constant_does_not_crystallise():
    assert np.isnan(FM.svl_ideal_C([0.5, 0.5], [None, None]))
    assert FM.svl_ideal_C([0.5, 0.5], [C16, None]) == pytest.approx(FM.svl_ideal_C([0.5], [C16]))


def test_dedup_removes_oakley_copy_of_thermoml_point():
    import importlib.util
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "scripts/v5/base/delta_ws2_freeze_multicomp.py"
    spec = importlib.util.spec_from_file_location("ws2", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    smi = ["CCCCCC", "Cc1ccc(C)cc1", "CCCCCCCCCCCCCCCC"]
    iks = [m.ik(s) for s in smi]
    sid = m.system_id(iks)
    d = pd.DataFrame([dict(source="thermoml", system_id=sid, inchikeys=iks, x=[0.7036, 0.2004, 0.096], T_K=260.80),
                      dict(source="oakley2017", system_id=sid, inchikeys=iks, x=[0.7036, 0.2004, 0.0960], T_K=260.80),
                      dict(source="ascent065A", system_id="other", inchikeys=iks, x=[0.1, 0.1, 0.8], T_K=250.0)])
    kept, dropped = m.dedup(d)
    assert list(kept.source) == ["thermoml", "ascent065A"]
    assert list(dropped.source) == ["oakley2017"]
