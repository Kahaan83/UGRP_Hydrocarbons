"""v5: Boehm et al. 2022 viscosity loader.  Fixture tests/data/boehm2022_sample.xlsx = 16 rows copied verbatim from the
real 'This Work' sheet (blocks A, series h, B, and the first rows of C incl. a POSF fuel block)."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.data import loaders_boehm2022 as L

FIX = Path(__file__).parent / "data" / "boehm2022_sample.xlsx"
CH = L._chem()
IK = {n: v[0] for n, v in CH.items()}


@pytest.fixture(scope="module")
def built():
    log = {}
    return L.build(FIX, log), log


def pick(d, prop, iks, T_K):
    iks = sorted(iks)
    m = d[(d.property == prop) & (np.abs(d.T_K - T_K) < 0.01) & (d.comp_inchikeys == json.dumps(iks))]
    assert len(m) == 1, (prop, iks, T_K, len(m))
    return m.iloc[0]


def test_schema_and_fuel_rows_dropped(built):
    d, log = built
    assert list(d.columns) == L.COLS
    assert log["rows_read"] == 16 and log["rows_kept"] == 12      # 4 rows touch A-1 (POSF 10264)
    assert log["dropped_fuel_or_blend_stock_components"] == ["A-1 (POSF 10264)"]
    assert d.is_hc.all() and not d.conflict.any() and (d.sources == "boehm2022").all()
    assert set(d.property) == {"viscosity_kin_cSt", "viscosity_dyn_mPas", "density_kg_m3"}
    # pure 2-methylpentane / ethylcyclohexane rows of the fuel block survive (fuel fraction is 0)
    assert IK["2-methylpentane"] in " ".join(d.comp_inchikeys)


def test_composition_conversion(built):
    d, _ = built
    for r in d.itertuples():
        assert abs(sum(json.loads(r.x_mole)) - 1) < 1e-5 and abs(sum(json.loads(r.w_mass)) - 1) < 1e-5
        assert r.system_id == hashlib.md5("|".join(sorted(json.loads(r.comp_inchikeys))).encode()).hexdigest()[:12]
    # 75/25 vol HMN/MCH with the lab's measured 15 C pure densities (0.78772, 0.77358)
    r = pick(d, "viscosity_kin_cSt", [IK["2,2,4,4,6,8,8-heptamethylnonane"], IK["methylcyclohexane"]], 233.15)
    iks = json.loads(r.comp_inchikeys)
    w = np.array(json.loads(r.w_mass)); x = np.array(json.loads(r.x_mole))
    m = {IK["2,2,4,4,6,8,8-heptamethylnonane"]: 0.75 * 0.78772, IK["methylcyclohexane"]: 0.25 * 0.77358}
    wexp = np.array([m[k] for k in iks]); wexp /= wexp.sum()
    mw = {IK[n]: CH[n][1] for n in CH}
    xexp = wexp / np.array([mw[k] for k in iks]); xexp /= xexp.sum()
    assert np.allclose(w, wexp, atol=1e-5) and np.allclose(x, xexp, atol=1e-5)


def test_ten_value_spot_check(built):
    d, _ = built
    hmn, mch, hexb = IK["2,2,4,4,6,8,8-heptamethylnonane"], IK["methylcyclohexane"], IK["hexylbenzene"]
    eb, bch, hpb = IK["ethylbenzene"], IK["butylcyclohexane"], IK["heptylbenzene"]
    checks = [
        ("viscosity_dyn_mPas", [hmn, mch], 233.15, 15.77 * 0.82421),
        ("viscosity_kin_cSt", [hmn, mch], 233.15, 15.77),
        ("density_kg_m3", [hmn, mch], 233.15, 824.21),
        ("viscosity_dyn_mPas", [hmn, mch], 288.15, 3.0183 * 0.78441),
        ("viscosity_kin_cSt", [mch], 288.15, 0.9767),
        ("viscosity_dyn_mPas", [hexb], 253.15, 5.1112 * 0.88823),
        ("viscosity_dyn_mPas", [eb, bch, hpb], 288.15, 0.91158 * 0.86073),
        ("density_kg_m3", [eb, bch, hpb], 253.15, 889.68),
        ("viscosity_kin_cSt", [IK["2-methylpentane"]], 233.15, 0.80418),
        ("viscosity_kin_cSt", [eb], 273.15, 0.92161),
    ]
    for prop, iks, T, exp in checks:
        assert abs(pick(d, prop, iks, T).value - exp) < 1e-4 * max(1, abs(exp)), (prop, T, exp)
