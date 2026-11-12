"""S7 viscosity pairwise bootstrap / conformal tests (small)."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("s7", ROOT / "scripts" / "v5" / "s7_viscosity_pairwise.py")
s7 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s7)


def test_bootstrap_self_difference_is_zero():
    rng = np.random.default_rng(1)
    e = rng.normal(size=60)
    fid = np.repeat(np.arange(20), 3)
    d, lo, hi = s7.paired_boot(e, e, fid, n=200, seed=0)
    assert d == 0.0 and lo == 0.0 and hi == 0.0


def test_conformal_factor_is_order_statistic():
    rng = np.random.default_rng(2)
    for n in (9, 19, 45, 100):
        sc = rng.uniform(size=n)
        k = math.ceil((n + 1) * 0.9)
        assert s7.conformal_factor(sc) == np.sort(sc)[k - 1]
    assert s7.conformal_factor(np.arange(1.0, 6.0)) == np.inf       # n=5 -> k=6 > n


def test_locked_fuels_never_enter_factor():
    rng = np.random.default_rng(3)
    df = pd.DataFrame(dict(fuel_id=[f"f{i}" for i in range(30) for _ in range(3)],
                           T_C=[20.0, -20.0, -40.0] * 30, r=rng.uniform(size=90)))
    df["locked"] = df.fuel_id.isin(["f0", "f1", "f2"])
    base = {c: s7.dev_factor(df, "r", c) for c in (None, 20.0, -40.0)}
    df2 = df.copy()
    df2.loc[df2.locked, "r"] = 1e6
    assert {c: s7.dev_factor(df2, "r", c) for c in (None, 20.0, -40.0)} == base
    assert base[None][1] == 27
