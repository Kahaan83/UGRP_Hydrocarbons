"""v5 delta_smoke_cert 1.4: real-fuel smoke margin (R-SC-2) and class-mix gate. Synthetic frames test the margin math."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.design.smoke_margin_v5 import class_mix_in_real_fuel_range, real_fuel_class_ranges

_p = Path(__file__).resolve().parents[1] / "scripts" / "v5" / "base" / "delta_smoke_cert.py"
_spec = importlib.util.spec_from_file_location("delta_smoke_cert", _p)
dsc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dsc)
fit_smoke_margin = dsc.fit_smoke_margin


def _synth(sizes, noise=None, bias=None, seed=0):
    rng = np.random.default_rng(seed)
    rows, fid = [], 0
    for g, n in sizes.items():
        s = (noise or {}).get(g, 1.0)
        b = (bias or {}).get(g, 0.0)
        for _ in range(n):
            a = rng.uniform(9, 22)
            rows.append(dict(fuel_id=fid, source_group=g, aromatics_vol_pct=a,
                             smoke_point_mm=26 - 0.3 * a + b + rng.normal(0, s), smoke_qualifier=""))
            fid += 1
    return pd.DataFrame(rows)


def test_censored_excluded():
    df = _synth({f"G{i}": 60 for i in range(5)})
    cens = df.iloc[:7].copy()
    cens["smoke_qualifier"] = ">"
    cens["smoke_point_mm"] = np.nan
    cens["fuel_id"] += 10_000
    out = fit_smoke_margin(pd.concat([df, cens], ignore_index=True))
    assert out["n_censored_excluded"] == 7
    assert out["n_rows"] == 300


def test_gate_too_few():
    d150 = _synth({f"G{i}": 30 for i in range(5)})            # 150 fuels
    assert fit_smoke_margin(d150)["margin"] == math.inf
    d4 = _synth({f"G{i}": 80 for i in range(4)} | {"S": 20})  # only 4 groups with >= 40
    assert fit_smoke_margin(d4)["margin"] == math.inf


def test_direction_over():
    df = _synth({f"G{i}": 60 for i in range(5)})
    out = fit_smoke_margin(df, alpha=0.10)
    r = np.sort(np.asarray(out["oof_residuals"]))
    n = len(r)
    k = math.ceil((n + 1) * 0.9)
    assert out["direction"] == "over"
    assert out["fuel_level_margin"] == r[k - 1]
    # residual = pred - measured: a model that under-predicts by construction would give negative residuals
    df2 = df.copy()
    df2["smoke_point_mm"] = df2["smoke_point_mm"] + 5.0  # measured higher -> pred - meas shifts down
    assert fit_smoke_margin(df2)["fuel_level_margin"] < out["fuel_level_margin"] + 1e-9


def test_certifying_margin_is_max():
    df = _synth({"A": 80, "B": 80, "C": 80, "D": 80, "E": 80, "F": 10},
                noise={"E": 3.0, "F": 9.0})
    out = fit_smoke_margin(df)
    elig = {g: m for g, m in out["per_group_margin"].items() if out["per_group_n"][g] >= 40}
    assert "F" not in elig and len(elig) == 5
    assert out["margin"] == max(out["fuel_level_margin"], max(elig.values()))
    assert max(elig.values()) > out["fuel_level_margin"]  # group E dominates here
    assert 0.0 <= out["coverage_oof"] <= 1.0


def test_class_mix_gate():
    lo_hi = real_fuel_class_ranges()
    assert lo_hi["cyclo"][1] < 0.75
    ok, reason = class_mix_in_real_fuel_range(0.75, 0.10, 0.15)
    assert not ok and reason == "class_mix_outside_real_fuel_range"
    mid = {k: (v[0] + v[1]) / 2 for k, v in lo_hi.items()}
    ok, reason = class_mix_in_real_fuel_range(mid["cyclo"], mid["nparaffin"], mid["isoparaffin"])
    assert ok and reason is None
