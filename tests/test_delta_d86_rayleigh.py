"""R-DR-1 (results/v5/BASE/delta_d86_rayleigh/README.md): Rayleigh D86 model unit tests."""
import math

import numpy as np
import pytest

from fuelblend.baselines.distillation_rayleigh import bubble_T_K, rayleigh_d86, rayleigh_quantities

P = 101.325


def _tb(A, B, C):
    return B / (A - math.log10(P)) - C - 273.15


def test_single_compound_all_cuts_equal_tb():
    A, B, C = 6.44, 1707.53, -60.28  # n-decane, flash registry
    r = rayleigh_d86([A], [B], [C], [1.0], [195.0])
    for k in ("IBP", "T10", "T50", "T90", "T99"):
        assert r[k] == pytest.approx(_tb(A, B, C), abs=1e-6)


def test_binary_constant_relative_volatility_matches_analytic_rayleigh():
    # same B and C -> P2/P1 = 10**(A2-A1) constant: alpha = 10**(A1-A2) for component 1 (light)
    A1, A2, Bc, Cc = 6.8, 6.5, 1800.0, -60.0
    alpha = 10 ** (A1 - A2)
    x0 = 0.4
    r = rayleigh_d86([A1, A2], [Bc, Bc], [Cc, Cc], [x0, 1 - x0], [1.0, 1.0], cuts_pct=(10.0, 50.0, 90.0))
    # analytic: ln(L/L0) = 1/(alpha-1) [ ln(x/x0) + alpha ln((1-x0)/(1-x)) ];  equal molar volumes -> V_d/V0 = 1 - L/L0
    from scipy.optimize import brentq

    def L_of_x(x):
        return math.exp((math.log(x / x0) + alpha * math.log((1 - x0) / (1 - x))) / (alpha - 1))

    for q in (10.0, 50.0, 90.0):
        x = brentq(lambda xx: (1 - L_of_x(xx)) - q / 100.0, 1e-12, x0 - 1e-12)
        T = bubble_T_K(np.array([x, 1 - x]), np.array([A1, A2]), np.array([Bc, Bc]), np.array([Cc, Cc])) - 273.15
        assert r[f"T{int(q)}"] == pytest.approx(T, abs=0.05)


def test_cuts_monotone_multicomponent():
    A = [5.9603, 6.44, 6.7481, 6.9108]
    B = [1270.87, 1707.53, 1996.0, 2153.31]
    C = [-50.98, -60.28, -63.17, -62.46]
    q = rayleigh_quantities(A, B, C, [0.05, 0.4, 0.3, 0.25], [114.2, 142.3, 170.3, 184.4], [696.0, 734.7, 753.7, 759.9])
    assert q["IBP"] <= q["T10"] <= q["T50"] <= q["T90"] <= q["FBP"]
    assert q["T90_T10"] > 0


def test_add_d86_models_default_is_engler_unchanged():
    from fuelblend.design import blendstock_v5 as B
    import inspect
    sig = inspect.signature(B.add_d86_models)
    assert sig.parameters["d86_model"].default == "engler"


def test_predict_fuel_reports_d86_block():
    import warnings
    from fuelblend.api.predict_route import predict_fuel
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = predict_fuel({"CCCCCCCCCC": 0.4, "CCCCCCCCCCCC": 0.3, "Cc1ccc(C)c(C)c1": 0.2, "CC(C)CC(C)(C)C": 0.1})
    d = r["distillation_D86_C"]
    assert d["status"] == "ok" and d["method"] == "engler_tbp_offsets" and d["adopted_in_uncertainty_yaml"] is False
    v = d["value"]
    assert v["T10"] <= v["T50"] <= v["T90"]
    rv = d["variants"]["rayleigh"]["value"]
    assert rv["T10"] <= rv["T50"] <= rv["T90"] <= rv["FBP"]
    assert any(w.startswith("d86_models_disagree") for w in d["warnings"])
