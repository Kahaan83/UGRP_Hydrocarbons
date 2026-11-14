"""WS5 (GH-GNN activity coefficients in Le Chatelier flash): unit tests of the wrapper pieces, written before the driver."""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "v5" / "ws5"))

import wilson as W  # noqa: E402


def test_wilson_inversion_round_trip():
    rng = np.random.default_rng(0)
    a = rng.uniform(-0.6, 1.5, 200)
    b = rng.uniform(-0.6, 1.5, 200)
    u, v, ok = W.wilson_invert(a, b)
    assert ok.mean() > 0.9
    f1, f2 = 1 - u - np.exp(v), 1 - v - np.exp(u)
    assert np.allclose(f1[ok], a[ok], atol=1e-7) and np.allclose(f2[ok], b[ok], atol=1e-7)


def test_wilson_failure_rule_sets_lambda_one():
    # a, b chosen with no real Wilson solution: very large negative inf-dil coefficients in both directions
    u, v, ok = W.wilson_invert(np.array([-6.0]), np.array([-6.0]))
    if not ok[0]:
        assert u[0] == 0.0 and v[0] == 0.0


def test_wilson_recovers_infinite_dilution_and_ideal():
    a, b = 0.35, -0.12
    u, v, ok = W.wilson_invert(np.array([a]), np.array([b]))
    assert ok[0]
    lnL = np.array([[0.0, u[0]], [v[0], 0.0]])
    g1 = W.wilson_gamma([1e-9, 1 - 1e-9], lnL)[0]          # component 1 at infinite dilution in 2
    g2 = W.wilson_gamma([1 - 1e-9, 1e-9], lnL)[1]
    assert np.log(g1) == pytest.approx(a, abs=1e-6) and np.log(g2) == pytest.approx(b, abs=1e-6)
    assert np.allclose(W.wilson_gamma([0.3, 0.7], np.zeros((2, 2))), 1.0)
    # pure component limits: gamma_i -> 1 as x_i -> 1
    assert W.wilson_gamma([1 - 1e-9, 1e-9], lnL)[0] == pytest.approx(1.0, abs=1e-6)


def test_margules_limits():
    g = W.margules_gamma_binary([1e-9, 1 - 1e-9], 0.4, -0.2)
    assert np.log(g[0]) == pytest.approx(0.4, abs=1e-6)
    g = W.margules_gamma_binary([1 - 1e-9, 1e-9], 0.4, -0.2)
    assert np.log(g[1]) == pytest.approx(-0.2, abs=1e-6)


@pytest.fixture(scope="module")
def lib():
    pytest.importorskip("torch")
    if not (ROOT / "data" / "raw" / "ghgnn_code" / "repo" / "notebooks" / "GHGNN" / "GHGNN.pth").exists():
        pytest.skip("GH-GNN clone not present (data/raw/ghgnn_code, see its SOURCES.md)")
    warnings.simplefilter("ignore")
    import ghgnn_lib
    return ghgnn_lib


def test_network_is_A_plus_B_over_T(lib):
    sol = ["c1ccccc1", "CCCCCCCCCC", "Cc1cc(C)cc(C)c1", "C1CCCCC1"]
    slv = ["CCCCCCCCCCCC", "c1ccccc1", "CCCCCCCCCC", "CCCCCCCCCCCCCCCC"]
    l25 = lib.predict_ln_gamma_inf(sol, slv, [25.0] * 4)
    l75 = lib.predict_ln_gamma_inf(sol, slv, [75.0] * 4)
    l50 = lib.predict_ln_gamma_inf(sol, slv, [50.0] * 4)
    T1, T2, T3 = 298.15, 348.15, 323.15
    B = (l25 - l75) / (1 / T1 - 1 / T2)
    A = l25 - B / T1
    assert np.allclose(A + B / T3, l50, atol=1e-4)


def test_ad_table_stores_real_values_for_interpolation_pairs(lib):
    r = lib.ad_pair("c1ccccc1", "CCCCCCCCCCCC")
    assert r["interpolation"] is True
    assert r["tanimoto10"] < 9000 and r["n_class_pair"] < 9000 and r["n_class_pair"] > 0
    assert r["ad_feasible"] is True


def test_gamma_one_wrapper_equals_repo_ideal_flash():
    import ws5_flash as F
    from fuelblend.baselines.flash_unifac import flash_point_unifac_detail
    lab = pd.read_csv(ROOT / "results" / "v4" / "F" / "flash_lab_loso.csv").iloc[::17].head(8)
    for r in lab.itertuples():
        smi, x = json.loads(r.comp_smiles), json.loads(r.x_mole)
        ref = flash_point_unifac_detail(smi, x, use_unifac=False)["flash_C"]
        got, status = F.flash_with_gamma(smi, x, None, gamma_fn=None)
        assert status == "ok" and got == pytest.approx(ref, abs=1e-9)
