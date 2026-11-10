"""v4 Task 16: set encoder invariances, zero-init = physics, structured T-form, locked-fuel contract, CPU smoke run."""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.mixture.set_encoder_v4 import T_REF, assert_no_locked, dev_group_folds, make_encoder  # noqa: E402

D = 7


def _batch(B=5, N=6, seed=0):
    g = torch.Generator().manual_seed(seed)
    comp = torch.randn(B, N, D, generator=g)
    w = torch.rand(B, N, generator=g) + 0.05
    w = w / w.sum(1, keepdim=True)
    T = 230.0 + 80.0 * torch.rand(B, generator=g)
    ctx = torch.ones(B, 1)
    return comp, w, T, ctx


def _randomise(model):
    """Give the zero-initialised heads non-trivial weights."""
    torch.manual_seed(1)
    for h in model.heads.values():
        torch.nn.init.normal_(h.out.weight, std=0.5)
        torch.nn.init.normal_(h.out.bias, std=0.5)
    return model


@pytest.mark.parametrize("pool", ["attention", "sum"])
def test_permutation_invariance(pool):
    torch.manual_seed(0)
    m = _randomise(make_encoder(D, pool, hidden=16, d_context=1)).eval()
    comp, w, T, ctx = _batch()
    base = torch.zeros(5)
    ref = m("viscosity", comp, w, T, base, ctx)
    perm = torch.randperm(comp.shape[1])
    out = m("viscosity", comp[:, perm], w[:, perm], T, base, ctx)
    assert torch.allclose(ref, out, atol=1e-6)
    assert ref.abs().max() > 1e-4                         # the test is not vacuous


@pytest.mark.parametrize("pool", ["attention", "sum"])
def test_zero_fraction_component_no_effect(pool):
    torch.manual_seed(0)
    m = _randomise(make_encoder(D, pool, hidden=16, d_context=1)).eval()
    comp, w, T, ctx = _batch()
    base = torch.zeros(5)
    ref = m("flash", comp, w, T, base, ctx)
    # append a wild extra component with zero weight
    comp2 = torch.cat([comp, 100.0 * torch.randn(5, 2, D)], dim=1)
    w2 = torch.cat([w, torch.zeros(5, 2)], dim=1)
    out = m("flash", comp2, w2, T, base, ctx)
    assert torch.allclose(ref, out, atol=1e-6)
    # zeroing an existing component equals deleting it (after renormalisation is NOT applied: weights are used as given)
    w3 = w.clone(); w3[:, 0] = 0.0
    a = m("viscosity", comp, w3, T, base, ctx)
    b = m("viscosity", comp[:, 1:], w3[:, 1:], T, base, ctx)
    assert torch.allclose(a, b, atol=1e-6)


@pytest.mark.parametrize("pool", ["attention", "sum"])
def test_residual_zero_init_returns_physics(pool):
    torch.manual_seed(3)
    m = make_encoder(D, pool, hidden=16, d_context=1).eval()
    comp, w, T, ctx = _batch()
    base = torch.randn(5) * 3.0
    for head in ("viscosity", "flash"):
        assert torch.equal(m(head, comp, w, T, base, ctx), base), head
        assert float(m.residual(head, comp, w, T, ctx).abs().max()) == 0.0


def test_head_output_linear_in_1000_over_T():
    torch.manual_seed(2)
    m = _randomise(make_encoder(D, "attention", hidden=16, d_context=1)).eval()
    comp, w, _, ctx = _batch(B=4)
    invT = torch.linspace(2.5, 5.5, 7)                    # 400 K ... 182 K, inc. well outside the training range
    outs = []
    for it in invT:
        T = torch.full((4,), 1000.0 / float(it))
        outs.append(m.residual("viscosity", comp, w, T, ctx).detach())
    outs = torch.stack(outs)                              # [7,4]
    second = outs[2:] - 2 * outs[1:-1] + outs[:-2]
    assert float(second.abs().max()) < 1e-5 * max(1.0, float(outs.abs().max()))
    slope = (outs[-1] - outs[0]) / (invT[-1] - invT[0])
    assert float(slope.abs().max()) > 1e-3                # genuinely T-dependent
    a0, a1 = m.heads["viscosity"].coefficients(m.mixture_state(comp, w, ctx))
    pred = m.heads["viscosity"].scale * (a0 + a1 * (1000.0 / torch.full((4,), 300.0) - T_REF))
    assert torch.allclose(pred, m.residual("viscosity", comp, w, torch.full((4,), 300.0), ctx), atol=1e-6)
    # scalar head has no T dependence at all
    f = [m.residual("flash", comp, w, torch.full((4,), t), ctx) for t in (200.0, 300.0, 400.0)]
    assert torch.allclose(f[0], f[1]) and torch.allclose(f[1], f[2])


def test_locked_fuels_never_in_training_or_calibration_indices():
    fuels = np.array([f"F{i}" for i in range(30) for _ in range(3)])
    locked = {"F3", "F7", "F11", "F19", "F25"}
    folds = dev_group_folds(fuels, locked, n_folds=7, seed=1)
    seen = []
    for tr, te in folds:
        assert_no_locked({"train": tr, "test": te}, fuels, locked)
        assert not (set(fuels[tr]) & set(fuels[te]))
        seen.extend(te.tolist())
    assert sorted(seen) == sorted(np.where(~np.isin(fuels, list(locked)))[0].tolist())   # every dev row tested once
    with pytest.raises(AssertionError):
        assert_no_locked({"calibration": np.where(fuels == "F3")[0]}, fuels, locked)
    # real split file: the 12 locked CRC fuels are excluded for the real inventory
    sp = json.load(open(ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"))
    assert len(sp["locked"]) >= 12


def test_cpu_smoke_run(tmp_path):
    if not (ROOT / "data" / "unified" / "compounds.csv").exists() or not (ROOT / "results" / "v4" / "F" / "flash_lab_loso.csv").exists():
        pytest.skip("local data not present")
    t0 = time.time()
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "v4" / "train_set_encoder.py"), "--epochs", "2", "--seeds", "1",
                        "--lab-folds", "0", "--pools", "attention", "--out", str(tmp_path), "--device", "cpu"],
                       capture_output=True, text=True, timeout=240)
    dt = time.time() - t0
    assert r.returncode == 0, r.stderr[-2000:]
    assert dt < 120, f"smoke run took {dt:.0f}s"
    audit = json.load(open(tmp_path / "split_audit.json"))
    assert audit["locked_in_any_train_or_calibration_index"] == 0
    assert audit["n_outer_folds"] == 9 and audit["n_real_visc_rows_locked"] == 36
    for f in ("viscosity_encoder_vs_task11.csv", "conformal_calibration.csv", "flash_encoder_real_d93.csv", "promotion_decision.json"):
        assert (tmp_path / f).exists()
