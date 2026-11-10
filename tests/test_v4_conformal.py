"""v4 Task 10: conformal margins from group-level residuals."""
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from fuelblend.uq.conformal import coverage, one_sided_margin, two_sided_margin

ROOT = Path(__file__).resolve().parents[1]


def _nine_groups():
    # group g has rows with residuals up to g; the max per group is 1..9
    res, grp = [], []
    for g in range(1, 10):
        res += [g - 0.5, float(g), g - 1.0]
        grp += [f"g{g}"] * 3
    return np.array(res), np.array(grp)


def test_margin_matches_hand_example():
    r, g = _nine_groups()
    assert one_sided_margin(r, g, alpha=0.1) == 9.0   # n=9 -> k=ceil(10*0.9)=9 -> 9th smallest


def test_margin_picks_quantile_rank():
    # n=19 groups, alpha=0.1 -> k=ceil(20*0.9)=18 -> 18th smallest = 18
    r = np.arange(1, 20, dtype=float)
    g = np.array([f"g{i}" for i in range(19)])
    assert one_sided_margin(r, g, alpha=0.1) == 18.0


@pytest.mark.parametrize("n", [0, 1, 5, 8])
def test_margin_infinite_when_too_few_groups(n):
    r = np.arange(1, n + 1, dtype=float)
    g = np.array([f"g{i}" for i in range(n)])
    with pytest.warns(RuntimeWarning, match="too few"):
        assert math.isinf(one_sided_margin(r, g, alpha=0.1))


def test_no_warning_when_enough_groups():
    r, g = _nine_groups()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        one_sided_margin(r, g, alpha=0.1)


def test_margin_uses_groups_not_rows():
    r, g = _nine_groups()
    base = one_sided_margin(r, g, alpha=0.1)
    mask = g == "g3"
    r2 = np.concatenate([r] + [r[mask]] * 50)
    g2 = np.concatenate([g] + [g[mask]] * 50)
    assert one_sided_margin(r2, g2, alpha=0.1) == base
    # an outlier repeated in many rows of one group still counts as one group
    r3 = np.r_[r, [100.0] * 40]
    g3 = np.r_[g, ["g9"] * 40]
    assert one_sided_margin(r3, g3, alpha=0.1) == 100.0   # n=9 -> k=9 -> the largest group score


def test_direction_over_vs_under():
    r = np.array([5.0, -2.0, 3.0, -7.0, 1.0, -1.0, 2.0, -3.0, 4.0, -4.0])
    g = np.array([f"g{i}" for i in range(10)])
    assert one_sided_margin(r, g, alpha=0.1, direction="over") == 5.0
    assert one_sided_margin(r, g, alpha=0.1, direction="under") == 7.0
    with pytest.raises(ValueError):
        one_sided_margin(r, g, direction="sideways")


def test_two_sided_uses_max_abs_per_group():
    r = np.array([-9.0, 1.0] + [float(i) for i in range(1, 9)])
    g = np.array(["a", "a"] + [f"g{i}" for i in range(8)])
    assert two_sided_margin(r, g, alpha=0.1) == 9.0   # group a max|r| = 9; n=9, k=9
    assert one_sided_margin(r, g, alpha=0.1, direction="over") == 8.0


def test_margin_can_be_negative_when_never_unsafe():
    r = -np.arange(1, 21, dtype=float)
    g = np.array([f"g{i}" for i in range(20)])
    assert one_sided_margin(r, g, alpha=0.1, direction="over") < 0


def test_coverage_counts_groups():
    r, g = _nine_groups()
    assert coverage(r, g, 5.0) == (5, 9)


def _load_yaml():
    import yaml
    p = ROOT / "configs" / "uncertainty_v4.yaml"
    if not p.exists():
        pytest.skip("configs/uncertainty_v4.yaml not generated")
    text = p.read_text(encoding="utf-8")
    return text, yaml.safe_load(text)


def test_generated_yaml_header_and_consistency():
    text, doc = _load_yaml()
    assert text.splitlines()[0].startswith("# GENERATED - do not edit")
    assert "CONSUMERS MUST read applied_margin" in text and "NO VALID MARGIN" in text
    assert doc["generated_by"] == "scripts/v4/build_uncertainty_v4.py"
    for key, e in doc["margins"].items():
        for fld in ("value", "n_groups", "residual_file", "model", "direction", "alpha"):
            assert fld in e, (key, fld)
        k = math.ceil(round((e["n_groups"] + 1) * (1 - e["alpha"]), 9))
        if k > e["n_groups"]:
            assert e["value"] == "inf" and e.get("reason"), key   # no hand-typed substitute
        else:
            assert isinstance(e["value"], float), key
        assert "dvpe" not in key.lower()
        assert (ROOT / e["residual_file"]).exists(), key
        assert "recompute" in e or "residual_derivation" in e or key.startswith("flash_point_C.real_fuels_D56"), key


def test_applied_margin_contract():
    _, doc = _load_yaml()
    for key, e in doc["margins"].items():
        if e["value"] == "inf":
            assert "applied_margin" not in e and e.get("reason"), key
            continue
        assert "applied_margin" in e, key
        if e["direction"] == "two_sided":
            assert e["applied_margin"] == e["value"] >= 0, key
        else:
            assert e["applied_margin"] >= 0, key
            assert e["applied_margin"] == pytest.approx(max(0.0, e["value"])), key
            if e["value"] < 0:
                assert e["applied_margin"] == 0.0 and "applied_margin_note" in e, key


def test_flash_guidance_present():
    _, doc = _load_yaml()
    for k in ("flash_point_C.lab_mixtures.ideal", "flash_point_C.real_fuels_D93_reconstructed.ideal"):
        assert "NOT use 0" in doc["margins"][k]["recommendation"]


def test_yaml_values_recomputed_from_residual_files():
    import pandas as pd
    _, doc = _load_yaml()
    n_checked = 0
    for key, e in doc["margins"].items():
        spec = e.get("recompute")
        if not spec:
            continue
        df = pd.read_csv(ROOT / e["residual_file"])
        for col, op, val in spec["where"]:
            if op == "eq":
                df = df[df[col] == val]
            elif op == "ne":
                df = df[df[col] != val]
            elif op == "close":
                df = df[np.isclose(df[col], val)]
            else:
                raise AssertionError(op)
        r, g = df[spec["column"]].to_numpy(float), df[spec["group"]].to_numpy()
        if e["direction"] == "two_sided":
            got = two_sided_margin(r, g, e["alpha"])
        else:
            got = one_sided_margin(r, g, e["alpha"], e["direction"])
        want = e["value"]
        if want == "inf":
            assert math.isinf(got), key
        else:
            assert got == pytest.approx(want, rel=1e-6, abs=1e-6), key
        n_checked += 1
    assert n_checked >= 20
