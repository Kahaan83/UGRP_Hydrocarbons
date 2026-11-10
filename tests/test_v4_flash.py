"""v4 Task 9: UNIFAC-corrected Le Chatelier flash point (src/fuelblend/baselines/flash_unifac.py)."""
import warnings

import numpy as np
import pytest

from fuelblend.baselines.flash_unifac import (
    UnifacFallbackWarning,
    flash_point_ideal,
    flash_point_unifac,
    flash_point_unifac_detail,
    unifac_subgroups,
)

DECANE, DODECANE = "CCCCCCCCCC", "CCCCCCCCCCCC"
TOLUENE = "Cc1ccccc1"


def test_unifac_reduces_to_ideal_for_alkane_pair():
    x = np.array([0.5, 0.5])
    ideal = flash_point_ideal([DECANE, DODECANE], x)
    uni = flash_point_unifac([DECANE, DODECANE], x)
    assert np.isfinite(ideal) and np.isfinite(uni)
    assert abs(uni - ideal) < 0.5


def test_unifac_lowers_flash_for_aromatic_alkane():
    x = np.array([0.5, 0.5])
    d = flash_point_unifac_detail([TOLUENE, DODECANE], x)
    # premise: toluene/alkane is a positive-deviation pair (gamma > 1 for both at the flash point)
    assert d["gamma_at_flash"].min() > 1.0
    assert d["flash_C"] <= flash_point_ideal([TOLUENE, DODECANE], x)


def test_group_assignment_original_unifac():
    assert unifac_subgroups(DECANE) == {1: 2, 2: 8}
    assert unifac_subgroups(TOLUENE) == {9: 5, 11: 1}
    assert unifac_subgroups("CCc1ccccc1") == {1: 1, 9: 5, 12: 1}
    assert unifac_subgroups("CC(C)c1ccccc1") == {1: 2, 9: 5, 13: 1}
    assert unifac_subgroups("C1CCCCC1") == {2: 6}
    assert unifac_subgroups("c1ccc2ccccc2c1") == {9: 8, 10: 2}
    assert unifac_subgroups("CC(C)(C)C") == {1: 4, 4: 1}


def test_deuterated_alkane_gets_alkane_groups():
    # the v4 CRC647 reconstructed compositions contain perdeuterated alkane SMILES; H and D group identically
    assert unifac_subgroups("[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])[2H]") == {1: 2, 2: 1}


def test_no_groups_falls_back_to_gamma_one_and_is_flagged():
    # alkyne has no subgroup in the hydrocarbon fragmenter -> documented fallback, never silent
    assert unifac_subgroups("C#CCCCC") is None
    x = np.array([0.5, 0.5])
    with pytest.warns(UnifacFallbackWarning):
        d = flash_point_unifac_detail(["C#CCCCC", DODECANE], x)
    assert d["fallback_components"] == [0]
    assert np.isfinite(d["flash_C"]) or d["status"] != "ok"


def test_ideal_single_component_equals_unifac():
    x = np.array([1.0])
    assert flash_point_unifac([DODECANE], x) == pytest.approx(flash_point_ideal([DODECANE], x), abs=1e-6)
