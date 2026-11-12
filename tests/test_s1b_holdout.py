"""v5 S1b: synthetic checks of the held-out test helpers (no real scoring)."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("s1b", ROOT / "scripts" / "v5" / "eval_lab_thermoml2020_holdout_v5.py")
s1b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s1b)


def test_binomial_threshold_n17_is_13():
    assert s1b.binom_threshold(17) == 13


def test_binomial_threshold_is_smallest_k_with_cdf_at_least_5_percent():
    from scipy.stats import binom
    for n in (9, 10, 25, 40):
        k = s1b.binom_threshold(n)
        assert binom.cdf(k, n, 0.9) >= 0.05
        assert k == 0 or binom.cdf(k - 1, n, 0.9) < 0.05


def test_verdict_not_evaluable_below_9_systems():
    k, verdict = s1b.coverage_verdict(8, 8)
    assert verdict == "NOT EVALUABLE" and k is None
    assert s1b.coverage_verdict(0, 0)[1] == "NOT EVALUABLE"


def test_verdict_holds_and_fails_at_threshold():
    assert s1b.coverage_verdict(13, 17) == (13, "HOLDS")
    assert s1b.coverage_verdict(17, 17) == (13, "HOLDS")
    assert s1b.coverage_verdict(12, 17) == (13, "FAILS")


def test_select_test_systems_excludes_ids_in_the_split():
    split = {"dev": ["a", "b"], "locked": ["c"]}
    test, not_indep = s1b.select_test_systems(["a", "c", "n1", "n2", "n1"], split)
    assert test == ["n1", "n2"]
    assert not_indep == ["a", "c"]
    assert not set(test) & (set(split["dev"]) | set(split["locked"]))


def test_adamantane_detector():
    from rdkit import Chem
    for smi in ("C1C2CC3CC1CC(C2)C3", "CC12CC3CC(CC(C3)C1)C2"):
        assert Chem.MolFromSmiles(smi) is not None
        assert s1b.has_adamantane(smi)
    assert not s1b.has_adamantane("C1CCC2CCCCC2C1")     # decalin
    assert not s1b.has_adamantane("CCCCCCCCCC")          # n-decane
    assert not s1b.has_adamantane("not-a-smiles")
    assert not s1b.has_adamantane(None)
