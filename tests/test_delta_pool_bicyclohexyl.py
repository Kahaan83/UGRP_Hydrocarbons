"""R-BC-1 (results/v5/BASE/delta_pool_bicyclohexyl/README.md): the flash-registry additions overlay (on by default since decision 3, =0 restores the v2 registry)."""
import os

import pandas as pd

from fuelblend.baselines import flash_unifac as fu


import pytest


@pytest.fixture
def clean_registry_cache():
    """The registry is lru_cached: clear it before and after so the env var under test takes effect and does not leak."""
    fu._default_registry.cache_clear()
    yield
    fu._default_registry.cache_clear()


def test_overlay_on_by_default(monkeypatch, clean_registry_cache):
    """Decision 3 (2026-11-15): the additions overlay is merged unless the env var is "0"."""
    monkeypatch.delenv(fu.REGISTRY_ADDITIONS_ENV, raising=False)
    fu._default_registry.cache_clear()
    reg = fu._default_registry()
    src = reg["antoine_source"].astype(str)
    from rdkit import Chem
    for smi in ("C1CCC(C2CCCCC2)CC1", "C1CC[C@H]2CCCC[C@H]2C1"):
        can = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        assert src.loc[can] == "anchored_thermo_ambient_v5add"
        assert fu._lookup(smi, None) is not None


def test_overlay_off_equals_base(monkeypatch, clean_registry_cache):
    """Review Focus 2: with =0 the registry is the v2 base build, frame for frame."""
    monkeypatch.setenv(fu.REGISTRY_ADDITIONS_ENV, "0")
    fu._default_registry.cache_clear()
    reg = fu._default_registry()
    from rdkit import Chem
    base = pd.read_csv(fu.REGISTRY_CSV)
    base["can"] = base["smiles"].map(lambda t: Chem.MolToSmiles(Chem.MolFromSmiles(t)))
    assert reg.equals(base.drop_duplicates("can").set_index("can"))
    assert not (reg["antoine_source"].astype(str) == "anchored_thermo_ambient_v5add").any()
    assert fu._lookup("C1CCC(C2CCCCC2)CC1", None) is None
    assert fu._lookup("C1CC[C@H]2CCCC[C@H]2C1", None) is None


def test_overlay_on_adds_entries(monkeypatch, clean_registry_cache):
    monkeypatch.setenv(fu.REGISTRY_ADDITIONS_ENV, "1")
    fu._default_registry.cache_clear()
    r = fu._lookup("C1CCC(C2CCCCC2)CC1", None)
    assert r is not None and r[1] > 0
    assert fu._lookup("CCCCCCCCCC", None) is not None


def test_strict_constraints_guard_raises_on_dropped_constraint(clean_registry_cache):
    """R-BC-2 finding: a constraint whose model is unavailable for a pool member is dropped by _Sides; the guard refuses it.
    Pool member isododecane has rho15 status "disputed", so density / aromatics / smoke are unavailable (the original decalin_cis
    trigger, MW unverified, was resolved by decision 7)."""
    import pytest
    from fuelblend.design import blendstock_v5 as B
    from fuelblend.design import engine_v4 as E
    bp = B.BlendPool([B.pure_component(n) for n in ("isododecane", "n-decane", "tetralin")])
    models = B.build_models(bp)
    with pytest.raises(E.ModelUnavailable, match="dropped from the optimiser"):
        B.design_pool("astm:jet_a", bp, models, None, n_starts=1, seed=0, strict_constraints=True)


def test_strict_constraints_is_default():
    """Decision 6 (2026-11-15): a pool whose spec constraint is dropped by _Sides is refused unless the caller opts out."""
    import inspect
    from fuelblend.design import blendstock_v5 as B
    assert inspect.signature(B.design_pool).parameters["strict_constraints"].default is True
