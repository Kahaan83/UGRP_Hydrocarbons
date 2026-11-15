"""R-BC-1 (results/v5/BASE/delta_pool_bicyclohexyl/README.md): the flash-registry additions overlay is off by default."""
import os

import pandas as pd

from fuelblend.baselines import flash_unifac as fu


def test_overlay_off_by_default_and_no_entry():
    assert os.environ.get(fu.REGISTRY_ADDITIONS_ENV) != "1"
    fu._default_registry.cache_clear()
    reg = fu._default_registry()
    from rdkit import Chem
    base = pd.read_csv(fu.REGISTRY_CSV)
    base["can"] = base["smiles"].map(lambda t: Chem.MolToSmiles(Chem.MolFromSmiles(t)))
    assert reg.equals(base.drop_duplicates("can").set_index("can"))   # off: the frame is exactly the base build
    assert not (reg["antoine_source"].astype(str) == "anchored_thermo_ambient_v5add").any()
    assert fu._lookup("C1CCC(C2CCCCC2)CC1", None) is None
    assert fu._lookup("C1CC[C@H]2CCCC[C@H]2C1", None) is None


def test_overlay_on_adds_entries(monkeypatch):
    monkeypatch.setenv(fu.REGISTRY_ADDITIONS_ENV, "1")
    fu._default_registry.cache_clear()
    try:
        r = fu._lookup("C1CCC(C2CCCCC2)CC1", None)
        assert r is not None and r[1] > 0
        assert fu._lookup("CCCCCCCCCC", None) is not None
    finally:
        monkeypatch.delenv(fu.REGISTRY_ADDITIONS_ENV)
        fu._default_registry.cache_clear()


def test_strict_constraints_guard_raises_on_dropped_constraint(monkeypatch):
    """R-BC-2 finding: a constraint whose model is unavailable for a pool member is dropped by _Sides; the opt-in guard refuses it."""
    import pytest
    from fuelblend.design import blendstock_v5 as B
    from fuelblend.design import engine_v4 as E
    monkeypatch.setenv(fu.REGISTRY_ADDITIONS_ENV, "1")
    fu._default_registry.cache_clear()
    try:
        bp = B.BlendPool([B.pure_component(n) for n in ("decalin_cis", "n-decane", "tetralin")])
        models = B.build_models(bp)
        with pytest.raises(E.ModelUnavailable, match="dropped from the optimiser"):
            B.design_pool("astm:jet_a", bp, models, None, n_starts=1, seed=0, strict_constraints=True)
    finally:
        monkeypatch.delenv(fu.REGISTRY_ADDITIONS_ENV)
        fu._default_registry.cache_clear()
