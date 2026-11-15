"""U7 (CLM-67): the optimiser smoke stand-in (``alt_standin``) is on by default in ``design_pool``; passing False gives the old behaviour."""
import inspect

import pytest

from fuelblend.design import blendstock_v5 as B
from fuelblend.design import engine_v4 as E


def test_alt_standin_default_is_true():
    assert inspect.signature(B.design_pool).parameters["alt_standin"].default is True


class _Stop(Exception):
    pass


@pytest.mark.parametrize("passed, expected", [({}, True), ({"alt_standin": True}, True), ({"alt_standin": False}, False)])
def test_alt_standin_reaches_sides(monkeypatch, passed, expected):
    """design_pool forwards the flag to engine_v4._Sides; False restores the previous (stand-in off) behaviour."""
    seen = {}

    def fake_sides(*a, **kw):
        seen.update(kw)
        raise _Stop

    monkeypatch.setattr(E, "_Sides", fake_sides)
    bp = B.BlendPool([B.pure_component(n) for n in ("n-decane", "n-dodecane", "tetralin")])
    models = B.build_models(bp)
    with pytest.raises(_Stop):
        B.design_pool("astm:jet_a", bp, models, None, n_starts=1, seed=0, **passed)
    assert seen["alt_point_standin"] is expected


def test_engine_default_unchanged():
    """The engine's own default stays off, so direct _Sides callers (judge, tests) are unaffected."""
    assert inspect.signature(E._Sides.__init__).parameters["alt_point_standin"].default is False
