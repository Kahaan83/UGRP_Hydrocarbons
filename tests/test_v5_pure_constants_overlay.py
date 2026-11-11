"""v5 pure-constant overlay: measured-density reference and primary listings (scripts/v5/audit_pure_constants_v5.py)."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("audit_v5", ROOT / "scripts/v5/audit_pure_constants_v5.py")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)


def _meas(rows):
    return pd.DataFrame([dict(comp_inchikeys='["K"]', T_K=t, T_key=t, P_kPa=101.325, value=v, ref=r, conflict=c)
                         for t, v, r, c in rows])


def test_measured_rho15_linear_fit_recovers_value():
    m = _meas([(283.15, 803.75, "a", False), (293.15, 796.25, "b", False), (303.15, 788.75, "b", False)])
    val, how = A.measured_rho15(m, "K")
    assert val == pytest.approx(800.0, abs=1e-6) and "interpolated" in how


@pytest.mark.parametrize("rows", [
    [(293.15, 796.0, "a", False), (303.15, 789.0, "b", False)],                               # only 2 temperatures
    [(283.15, 803.0, "a", False), (293.15, 796.0, "a", False), (303.15, 789.0, "a", False)],  # one reference
    [(298.15, 792.0, "a", False), (303.15, 789.0, "b", False), (308.15, 785.0, "b", False)],  # 10 K extrapolation
    [(283.15, 780.0, "a", False), (293.15, 796.0, "b", False), (303.15, 800.0, "b", False)],  # density rising with T
    [(283.15, 803.0, "a", True), (293.15, 796.0, "b", False), (303.15, 789.0, "b", False)],   # conflict row dropped
])
def test_measured_rho15_support_rule_rejects(rows):
    assert A.measured_rho15(_meas(rows), "K") is None


def test_other_stereoisomer_is_not_used():
    m = _meas([(283.15, 803.75, "a", False), (293.15, 796.25, "b", False), (303.15, 788.75, "b", False)])
    assert A.measured_rho15(m, "OTHER") is None


def test_overlay_on_real_tables():
    pp, ch = A.build()
    v4 = pd.read_csv(A.V4)
    assert list(pp.columns) == list(v4.columns) and len(pp) == len(v4)
    # R16: a measured reference never changes the listed value
    r = ch[ch.change == "rho15 reference -> measured"]
    assert len(r) and (pp.set_index(["compound", "prop"]).loc[list(zip(r.compound, r.prop)), "value"].to_numpy()
                       == v4.set_index(["compound", "prop"]).loc[list(zip(r.compound, r.prop)), "value"].to_numpy()).all()
    t = pp.set_index(["compound", "prop"])
    assert t.at[("decalin_trans", "rho15"), "status"] == "verified"
    assert abs(t.at[("decalin_trans", "rho15"), "delta"]) < 0.01
    assert t.at[("decalin_trans", "dHfus"), "status"] == "verified"
    assert t.at[("decalin_trans", "dHfus"), "value"] == pytest.approx(14.414)
    # untouched rows are identical to v4
    touched = set(zip(ch.compound, ch.prop))
    keep = [k not in touched for k in zip(pp.compound, pp.prop)]
    pd.testing.assert_frame_equal(pp[keep].reset_index(drop=True), v4[keep].reset_index(drop=True))
