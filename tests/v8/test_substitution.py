import numpy as np
import pytest

from fuelblend.v8.substitution import score_substitution


@pytest.mark.parametrize("prop", ["density", "freezing", "flash"])
def test_identity_when_ml_equals_measured(prop):
    df = score_substitution(prop, set_name="dev", _ml_override="measured")
    ok = df[df.status == "ok"]
    assert len(ok) > 0
    assert np.allclose(ok["mae_ml"], ok["mae_measured"])


def test_only_systems_with_measured_truth_scored():
    df = score_substitution("freezing", set_name="dev")
    ok = df[df.status == "ok"]
    assert (ok["n_substituted"] >= 1).all() and ok["mae_ml"].notna().all()


def test_flash_skips_components_without_measured_antoine():
    df = score_substitution("flash", set_name="dev")
    skipped = df[df["status"] == "skipped_no_measured_antoine"]
    assert skipped["mae_ml"].isna().all()                          # never scored with a substitute Antoine
    assert (df.loc[df["status"] == "ok", "n_substituted"] >= 1).all()
