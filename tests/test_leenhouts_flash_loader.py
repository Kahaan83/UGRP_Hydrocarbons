"""v6 flash arm loader test. Fixture tests/data/leenhouts_flash_mix_sample.csv = 6 rows copied verbatim from the real
flashpoint_mix_exp.csv (HC ternary, two rows of one oxygenate binary, oxygenate ternary and quaternary, HC binary)."""
import json
from pathlib import Path

import pandas as pd
import pytest

from fuelblend.data import loaders_leenhouts_flash as L

FIX = Path(__file__).parent / "data" / "leenhouts_flash_mix_sample.csv"


@pytest.fixture(scope="module")
def built():
    log = {}
    return L.build(FIX, log), log


def test_row_and_class_counts(built):
    d, log = built
    assert len(d) == 6 and log["rows_dropped"] == 0
    assert d.chem_class.value_counts().to_dict() == {"oxy": 4, "hc": 2}
    assert d.is_hc.sum() == 2
    assert set(d.n_comp) == {2, 3, 4}


def test_units_and_schema(built):
    d, _ = built
    assert list(d.columns) == L.COLS
    assert (d.property == "flash_point_C").all() and (d.P_kPa == 101.325).all() and d.T_K.isna().all()
    assert d.value.between(-30, 130).all()  # degC, not K
    assert d.loc[d.fuel_id == "leenhouts_flash_mix_0", "value"].iloc[0] == pytest.approx(71.85)


def test_fractions_normalised_and_sorted(built):
    d, _ = built
    for r in d.itertuples():
        x, w, iks = json.loads(r.x_mole), json.loads(r.w_mass), json.loads(r.comp_inchikeys)
        assert sum(x) == pytest.approx(1, abs=1e-5) and sum(w) == pytest.approx(1, abs=1e-5)
        assert iks == sorted(iks) and len(iks) == len(x) == len(w) == r.n_comp


def test_same_system_shares_id(built):
    d, _ = built
    # fixture rows 1 and 2 are source rows 65 and 66 (same oxygenate binary, x = 0.01 and 0.02)
    two = d[d.fuel_id.isin(["leenhouts_flash_mix_1", "leenhouts_flash_mix_2"])]
    assert two.system_id.nunique() == 1 and len(two) == 2
    assert d.system_id.nunique() == 5


def test_dedup_is_order_invariant_and_flags_conflict(tmp_path):
    raw = pd.read_csv(FIX)
    r = raw[raw.index == 0].iloc[0].copy()
    # same blend, components listed in a different order, value 3 degC higher -> one row, median, conflict flagged
    r2 = r.copy()
    for k, j in ((0, 2), (2, 0)):
        r2[f"component_inchi_{k}"], r2[f"fraction_{k}"] = r[f"component_inchi_{j}"], r[f"fraction_{j}"]
    r2["value"] = r["value"] + 3.0
    p = tmp_path / "dup.csv"
    pd.DataFrame([r, r2]).to_csv(p, index=False)
    d = L.build(p)
    assert len(d) == 1 and d.n_obs.iloc[0] == 2 and bool(d.conflict.iloc[0])
    assert d.value.iloc[0] == pytest.approx(r["value"] + 1.5)


def test_native_unit_string_is_asserted(tmp_path):
    raw = pd.read_csv(FIX)
    raw["unit"] = "Kelvin"
    p = tmp_path / "bad.csv"
    raw.to_csv(p, index=False)
    with pytest.raises(AssertionError):
        L.build(p)
