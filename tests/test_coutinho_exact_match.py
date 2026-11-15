"""Decision 5 (2026-11-15): exact aromatic matching in the Coutinho table for v5 callers; legacy substring default kept for v4."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.thermo.coutinho_wax import get_pure_thermal_properties

ROOT = Path(__file__).resolve().parents[1]
MISASSIGNED = ROOT / "results" / "v5" / "BASE" / "delta_coutinho_lookup_audit" / "misassigned.csv"

ROWS = [("BENZ", "c1ccccc1", 6), ("PXYL", "Cc1ccc(C)cc1", 8), ("NAPH", "c1ccc2ccccc2c1", 10),
        ("TOL", "Cc1ccccc1", 7), ("FLUO", "c1ccc2c(c1)Cc1ccccc1-2", 13), ("PHEN", "c1ccc2c(c1)ccc1ccccc12", 14)]


@pytest.fixture
def df6():
    return pd.DataFrame({"smiles": [r[1] for r in ROWS], "cls": "aromatic", "n_carbon": [float(r[2]) for r in ROWS]},
                        index=[r[0] for r in ROWS])


def test_exact_keeps_specials_and_drops_others(df6):
    t = get_pure_thermal_properties(df6, aromatic_match="exact")
    assert t["BENZ"]["tm_K"] == 278.7 and t["PXYL"]["tm_K"] == 286.4 and t["NAPH"]["tm_K"] == 353.4
    assert t["BENZ"]["dh_fus_J_mol"] == 9870.0 and t["PXYL"]["dh_fus_J_mol"] == 17100.0 and t["NAPH"]["dh_fus_J_mol"] == 19000.0
    assert all(np.isnan(t[k]["tm_K"]) and np.isnan(t[k]["dh_fus_J_mol"]) for k in ("TOL", "FLUO", "PHEN"))


def test_exact_matches_noncanonical_smiles_and_durene():
    df = pd.DataFrame({"smiles": ["C1=CC=CC=C1", "Cc1cc(C)c(C)cc1C", "c1ccc2ccccc2c1.O"], "cls": "aromatic", "n_carbon": [6.0, 10.0, 10.0]},
                      index=["BZ", "DUR", "NAPHSALT"])
    t = get_pure_thermal_properties(df, aromatic_match="exact")
    assert t["BZ"]["tm_K"] == 278.7 and t["DUR"]["tm_K"] == 352.4 and t["DUR"]["dh_fus_J_mol"] == 21000.0
    assert np.isnan(t["NAPHSALT"]["tm_K"])          # a multi-fragment SMILES is not naphthalene


def test_substring_default_unchanged(df6):
    t = get_pure_thermal_properties(df6)
    assert t["TOL"]["tm_K"] == 278.7 and t["FLUO"]["tm_K"] == 278.7      # legacy defect kept for v4 reproduction
    assert t == get_pure_thermal_properties(df6, aromatic_match="substring")


def test_bad_mode_raises(df6):
    with pytest.raises(ValueError):
        get_pure_thermal_properties(df6, aromatic_match="fuzzy")


def test_non_aromatic_branches_unchanged():
    df = pd.DataFrame({"smiles": ["CCCCCCCCCC", "C1CCCCC1"], "cls": ["n_paraffin", "cycloparaffin"], "n_carbon": [10.0, 6.0]},
                      index=["C10", "CH"])
    assert get_pure_thermal_properties(df, aromatic_match="exact").keys() == get_pure_thermal_properties(df).keys()
    a, b = get_pure_thermal_properties(df, aromatic_match="exact"), get_pure_thermal_properties(df)
    for k in a:
        for f in a[k]:
            assert (a[k][f] == b[k][f]) or (a[k][f] != a[k][f] and b[k][f] != b[k][f])


def test_no_misassigned_under_exact():
    # every row of the 2026-11 lookup audit (112 aromatics given another molecule's constants) must come out NaN under "exact"
    m = pd.read_csv(MISASSIGNED)
    assert len(m) == 112
    df = pd.DataFrame({"smiles": m.smiles.values, "cls": "aromatic", "n_carbon": 10.0}, index=m.inchikey.values)
    t = get_pure_thermal_properties(df, aromatic_match="exact")
    kept = [ik for ik in m.inchikey if np.isfinite(t[ik]["tm_K"]) or np.isfinite(t[ik]["dh_fus_J_mol"])]
    assert kept == []
    # and the legacy mode still reproduces the audit's assigned value for each of them
    t0 = get_pure_thermal_properties(df)
    assert all(t0[r.inchikey]["tm_K"] == r.assigned_tm_K for r in m.itertuples())


def test_build_tmap_defaults_to_exact():
    import inspect
    from fuelblend.eval.composition_route import build_tmap
    assert inspect.signature(build_tmap).parameters["aromatic_match"].default == "exact"
    from fuelblend.thermo.coutinho_wax import get_pure_thermal_properties as g
    assert inspect.signature(g).parameters["aromatic_match"].default == "substring"


def test_predict_fuel_alkylbenzene_only_freezing_not_silently_certified():       # Review Focus 4
    from fuelblend.api.predict_route import predict_fuel
    # ethylbenzene has no verified constants and (under "exact") no Coutinho entry: no value, uncertified, never a benzene constant
    fp = predict_fuel({"CCc1ccccc1": 1.0})["freezing_point_C"]
    assert fp["value"] is None and fp["status"] == "non_finite" and fp["uncertified"] is True
    assert "freezing_uncertified" in fp["warnings"]
    assert fp["variants"]["coutinho_constants"]["value"] is None
    # toluene (verified) + ethylbenzene (none): a value from the verified component only, still uncertified (ethylbenzene has no constants)
    fp = predict_fuel({"Cc1ccccc1": 0.5, "CCc1ccccc1": 0.5})["freezing_point_C"]
    assert fp["status"] == "ok" and fp["uncertified"] is True and "freezing_uncertified" in fp["warnings"]
    assert any(w.startswith("components_without_sle_constants") for w in fp["warnings"])
    assert fp["variants"]["coutinho_constants"]["value"] is None            # no Coutinho fallback value for either aromatic
    # verified toluene alone stays certified (unchanged)
    fp = predict_fuel({"Cc1ccccc1": 1.0})["freezing_point_C"]
    assert fp["status"] == "ok" and fp["uncertified"] is False
