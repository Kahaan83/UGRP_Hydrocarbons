"""Tests for fuelblend.data.loaders_thermoml2020.build (ThermoML JCED 2020 extractor output -> measurements_dedup schema).

The fixture mirrors results/v5/thermoml_2020/thermoml_hc_2020.csv: JSON lists of InChI strings / names / mole fractions,
`file` like 10.1021/acs.jced.0c00581.xml, native ThermoML property strings and units.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from fuelblend.data import loaders_thermoml2020 as L
from fuelblend.data.loaders_boehm2022 import COLS

ROOT = Path(__file__).resolve().parents[1]
REAL = ROOT / "results" / "v5" / "thermoml_2020" / "thermoml_hc_2020.csv"
REAL_LOG = ROOT / "data" / "real_fuels" / "v5" / "thermoml2020_build_log.json"

# real InChIs copied from the extractor output
INCHI_A = "InChI=1S/C12H26/c1-10(8-11(2,3)4)9-12(5,6)7/h10H,8-9H2,1-7H3"  # 2,2,4,6,6-pentamethylheptane
INCHI_B = "InChI=1S/C9H12/c1-2-6-9-7-4-3-5-8-9/h3-5,7-8H,2,6H2,1H3"       # propylbenzene
NAME_A, NAME_B = "2,2,4,6,6-pentamethylheptane", "propylbenzene"

VISC = "Viscosity, Pa*s"
KIN = "Kinematic viscosity, m2/s"
RHO = "Mass density, kg/m3"
F1 = "10.1021/acs.jced.0c00020.xml"
F2 = "10.1021/acs.jced.0c00581.xml"


def row(prop, value, T=298.15, x=(0.4, 0.6), f=F1, comp_type="mole fraction", pair=True, P=102.0):
    inchis, names = ([INCHI_A, INCHI_B], [NAME_A, NAME_B]) if pair else ([INCHI_A], [NAME_A])
    return dict(file=f, prop=prop, prop_phase="Liquid", value=value, T_K=T, P_kPa=P, n_comp=len(inchis),
                inchis=json.dumps(inchis), names=json.dumps(names), comp_type=comp_type, comp_values=json.dumps(list(x)))


def run(tmp_path, rows):
    p = tmp_path / "fixture.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    log = {}
    return L.build(p, log), log


def keys_and_mw():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    out = {}
    for inchi in (INCHI_A, INCHI_B):
        m = Chem.MolFromInchi(inchi)
        out[inchi] = (Chem.MolToInchiKey(m), Descriptors.MolWt(m))
    return out


def test_pure_endpoint_and_binary_ordering_and_mass_fraction(tmp_path):
    out, _ = run(tmp_path, [row(RHO, 700.0, x=(0.0, 1.0)), row(RHO, 800.0, x=(0.4, 0.6))])
    pure = out[out.kind == "pure"]
    assert len(pure) == 1 and pure.iloc[0].n_comp == 1
    assert len(json.loads(pure.iloc[0].comp_inchikeys)) == 1
    assert json.loads(pure.iloc[0].comp_inchikeys) == [keys_and_mw()[INCHI_B][0]]  # fraction 0 removed -> only B remains

    b = out[out.kind == "binary"].iloc[0]
    assert b.n_comp == 2
    kw = keys_and_mw()
    iks = json.loads(b.comp_inchikeys)
    assert iks == sorted(iks)
    # independent expected values, ordered by InChIKey
    comp = {kw[INCHI_A][0]: (0.4, kw[INCHI_A][1]), kw[INCHI_B][0]: (0.6, kw[INCHI_B][1])}
    exp_x = [comp[k][0] for k in iks]
    mass = [comp[k][0] * comp[k][1] for k in iks]
    exp_w = [m / sum(mass) for m in mass]
    assert json.loads(b.x_mole) == pytest.approx(exp_x, abs=1e-5)
    assert json.loads(b.w_mass) == pytest.approx(exp_w, abs=1e-5)


def test_unit_conversions(tmp_path):
    out, _ = run(tmp_path, [row(VISC, 0.001), row(KIN, 1e-6, T=303.15), row(RHO, 745.16, T=313.15)])
    g = lambda prop: out[out.property == prop].iloc[0].value
    assert g("viscosity_dyn_mPas") == pytest.approx(1.0)
    assert g("viscosity_kin_cSt") == pytest.approx(1.0)
    assert g("density_kg_m3") == pytest.approx(745.16)


def test_dynamic_from_kinematic_with_matching_density(tmp_path):
    rows = [row(KIN, 2e-6, T=298.15), row(RHO, 800.0, T=298.15)]
    out, log = run(tmp_path, rows)
    dyn = out[out.property == "viscosity_dyn_mPas"]
    assert len(dyn) == 1
    assert dyn.iloc[0].value == pytest.approx(2.0 * 800.0 / 1000.0)
    assert len(out[out.property == "viscosity_kin_cSt"]) == 1  # kinematic row kept as well
    assert log["dynamic_from_kinematic"] == 1


def test_no_dynamic_without_matching_density(tmp_path):
    rows = [row(KIN, 2e-6, T=298.15),
            row(RHO, 800.0, T=313.15),                 # different temperature
            row(RHO, 800.0, T=298.15, f=F2),           # different file
            row(RHO, 800.0, T=298.15, x=(0.5, 0.5))]   # different composition
    out, log = run(tmp_path, rows)
    assert (out.property == "viscosity_dyn_mPas").sum() == 0
    assert log["kinematic_rows"] == 1 and log["dynamic_from_kinematic"] == 0


def test_repeats_collapse_to_median_and_conflict_flag(tmp_path):
    out, _ = run(tmp_path, [row(RHO, 800.0), row(RHO, 850.0)])  # 50 > 5 % of median 825
    assert len(out) == 1
    r = out.iloc[0]
    assert r.n_obs == 2 and r.value == pytest.approx(825.0) and bool(r.conflict) is True

    out, _ = run(tmp_path, [row(RHO, 800.0), row(RHO, 801.0)])
    assert len(out) == 1
    r = out.iloc[0]
    assert r.n_obs == 2 and r.value == pytest.approx(800.5) and bool(r.conflict) is False


def test_bad_composition_and_unsupported_scale_dropped(tmp_path):
    rows = [row(RHO, 800.0),
            row(RHO, 810.0, T=303.15, x=(0.4, 0.5)),                         # sums to 0.9
            row(RHO, 820.0, T=308.15, comp_type="mass fraction")]
    out, log = run(tmp_path, rows)
    assert len(out) == 1 and out.iloc[0].value == pytest.approx(800.0)
    assert log["dropped"]["bad composition"] == 1
    assert log["dropped"]["unsupported composition scale"] == 1
    assert log["rows_read"] == 3


def test_unsupported_property_dropped(tmp_path):
    out, log = run(tmp_path, [row(RHO, 800.0), row("Heat capacity at constant pressure, J/K/mol", 200.0, T=303.15)])
    assert len(out) == 1
    assert log["dropped"]["unsupported property"] == 1


def test_columns_and_source(tmp_path):
    out, _ = run(tmp_path, [row(RHO, 800.0), row(VISC, 0.001)])
    assert list(out.columns) == list(COLS)
    assert (out.sources == "thermoml_jced_2020_unlisted").all()


def test_real_file_matches_committed_build_log():
    if not REAL.exists():
        pytest.skip("real ThermoML 2020 extractor output not present")
    exp = json.loads(REAL_LOG.read_text(encoding="utf-8"))
    log = {}
    out = L.build(REAL, log)
    assert len(out) == exp["rows_out"] == 3399
    assert log["mixture_viscosity_systems"] == exp["mixture_viscosity_systems"] == 19
    assert log["conflicts"] == exp["conflicts"] == 0
    assert int(out.conflict.sum()) == 0
