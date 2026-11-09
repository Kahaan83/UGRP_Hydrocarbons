"""v4 Task 4: pure-component constant audit (F1), policy R16 (never auto-substitute).

reference_value(smiles, prop) -> (value, source) | None, plus the audit of the listed tables
(purchasable_compounds.csv, CANDIDATE_REGISTRY) against chemicals/thermo and the Co-Optima/NREL files.
"""
from pathlib import Path
import importlib.util
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
pytest.importorskip("chemicals")
pytest.importorskip("rdkit")
_spec = importlib.util.spec_from_file_location("v4_apc", ROOT / "scripts" / "v4" / "audit_pure_constants.py")
A = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A)


@pytest.fixture(scope="module")
def res():
    return A.run(write=False)


def test_known_reference_values(res):
    tm = A.reference_value("CCCCCCCCCCCC", "Tm")
    assert tm is not None
    assert abs(tm[0] - (-9.6)) <= 0.5
    assert isinstance(tm[1], str) and tm[1]
    # the Co-Optima reference file must be on disk; n-octane RON listed -10 must be flagged
    assert A.COOPTIMA.exists(), f"missing {A.COOPTIMA}"
    ref = A.reference_value("CCCCCCCC", "RON")
    assert ref is not None and abs(ref[0] - (-10.0)) > 3
    df = res["diffs"]
    row = df[(df.smiles == "CCCCCCCC") & (df.prop == "RON") & (df.listed_value == -10.0)]
    assert len(row) >= 1 and (row.flag == "flagged").all()


def test_unknown_returns_none():
    assert A.reference_value("not_a_smiles", "Tm") is None
    assert A.reference_value("CCCC", "bogus_prop") is None


def test_threshold_logic():
    assert A.is_flagged("Tm", 3.0, None) is True
    assert A.is_flagged("Tb", 1.9, None) is False
    assert A.is_flagged("rho15", None, 0.011) is True
    assert A.is_flagged("rho15", None, 0.009) is False
    assert A.is_flagged("DCN", 3.1, None) is True
    assert A.is_flagged("Psat37_8", None, 0.06, 0.5) is True
    # relative Psat is noise below 0.01 kPa
    assert A.is_flagged("Psat37_8", None, 3.39, 0.001) is False


def test_decalin_trans_density_stays_listed_and_disputed(res):
    p = res["props"]
    r = p[(p.compound == "decalin_trans") & (p.prop == "rho15")]
    assert len(r) == 1
    r = r.iloc[0]
    assert r.value == pytest.approx(869.5) and r.listed_value == pytest.approx(869.5)
    assert r.status == "disputed"
    assert r.ref_value == pytest.approx(846.9, abs=0.5)


def test_heptamethylnonane_tm_disputed_not_hidden(res):
    p = res["props"]
    r = p[(p.compound == "heptamethylnonane") & (p.prop == "Tm")].iloc[0]
    assert r.kind == "pure" and r.status == "disputed"
    assert r.value == pytest.approx(-65.0)
    assert not (res["props"][(res["props"].prop == "Psat37_8") & (res["props"].status == "reference_only")]
                .ref_value < A.PSAT_MIN_REF_KPA).any()


def test_unknown_prop_lookup_is_error_not_dcn():
    assert A._lookup("CCCC", "bogus")[2] == "no_reference"  # unlisted prop never reaches the DCN branch


def test_never_substituted(res):
    p = res["props"]
    listed = p[p.status != "reference_only"]
    assert (listed.value == listed.listed_value).all()
    ro = p[p.status == "reference_only"]
    assert (ro.kind == "pure").all() and ro.listed_value.isna().all()
    assert set(p.status) <= {"verified", "disputed", "unverified", "not_comparable", "lookup_error",
                             "reference_only"}
    assert not p.duplicated(["smiles", "prop"]).any()


def test_mixture_rows_not_comparable(res):
    d, p = res["diffs"], res["props"]
    for name in ("isododecane_commercial", "heptamethylnonane_commercial", "decalin_mixed", "decalin"):
        rows = d[d.compound == name]
        assert len(rows) > 0
        assert (rows.kind == "mixture_or_surrogate").all()
        assert (rows.flag == "not_comparable").all()
        assert (p[p.compound == name].status == "not_comparable").all()
    assert A.is_mixture("foo", "CCCCCCCCCCCC", 200.0)  # MW 170.3 vs 200: > 3 %
    assert not A.is_mixture("n-dodecane", "CCCCCCCCCCCC", 170.33)


def test_lookup_error_distinct_from_no_reference(monkeypatch):
    import chemicals
    # no data: a compound that search_chemical cannot resolve
    assert A._lookup("not_a_smiles", "Tm")[2] == "no_reference"

    def boom(cas, *a, **k):
        raise KeyError(cas)

    monkeypatch.setattr(chemicals, "Tm", boom)
    val, src, state = A._lookup("CCCCCCCCCCCC", "Tm")
    assert val is None and state == "lookup_error:KeyError" and src == "lookup_error:KeyError"
    assert A.reference_value("CCCCCCCCCCCC", "Tm") is None


def test_outputs_written_and_consistent(tmp_path):
    d_path, p_path = tmp_path / "diffs.csv", tmp_path / "props.csv"
    A.run(write=True, out_diffs=d_path, out_props=p_path)
    d = pd.read_csv(d_path)
    assert list(d.columns[:9]) == ["compound", "smiles", "prop", "listed_value", "listed_source",
                                   "ref_value", "ref_source", "delta", "flag"]
    assert "kind" in d.columns
    assert set(d.flag) <= {"ok", "flagged", "no_reference", "not_comparable", "not_assessed", "lookup_error"}
    assert (d.flag == "flagged").sum() > 0
    p = pd.read_csv(p_path)
    assert {"smiles", "prop", "value", "source", "status", "listed_value", "listed_source", "ref_value",
            "ref_source", "delta", "listed_value_alt", "kind"} <= set(p.columns)
    assert p.source.notna().all()
