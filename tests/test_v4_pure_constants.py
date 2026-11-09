"""v4 Task 4: pure-component constant audit (F1).

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

DIFFS = ROOT / "results" / "v4" / "P1" / "pure_constant_diffs.csv"
PROPS = ROOT / "data" / "processed" / "v4" / "pure_props_v4.csv"


def test_known_reference_values():
    tm = A.reference_value("CCCCCCCCCCCC", "Tm")
    assert tm is not None
    assert abs(tm[0] - (-9.6)) <= 0.5
    assert isinstance(tm[1], str) and tm[1]
    # n-octane RON: listed -10 must be flagged against Co-Optima when the file exists
    ref = A.reference_value("CCCCCCCC", "RON")
    df = A.run(write=False)["diffs"]
    row = df[(df.smiles == "CCCCCCCC") & (df.prop == "RON") & (df.listed_value == -10.0)]
    assert len(row) >= 1
    if ref is not None:
        assert abs(ref[0] - (-10.0)) > 3
        assert (row.flag == "flagged").all()
    else:
        assert (row.flag == "no_reference").all()


def test_unknown_returns_none():
    assert A.reference_value("not_a_smiles", "Tm") is None
    assert A.reference_value("CCCC", "bogus_prop") is None


def test_threshold_logic():
    assert A.is_flagged("Tm", 3.0, None) is True
    assert A.is_flagged("Tb", 1.9, None) is False
    assert A.is_flagged("rho15", None, 0.011) is True
    assert A.is_flagged("rho15", None, 0.009) is False
    assert A.is_flagged("DCN", 3.1, None) is True
    assert A.is_flagged("Psat37_8", None, 0.06) is True


def test_outputs_written_and_consistent():
    A.run(write=True)
    d = pd.read_csv(DIFFS)
    assert list(d.columns[:9]) == ["compound", "smiles", "prop", "listed_value", "listed_source",
                                   "ref_value", "ref_source", "delta", "flag"]
    assert set(d.flag) <= {"ok", "flagged", "no_reference"}
    assert (d.flag == "flagged").sum() > 0
    p = pd.read_csv(PROPS)
    assert {"smiles", "prop", "value", "source"} <= set(p.columns)
    assert not p.duplicated(["smiles", "prop"]).any()
    assert p.source.notna().all()
    # a flagged row with a reference gets the reference value in v4
    f = d[d.flag == "flagged"].iloc[0]
    v = p[(p.smiles == f.smiles) & (p.prop == f["prop"])].iloc[0]
    assert v.value == pytest.approx(f.ref_value)
