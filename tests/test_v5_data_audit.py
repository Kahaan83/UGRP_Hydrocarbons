"""Rules used by the v5 data-audit scripts (scripts/v5/eval_sle_extended_v5.py)."""
import sys
from pathlib import Path

import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "v5"))
sle = pytest.importorskip("eval_sle_extended_v5")


def test_ajftd_gcxgc_metric_labels():
    from fuelblend.data.loaders_ajftd_gcxgc import parse_metric
    assert parse_metric("n-C10") == ("n_paraffin", 10, False)
    assert parse_metric("n-C07 and lower") == ("n_paraffin", 7, True)
    assert parse_metric("C3-benzene (C09)") == ("alkylbenzene", 9, False)          # carbon number is the bracketed one
    assert parse_metric("C10+-benzene (C16+)") == ("alkylbenzene", 16, True)
    assert parse_metric("C07-monocyclocycloparaffins") == ("monocycloparaffin", 7, False)
    assert parse_metric("C17+-dicycloparaffins") == ("dicycloparaffin", 17, True)
    assert parse_metric("C1-naphthalene (C11)") == ("diaromatic", 11, False)
    assert parse_metric("C07 and lower-iso") == ("iso_paraffin", 7, True)
    assert parse_metric("aromatics") is None and parse_metric("n-paraffins") is None   # class totals are not bins


def test_ajftd_gcxgc_load_mirrors_detail_file_layout(tmp_path):
    """Fixture has the real long layout: one cell per line, table 0 = composition tree, totals mixed with leaf rows."""
    import pandas as pd
    from fuelblend.data.loaders_ajftd_gcxgc import fuel_table, load
    cells = [(1, "aromatics", "20.0"), (2, "Metric Weight % toluene (C07) 20.0", None), (3, "toluene (C07)", "20.0"),
             (4, "n-paraffins", "30"), (5, "n-C10", "30.0"), (6, "C11-isoparaffins", "49.5"), (7, "average_molecular_formula", None)]
    rows = []
    for r, metric, wt in cells:
        rows.append((7, "99999", 0, r, "Metric", metric))
        rows.append((7, "99999", 0, r, "Weight %", wt))
    rows.append((7, "99999", 3, 1, "Metric", "n-C12"))          # another table of the page must be ignored
    rows.append((7, "99999", 3, 1, "Weight %", "55"))
    p = tmp_path / "details.csv"
    pd.DataFrame(rows, columns=["fuel_id", "fuel_label", "table_index", "row", "column", "value"]).to_csv(p, index=False)
    comp, unparsed = load(p)
    assert sorted(comp.hc_class) == ["alkylbenzene", "iso_paraffin", "n_paraffin"]
    assert comp.mass_pct.sum() == pytest.approx(99.5)
    assert "aromatics" in unparsed and "n-paraffins" in unparsed
    f = fuel_table(comp, wide=tmp_path / "absent.csv")
    assert bool(f.composition_complete.iloc[0]) and f.n_bins.iloc[0] == 3


def test_gcxgc_bin_mapping_counts_unmapped_mass():
    """Bins without a representative (tricycloparaffins, dienes, or no compound within the carbon gap) are unmapped mass."""
    import pandas as pd
    g1 = pytest.importorskip("eval_ajftd_gcxgc_v5")

    class FakeReps:
        table = {("n_paraffin", 10): 0, ("n_paraffin", 12): 1, ("cyclo_di", 10): 2}

        def get(self, kind, n, max_gap=3):
            ns = [k[1] for k in self.table if k[0] == kind]
            if not ns:
                return None
            b = min(ns, key=lambda k: (abs(k - n), k))
            return self.table[(kind, b)] if abs(b - n) <= max_gap else None

    bins = pd.DataFrame([("n_paraffin", 10, 40.0), ("n_paraffin", 11, 10.0), ("n_paraffin", 20, 5.0),
                         ("dicycloparaffin", 16, 25.0), ("tricycloparaffin", 11, 15.0), ("diene", 12, 5.0)],
                        columns=["hc_class", "n_carbon", "mass_pct"])
    agg, unmapped, lost = g1.map_fuel(bins, FakeReps())
    assert agg == {0: 50.0, 2: 25.0}                     # C11 goes to the nearest (C10); decalin-type uses any carbon number
    assert unmapped == pytest.approx(0.25)
    assert lost == pytest.approx({"no_rep_n_paraffin": 0.05, "no_rep_tricycloparaffin": 0.15, "no_rep_diene": 0.05})


def test_tm_consensus_agreeing_sources_give_median():
    tm, rule = sle.tm_consensus({"CRC_ORG": 263.6, "WEBBOOK": 263.5, "OPEN_NTBKM": 263.55})
    assert rule == "agree" and tm == pytest.approx(263.55)


def test_tm_consensus_drops_one_outlier_and_names_it():
    # 3-methylpentane: one compilation lists the glass transition (110 K) instead of the melting point (155 K)
    tm, rule = sle.tm_consensus({"CRC_ORG": 110.26, "WEBBOOK": 155.15, "OPEN_NTBKM": 155.15, "COMMON_CHEMISTRY": 155.15})
    assert tm == pytest.approx(155.15) and rule.startswith("agree_dropped_CRC_ORG")


def test_tm_consensus_two_camps_is_disputed():
    # 1-methylnaphthalene: sources split between about 243 K and 251 K, so no value is accepted
    tm, rule = sle.tm_consensus({"CRC_ORG": 242.72, "WEBBOOK": 244.0, "OPEN_NTBKM": 251.15, "COMMON_CHEMISTRY": 251.15})
    assert tm is None and rule == "disputed"


def test_tm_consensus_needs_two_sources():
    assert sle.tm_consensus({"CRC_ORG": 300.0}) == (None, "single_source")


def test_n_alkane_detection_and_unassigned_stereo():
    assert sle.n_alkane_cn(Chem.MolFromSmiles("CCCCCCCCCCCC")) == 12
    assert sle.n_alkane_cn(Chem.MolFromSmiles("CCC(C)CC")) is None
    assert sle.n_alkane_cn(Chem.MolFromSmiles("C1CCCCC1")) is None
    assert sle.unassigned_stereo(Chem.MolFromSmiles("C1CCC2CCCCC2C1"))          # decalin, ring fusion not specified
    assert not sle.unassigned_stereo(Chem.MolFromSmiles("C1CC[C@H]2CCCC[C@@H]2C1"))
