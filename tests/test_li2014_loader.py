"""Fixture test: 3 rows copied exactly as printed (Li et al. 2014 Tables 5, 6 and 10), plus layout checks."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fuelblend.data.loaders_li2014 import load_li2014, system_id_from_smiles  # noqa: E402


def test_three_rows_as_printed():
    d = load_li2014()
    li = d[d.source == "li2014"]
    r = li[(li.table == "Table 5")].iloc[0]       # octane 0.31, dodecane 0.69, exp 42.00
    assert r.measured_C == 42.0 and json.loads(r.x_mole) == [0.31, 0.69] and r.page == 3
    r = li[(li.table == "Table 6")].iloc[1]       # nonane 0.46, decane 0.54, exp 44.00
    assert r.measured_C == 44.0 and json.loads(r.x_mole) == [0.46, 0.54]
    r = li[(li.table == "Table 10")].iloc[0]      # nonane 0.19 decane 0.18 dodecane 0.63, exp 58.40
    assert r.measured_C == 58.4 and json.loads(r.x_mole) == [0.19, 0.18, 0.63] and r.page == 4 and r.doi == "10.1155/2014/239023"


def test_layout_and_ids():
    d = load_li2014()
    assert len(d) == 33 and d.n_comp.isin([2, 3]).all()
    assert all(abs(sum(json.loads(x)) - 1) < 1e-5 for x in d.x_mole)
    assert d.system_id.nunique() == 6
    # system_id is order independent
    assert system_id_from_smiles(["CCCCCCCC", "CCCCCCCCCCCC"]) == system_id_from_smiles(["CCCCCCCCCCCC", "CCCCCCCC"])
