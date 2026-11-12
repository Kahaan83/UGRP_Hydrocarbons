import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
_s = importlib.util.spec_from_file_location("delta_lab_flash_li2014", ROOT / "scripts/v5/base/delta_lab_flash_li2014.py")
m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(m)


def _row(smi, x, f, sid="s"):
    return pd.Series(dict(system_id=sid, comp_smiles=json.dumps(smi), x_mole=json.dumps(x), measured_C=f))


def test_is_duplicate_aligns_components_and_tolerances():
    old = pd.DataFrame([dict(system_id="s", comp_smiles=json.dumps(["A", "B"]), x_mole=json.dumps([0.3, 0.7]), measured_C=40.0)])
    assert m.is_duplicate(_row(["B", "A"], [0.7, 0.3], 40.4), old)          # component order reversed, within tolerance
    assert not m.is_duplicate(_row(["A", "B"], [0.31, 0.69], 40.0), old)    # dx 0.01 > 0.005
    assert not m.is_duplicate(_row(["A", "B"], [0.3, 0.7], 41.0), old)      # dflash 1.0 > 0.5
    assert not m.is_duplicate(_row(["A", "B"], [0.3, 0.7], 40.0, sid="t"), old)


def test_outputs_consistent():
    out = ROOT / "results/v5/BASE/delta_lab_flash_li2014"
    r = json.loads((out / "readings.json").read_text())
    pts = pd.read_csv(out / "points.csv")
    assert r["n_new"] == int((pts.label == "dev").sum()) and not (pts.label == "held_out").any()
    mt = pd.read_csv(out / "lab_margin.csv").set_index(["variant", "scope"]).margin
    assert abs(mt[("lfl_overlay", "old_31_systems")] - 3.698768343688585) < 1e-9     # reproduces the v4 margin
    assert abs(mt[("lfl_overlay", "old_minus_mesitylene")] - 0.3106446441712771) < 1e-9
