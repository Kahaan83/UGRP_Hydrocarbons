"""v5 BASE Task 14: freezing follow-ups."""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "v5" / "BASE" / "delta_freezing"


def test_run_and_f1x_reproduces_s5():
    subprocess.run([sys.executable, str(ROOT / "scripts/v5/base/delta_freezing.py")], cwd=ROOT, check=True)
    r = json.loads((D / "readings.json").read_text())
    assert abs(r["f1x_pooled_mae"] - 3.07) <= 0.02 and r["f1x_pooled_n"] == 49


def test_wax_rows_never_in_liquidus_cell():
    d = pd.read_csv(D / "delta.csv")
    wax = d[d.evidence_path.fillna("").str.contains("wax_dauphin")]
    assert len(wax) > 0 and wax.set.str.contains("wax_appearance").all()
    assert not d[~d.set.str.contains("wax_appearance")].evidence_path.fillna("").str.contains("wax").any()
    s = pd.read_csv(ROOT / "results/v5/BASE/scoreboard.csv")
    assert not s.evidence_path.fillna("").str.contains("wax_dauphin").any()      # nothing wax-derived in the frozen liquidus cells
