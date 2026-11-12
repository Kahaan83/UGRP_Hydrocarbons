"""v5 BASE Task 16a: vapour-pressure delta uses the frozen designated flash numbers."""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "v5" / "BASE" / "delta_vapor_pressure"


def test_flash_side_equals_frozen_scoreboard():
    subprocess.run([sys.executable, str(ROOT / "scripts/v5/base/delta_vapor_pressure.py")], cwd=ROOT, check=True)
    r = json.loads((D / "readings.json").read_text())
    s = pd.read_csv(ROOT / "results/v5/BASE/scoreboard.csv")
    b = s[(s.property == "flash_point_C") & (s.tier == "T2_measured") & (s.set == "av3322") & (s.model == "ideal_le_chatelier") & (s.variant == "designated")].iloc[0]
    assert r["flash_bias_av3322"] == pytest.approx(b.bias, abs=1e-6)
    v = pd.read_csv(D / "rows_vp.csv")
    assert v[v.status == "ok"].pred.notna().all() and (v[v.status != "ok"].err_ln.isna()).all()
