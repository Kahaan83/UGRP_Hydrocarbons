"""v5 BASE Task 16b: benchmark-split summary uses the frozen lab flash errors."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "v5" / "BASE" / "delta_leenhouts_splits"
Z = ROOT / "extras/fuelblend_data/data/raw/zenodo/leenhouts_benchmark/splits.zip"


@pytest.mark.skipif(not Z.exists(), reason="extras/ absent")
def test_all_rows_matched_and_designated_mae_equals_frozen():
    subprocess.run([sys.executable, str(ROOT / "scripts/v5/base/delta_leenhouts_splits.py")], cwd=ROOT, check=True)
    r = json.loads((D / "readings.json").read_text())
    assert r["n_matched"] == 126
    assert r["all_rows_mae"] == pytest.approx(3.7444853, abs=1e-6)          # frozen T1 lab flash MAE (scoreboard designated row)
