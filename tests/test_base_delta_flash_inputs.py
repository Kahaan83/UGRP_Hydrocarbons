"""v5 BASE Task 10: flash input package delta."""
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "v5" / "BASE" / "delta_flash_inputs"
REG = [ROOT / "data/processed/v5/pure_props_v5.csv", ROOT / "data/processed/v5/lfl_v5.csv"]


def _sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def ran():
    before = [_sha(p) for p in REG]
    subprocess.run([sys.executable, str(ROOT / "scripts/v5/base/delta_flash_inputs.py")], cwd=ROOT, check=True)
    return before


def test_v4_designated_reproduced(ran):
    d = pd.read_csv(D / "delta.csv")
    s = pd.read_csv(ROOT / "results/v5/BASE/scoreboard.csv")
    crc = s[(s.tier == "CRC647_r3") & (s.property == "flash_point_C") & (s.variant == "designated") & (s.label == "dev")].iloc[0]
    assert crc.bias == pytest.approx(0.44, abs=0.01)                                  # frozen scoreboard value, no overrides
    lab = d[(d.tier == "T1_lab") & (d.variant == "designated")].iloc[0]
    assert lab.margin90 == pytest.approx(3.698768, abs=1e-4)                          # v4 margin reproduced from re-predicted rows
    assert lab.bias == pytest.approx(-3.3007, abs=1e-3)


def test_overrides_do_not_touch_registry_files(ran):
    assert [_sha(p) for p in REG] == ran


def test_overlay_unsafe_rate_never_rises():
    d = pd.read_csv(D / "delta.csv")
    o = d[(d.variant == "lfl_overlay") & d.baseline_unsafe_rate.notna()]
    assert len(o) >= 6
    assert (o.unsafe_rate <= o.baseline_unsafe_rate + 1e-9).all()


def test_margin_table_has_leave_mesitylene_out():
    m = pd.read_csv(D / "lab_margin.csv")
    assert set(m.variant) == {"designated", "lfl_overlay", "mesitylene_fixed", "lfl_overlay+mesitylene_fixed"}
    assert set(m.scope) == {"all_31_systems", "leave_mesitylene_system_out"}
    assert np.allclose(m.v4_margin, 3.698768, atol=1e-4)
