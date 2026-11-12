"""v5 BASE Task 12: better representatives."""
import json
import subprocess
import sys
import warnings
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
D = ROOT / "results" / "v5" / "BASE" / "delta_representatives"


@pytest.fixture(scope="module")
def env():
    warnings.simplefilter("ignore")
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SUBKINDS, SubReps
    from fuelblend.data.representatives_v5 import SubRepsV5, candidate_table
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    vp = lambda s: _lookup(s, None)                                                          # noqa: E731
    tab = candidate_table(reg, prov)
    return dict(tab=tab, r0=SubReps(reg, prov, vp_lookup=vp), r1=SubRepsV5(reg, prov, vp_lookup=vp, table=tab), kinds=SUBKINDS)


def test_get_unchanged_outside_replaced_bins(env):
    from fuelblend.data.representatives_v5 import REPLACED_BINS
    for k in env["kinds"]:
        for n in range(4, 30):
            if (k, n) not in REPLACED_BINS:
                assert env["r1"].get(k, n) == env["r0"].get(k, n)
    assert env["r1"].get("alkylbenzene", 9) != env["r0"].get("alkylbenzene", 9)


def test_candidate_without_verified_constant_not_admitted(env):
    t = env["tab"]
    bad = t[~(t.has_antoine & t.has_lfl & t.has_viscosity & t.has_density)]
    assert len(bad) > 0 and not bad.admitted.any()
    assert not bool(t.loc[t.name == "1-methylnaphthalene", "admitted"].iloc[0])       # no Antoine in the registry


def test_delta_outputs():
    subprocess.run([sys.executable, str(ROOT / "scripts/v5/base/delta_representatives.py")], cwd=ROOT, check=True)
    d = pd.read_csv(D / "delta.csv")
    s = pd.read_csv(ROOT / "results/v5/BASE/scoreboard.csv")
    b = s[(s.property == "flash_point_C") & (s.tier == "T2_measured") & (s.set == "pooled") & (s.model == "ideal_le_chatelier") & (s.variant == "designated")].iloc[0]
    r = d[(d.property == "flash_point_C") & (d.set == "pooled") & (d.variant == "designated+subreps_v5")].iloc[0]
    assert r.baseline_bias == pytest.approx(b.bias, abs=1e-6)                       # frozen designated number reproduced
    assert set(json.loads((D / "readings.json").read_text())) >= {"reading_pass", "diesel_below_10pct_v5"}
    assert (D / "candidate_gaps.csv").exists()
