"""SAK blendstock (Virent SAK, Feldhausen et al. 2022): source arithmetic, mapping, and identity of default pools when SAK is not selected."""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.design import blendstock_v5 as B

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parents[1]
TABLE1 = ROOT / "data" / "raw" / "sak" / "dot_64804" / "feldhausen2022_table1.csv"
BINS = ROOT / "data" / "real_fuels" / "v5" / "sak_feldhausen2022_bins.csv"
NAME = "sak_virent_feldhausen2022"
DEFAULT = ["jetA_posf10325", "hefa_came", "atj_posf11498"]


def test_bins_are_the_published_group_values():
    b = pd.read_csv(BINS)
    assert abs(b.mass_pct.sum() - 100.0) < 1e-9
    ab = b[b.hc_class == "alkylbenzene"].set_index("n_carbon").mass_pct
    assert dict(ab) == {8: 3.7, 9: 62.0, 10: 27.5}                         # Feldhausen 2022 sec 3.1
    ca = b[b.hc_class == "indan_tetralin"].set_index("n_carbon").mass_pct
    assert dict(ca) == {9: 0.8, 10: 1.2}
    assert abs(b[b.hc_class.str.startswith("no_rep")].mass_pct.sum() - 4.8) < 1e-9


def test_isomer_table_reproduces_paper_totals():
    t = pd.read_csv(TABLE1)
    assert len(t) == 26 and abs(t.mass_pct.sum() - 93.581) < 1e-9            # paper: 26 isomers, 93.6 %m
    g = t.groupby("n_carbon").mass_pct.sum()
    assert abs(g[8] - 3.683) < 1e-9 and abs(g[10] - 27.068) < 1e-9
    assert abs(t[t.group_type == "C9 alkylbenzene"].mass_pct.sum() - 62.031) < 1e-9
    f = t[t.peak_label == "f*"].mass_pct.sum()
    assert abs(f - 31.631) < 1e-9                                            # text: 31.6 %m


def test_sak_blendstock_maps_and_is_optional():
    b = B.load_blendstock(NAME)
    assert abs(b.w.sum() - 1.0) < 1e-12 and (b.w > 0).all()
    assert abs(b.unmapped_fraction - 0.048) < 1e-9
    assert b.synthetic and b.represented
    assert set(b.class_fractions()) == {"aromatic"} and abs(b.class_fractions()["aromatic"] - 1.0) < 1e-12
    assert NAME not in B.BLENDSTOCKS and NAME in B.OPTIONAL_BLENDSTOCKS
    with pytest.raises(KeyError):
        B.load_blendstock("not_a_blendstock")


def test_default_blendstocks_and_pools_unchanged_when_sak_not_selected():
    assert list(B.BLENDSTOCKS) == DEFAULT
    before = B.BlendPool([B.load_blendstock(n) for n in DEFAULT])
    B.load_blendstock(NAME)                                                   # loading SAK must not perturb anything
    after = B.BlendPool([B.load_blendstock(n) for n in DEFAULT])
    assert before.comp_ids == after.comp_ids and before.names == after.names
    assert np.array_equal(before.A, after.A)
    u = np.array([0.5, 0.3, 0.2])
    assert np.array_equal(before.expand(u), after.expand(u))
    assert before.synthetic_mass(u) == after.synthetic_mass(u)


def test_sak_can_join_a_pool_without_touching_other_columns():
    base = B.BlendPool([B.load_blendstock(n) for n in DEFAULT])
    full = B.BlendPool([B.load_blendstock(n) for n in DEFAULT + [NAME]])
    ids = {c: i for i, c in enumerate(full.comp_ids)}
    for j in range(3):                                                        # existing members expand to identical compound masses
        for i, c in enumerate(base.comp_ids):
            assert full.A[ids[c], j] == base.A[i, j]
    u = np.array([0.6, 0.2, 0.1, 0.1])
    assert abs(full.expand(u).sum() - 1.0) < 1e-12 and abs(full.synthetic_mass(u) - 0.4) < 1e-12
