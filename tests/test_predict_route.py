"""Task 1 (v5 api route): route context and per-property component mapping."""
import os

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import inchi

from fuelblend.api import predict_route as pr
from fuelblend.baselines.flash_unifac import _lookup
from fuelblend.eval.composition_route import kind_bins, map_bins

DECANE = "CCCCCCCCCC"
DODECANE = "CCCCCCCCCCCC"
TOLUENE = "Cc1ccccc1"


def _ik(smi):
    return inchi.MolToInchiKey(Chem.MolFromSmiles(smi))


@pytest.fixture(scope="module")
def ctx():
    return pr.get_ctx()


def test_resolve_smiles_and_inchikey_agree(ctx):
    a = pr.resolve({DECANE: 1.0})
    b = pr.resolve({_ik(DECANE): 1.0})
    assert len(a) == len(b) == 1 and a[0] == b[0]
    assert a[0].w == pytest.approx(1.0) and a[0].is_hc and a[0].n_carbon == 10
    assert a[0].subkind == "n_paraffin" and a[0].reg_idx is not None


def test_resolve_merges_duplicates(ctx):
    c = pr.resolve([(DECANE, 1.0), (_ik(DECANE), 1.0), ("C(CCCCCCCC)C", 2.0), (TOLUENE, 0.0), (DODECANE, 4.0)])
    assert len(c) == 2
    d = {x.n_carbon: x.w for x in c}
    assert d[10] == pytest.approx(0.5) and d[12] == pytest.approx(0.5)


def test_resolve_bad_smiles_raises():
    with pytest.raises(ValueError):
        pr.resolve({"not_a_smiles(((": 1.0})
    with pytest.raises(ValueError):
        pr.resolve({"AAAAAAAAAAAAAA-BBBBBBBBBB-N": 1.0})


def test_map_for_own_data_no_substitution(ctx):
    comps = pr.resolve({DECANE: 0.4, DODECANE: 0.4, TOLUENE: 0.2})
    p = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert p.substituted_mass == 0 and p.unmapped_mass == 0 and p.own_mass == pytest.approx(1.0)
    assert p.w.sum() == pytest.approx(1.0)


def test_map_for_substitutes_missing(ctx):
    odd = "CC(C)CCCC(C)c1ccccc1"          # hydrocarbon absent from the registry, so no grid
    comps = pr.resolve({DECANE: 0.5, odd: 0.5})
    odd_c = [c for c in comps if c.ik == _ik(odd)][0]
    assert odd_c.reg_idx is None and odd_c.is_hc
    p = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert len(p.substitutions) == 1
    a, b, gap = p.substitutions[0]
    assert a == odd_c.ik and b != a and gap <= 3
    assert p.substituted_mass == pytest.approx(0.5) and p.own_mass == pytest.approx(0.5)


def test_non_hc_unmapped(ctx):
    comps = pr.resolve({DECANE: 0.8, "CCO": 0.2})
    p = pr.map_for(comps, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert p.unmapped_mass == pytest.approx(0.2) and p.w.sum() == pytest.approx(1.0)


def test_map_for_route_equivalence(ctx):
    root = pr.ROOT
    fu = pd.read_csv(root / "results/v5/BASE/t2_fuels.csv")
    fu = fu[~fu.low_coverage]
    V5 = root / "data/real_fuels/v5"
    srcs = {"av3322": "crc_av3322_composition.csv", "fuellib": "fuellib_gcxgc_compositions.csv", "ajftd": "ajftd_gcxgc_composition.csv"}
    comps = {}
    n = 0
    for r in fu.itertuples():
        if r.source not in srcs:
            continue
        if r.source not in comps:
            comps[r.source] = pd.read_csv(V5 / srcs[r.source])
        fid = int(r.fuel_id) if r.source == "ajftd" else r.fuel_id
        c = comps[r.source]
        mf = map_bins(kind_bins(r.source, c[c.fuel_id == fid]), ctx.reps)
        cs = pr.resolve({ctx.reg.ik[int(i)]: float(w) for i, w in zip(mf.idx, mf.w)})
        p = pr.map_for(cs, pr.predicates(ctx)["HAS_FLASH"], ctx.reps)
        assert np.array_equal(np.sort(p.idx), np.sort(mf.idx))
        order = np.argsort(p.idx)
        assert np.allclose(p.w[order], mf.w[np.argsort(mf.idx)])
        n += 1
        if n == 5:
            break
    assert n == 5


def test_provider_never_writes():
    path = pr.GRID_NPZ
    m0 = os.stat(path).st_mtime_ns
    pr.get_ctx.cache_clear()
    c = pr.get_ctx()
    c.prov.lookup(np.array([[0]]), 300.0)
    c.prov.ensure([0, 1, 2])
    assert os.stat(path).st_mtime_ns == m0


def test_provider_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        pr.ReadOnlyProvider(pr.get_ctx().reg, tmp_path / "nope.npz")


def test_non_hc_with_grid_is_unmapped(ctx):
    etbe = "CCOC(C)(C)C"
    c = pr.resolve({DECANE: 0.8, etbe: 0.2})
    e = [x for x in c if not x.is_hc][0]
    if e.reg_idx is None or not ctx.prov.valid[e.reg_idx]:
        pytest.skip("ETBE has no grid in this registry")
    p = pr.map_for(c, pr.predicates(ctx)["HAS_GRID"], ctx.reps)
    assert p.unmapped_mass == pytest.approx(0.2) and p.own_mass == pytest.approx(0.8)


def test_rep_without_own_data_is_unmapped(ctx):
    P = pr.predicates(ctx)
    cand = [i for i in range(len(ctx.reg.ik)) if P["HAS_GRID"](i) and P["HAS_FLASH"](i) and not P["HAS_SLE"](i)
            and pr.subkind(ctx.reg.smiles[i]) is not None]
    if not cand:
        pytest.skip("no grid+flash compound lacking SLE")
    c = pr.resolve({ctx.reg.ik[cand[0]]: 1.0})
    p = pr.map_for(c, P["HAS_SLE"], ctx.reps)
    assert p.substituted_mass == 0 and p.unmapped_mass == pytest.approx(1.0) and p.idx.size == 0


def test_resolve_mole_fractions(ctx):
    c = pr.resolve({DECANE: 0.5, DODECANE: 0.5}, mole_fractions=True)
    d = {x.n_carbon: x.w for x in c}
    m10, m12 = ctx.mw(c[0].smiles if c[0].n_carbon == 10 else c[1].smiles), ctx.mw(c[0].smiles if c[0].n_carbon == 12 else c[1].smiles)
    assert d[10] == pytest.approx(m10 / (m10 + m12))
