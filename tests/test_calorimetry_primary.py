"""Exactness / consistency checks for data/raw/calorimetry_primary (transcribed from four primary papers)."""
import json
import math
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

D = Path(__file__).resolve().parents[1] / "data" / "raw" / "calorimetry_primary"
VERIFIED = Path(__file__).resolve().parents[1] / "data" / "raw" / "sle_constants" / "tm_dhfus_verified.csv"


@pytest.fixture(scope="module")
def pc():
    return pd.read_csv(D / "phase_change_constants.csv")


def test_inchikeys_match_smiles(pc):
    for _, r in pc.iterrows():
        assert Chem.MolToInchiKey(Chem.MolFromSmiles(r.smiles)) == r.inchikey, r["name"]


def test_finke_exactness_against_repo_values(pc):
    """cal x 4.184 reproduces the repo's existing values to the last digit (misreads would break this)."""
    v = pd.read_csv(VERIFIED).set_index("inchikey")
    expect = {"n-decane": 28.715, "n-undecane": 22.179, "n-dodecane": 36.836, "n-tridecane": 28.501, "n-tetradecane": 45.070,
              "n-pentadecane": 34.593, "n-hexadecane": 53.359, "n-octane": 20.740}
    f = pc[pc.paper == "F54"].set_index("name")
    for n, val in expect.items():
        assert f.loc[n, "dhfus_kJ_mol"] == pytest.approx(val, abs=1e-9), n
        assert v.loc[f.loc[n, "inchikey"], "dhfus_kJ_mol"] == pytest.approx(val, abs=1e-9), n
    assert f.loc["n-undecane", "trans1_dH_kJ_mol"] == pytest.approx(6.858, abs=1e-9)     # 1639 cal at 236.6 K
    assert f.loc["n-tridecane", "trans1_dH_kJ_mol"] == pytest.approx(7.661, abs=1e-9)    # 1831 cal at 255.0 K
    assert f.loc["n-pentadecane", "trans1_dH_kJ_mol"] == pytest.approx(9.167, abs=1e-9)  # 2191 cal at 270.9 K
    assert f.loc["n-nonane", "trans1_dH_kJ_mol"] == pytest.approx(6.280, abs=1e-9)       # 1501 cal at 217.2 K
    assert f.loc["n-undecane", "dhfus_total_incl_trans_kJ_mol"] == pytest.approx(29.037, abs=1e-9)
    assert f.loc["n-tridecane", "dhfus_total_incl_trans_kJ_mol"] == pytest.approx(36.162, abs=1e-9)


def test_other_papers_against_repo(pc):
    v = pd.read_csv(VERIFIED).set_index("inchikey")
    for n, dh in {"toluene": 6.636, "o-xylene": 13.598, "2,2-dimethylbutane": 0.579}.items():
        r = pc[pc["name"] == n].iloc[0]
        assert r.dhfus_kJ_mol == pytest.approx(dh, abs=0.0011), n
        assert v.loc[r.inchikey, "dhfus_kJ_mol"] == pytest.approx(r.dhfus_kJ_mol, abs=0.0011), n


def test_3_methylpentane_left_blank(pc):
    r = pc[pc["name"] == "3-methylpentane"].iloc[0]
    assert pd.isna(r.tm_K) and pd.isna(r.dhfus_kJ_mol)


def test_triple_points_ordered_and_transitions_below_melt(pc):
    for _, r in pc.dropna(subset=["tm_K"]).iterrows():
        for i in (1, 2):
            t = r[f"trans{i}_T_K"]
            if not pd.isna(t):
                assert t < r.tm_K, (r["name"], i)
    assert (pc.dropna(subset=["dhfus_kJ_mol"]).dhfus_kJ_mol > 0).all()


def test_vapor_pressure_equations_reproduce_table_vii():
    """The paper's own fits (eqs 1-5) must reproduce its calc column to ~1 %; guards the transcription of both."""
    vp = pd.read_csv(D / "vapor_pressure_pitzer_scott1943.csv")
    eq = pd.read_csv(D / "vapor_pressure_equations_pitzer_scott1943.csv")
    for _, r in vp.iterrows():
        e = eq[(eq["name"] == r["name"]) & (eq.phase == r.phase)].iloc[0]
        T = r.T_K_published_basis
        p = 10 ** (e.a / T + e.b * math.log10(T) + e.c)
        assert p == pytest.approx(r.p_calc_mmHg, rel=0.012), (r["name"], r.T_C_published)
        assert r.p_obs_mmHg == pytest.approx(r.p_calc_mmHg, rel=0.06)
        assert r.p_obs_kPa == pytest.approx(r.p_obs_mmHg * 0.133322368, abs=6e-5)


def test_toluene_dhvap_equation_matches_table():
    dv = pd.read_csv(D / "dhvap.csv")
    t = dv[(dv["name"] == "toluene") & (dv.paper == "S62") & dv.T_K.notna()]
    for _, r in t.iterrows():
        eq2 = 11637 - 4.823 * r.T_K - 1.260e-2 * r.T_K ** 2
        assert eq2 == pytest.approx(r.dHvap_cal_mol, abs=40), r.T_K  # eq (2) fits to a few cal; 298-410 K validity


def test_xylene_eutectics_schema_and_physics():
    e = pd.read_csv(D / "xylene_eutectics_pitzer_scott1943.csv")
    for _, r in e.iterrows():
        x = json.loads(r.x_liquid)
        assert sum(x) == pytest.approx(1.0, abs=2e-3)
        assert len(x) == len(json.loads(r.inchikeys)) == r.n_comp
    ex = e[e.point_kind == "eutectic_measured"]
    assert len(ex) == 3 and ex.T_K.min() > 200 and ex.T_K.max() < 225
    ternary = e[(e.n_comp == 3) & (e.point_kind == "eutectic_measured")].T_K.iloc[0]
    binary = e[(e.n_comp == 2) & (e.point_kind == "eutectic_measured")].T_K
    assert ternary < binary.min()                      # ternary eutectic is below every binary eutectic
    pure = {"m-xylene": 225.27, "p-xylene": 286.39}
    assert ex.T_K.max() < min(pure.values())            # every eutectic is below the lowest pure melting point


def test_new_compounds_in_verified_and_integrated():
    v = pd.read_csv(VERIFIED).set_index("name")
    for n in ("n-nonane", "2,3-dimethylbutane", "m-xylene", "p-xylene"):
        assert n in v.index and pd.notna(v.loc[n, "tm_K"]) and pd.notna(v.loc[n, "dhfus_kJ_mol"])
    assert v.loc["n-nonane", "svl_dhfus_kJ_mol"] == pytest.approx(21.748, abs=1e-9)
    m = Path(__file__).resolve().parents[1] / "data" / "unified" / "measurements_dedup.csv"
    if m.exists():  # generated, gitignored
        d = pd.read_csv(m, low_memory=False)
        assert d.sources.astype(str).str.contains("calorimetry_primary").sum() >= 21
        assert not d.duplicated(["property", "system_id", "x_key", "T_key", "P_key"]).any()
