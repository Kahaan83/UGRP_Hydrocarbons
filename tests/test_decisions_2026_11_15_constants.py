import io
import subprocess

import pandas as pd

from fuelblend.eval.composition_route import build_tmap_verified


def test_cis_decalin_mw_verified_formula():
    pp = pd.read_csv("data/processed/v5/pure_props_v5.csv")
    r = pp[(pp.compound == "decalin_cis") & (pp.prop == "MW")].iloc[0]
    assert r.status == "verified" and abs(r.ref_value - 138.25) < 0.01 and "C10H18" in r.ref_source


def test_pure_props_only_cis_decalin_mw_changed():
    old = pd.read_csv(io.StringIO(subprocess.check_output(
        ["git", "show", "eee2217:data/processed/v5/pure_props_v5.csv"], text=True)))
    new = pd.read_csv("data/processed/v5/pure_props_v5.csv")
    diff = old.merge(new, how="outer", indicator=True).query("_merge != 'both'")
    assert set(zip(diff.compound, diff.prop)) == {("decalin_cis", "MW")}


def test_bicyclohexyl_svl_dhfus_total():
    t = build_tmap_verified()["WVIIMZNLDWSIRH-UHFFFAOYSA-N"]
    assert t["dh_fus_J_mol"] == 13940.0 and t["tm_K"] == 276.8
