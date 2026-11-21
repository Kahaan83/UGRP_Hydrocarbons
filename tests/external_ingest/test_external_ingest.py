import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

import kaust_registry as kr
from _common import STAGING_DIR, REPO, system_id

CANON = ["kaust_jameel_canonical_measurements.csv", "canonical_ai4fuel_measurements.csv",
         "canonical_leenhouts_mon_mix_measurements.csv"]


def _ik(smi):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smi))


def test_registry_inchikeys_match_rdkit_and_known_values():
    keys = kr.component_inchikeys()
    for name, c in kr.COMPONENTS.items():
        assert keys[name] == _ik(c["smiles"])
    # independent anchors (also present in data/unified/compounds.csv)
    comp = set(pd.read_csv(REPO / "data/unified/compounds.csv", usecols=["inchikey"])["inchikey"])
    for name in ["n-heptane", "i-octane", "toluene", "1,2,4-trimethylbenzene", "cyclopentane", "1-hexene"]:
        assert keys[name] in comp, name
    assert keys["1,2,4-trimethylbenzene"] == "GWHJZXXIDMPWGX-UHFFFAOYSA-N"
    assert keys["1-hexene"] == "LIKMAJRDDDTEIG-UHFFFAOYSA-N"
    assert keys["diisobutylene"] == "FXNDIJDIPNCZQJ-UHFFFAOYSA-N"


def test_table3_all_128_mapped_and_valid():
    assert sorted(kr.TABLE3_SMILES) == list(range(1, 129))
    for s in kr.TABLE3_SMILES.values():
        assert Chem.MolFromSmiles(s) is not None


def test_resolved_face_inchikeys_are_rdkit_keys():
    df = pd.read_csv(STAGING_DIR / "canonical_resolved_face_gasoline.csv")
    valid = set(kr.component_inchikeys().values())
    for s in df.loc[df.category == "surrogate_blend", "comp_inchikeys"]:
        assert set(json.loads(s)) <= valid


def test_kaust_pure_values_equal_plus0_rows():
    t4 = pd.read_csv(STAGING_DIR / "kaust_jameel_table4_blends.csv")
    pp = kr.pure_props()
    for comp, row in [("1,2,4-trimethylbenzene", "TMB+0"), ("cyclopentane", "CP+0"), ("1-hexene", "Hex+0")]:
        r = t4[t4.blend_name == row].iloc[0]
        assert pp[comp]["ron"] == r.ron and pp[comp]["mon"] == r.mon
    assert (pp["1,2,4-trimethylbenzene"]["ron"], pp["1,2,4-trimethylbenzene"]["mon"]) == (109.5, 108.0)


def test_dib_present_in_all_mc_blends():
    m = pd.read_csv(STAGING_DIR / "kaust_jameel_canonical_measurements.csv")
    dib = kr.component_inchikeys()["diisobutylene"]
    mc = m[m.fuel_id.astype(str).str.startswith("MC ") & (m.property == "ron")]
    assert len(mc) == 26
    assert all(dib in json.loads(s) for s in mc.comp_inchikeys)
    assert (mc.vol_sum_pct.sub(100).abs() < 0.2).all()


def test_kaust_kinds_and_oxygenate_tag():
    m = pd.read_csv(STAGING_DIR / "kaust_jameel_canonical_measurements.csv")
    assert (m[m.n_comp == 1].kind == "pure").all()
    assert not (m[m.fuel_id.astype(str).str.startswith("TMB+")].kind == "multi").any()
    eth = kr.component_inchikeys()["ethanol"]
    has_eth = m.comp_inchikeys.map(lambda s: eth in json.loads(s))
    assert not m[has_eth].is_hc.any()
    assert m[m.fuel_id.astype(str).str.startswith("FACE")].shape[0] == 60  # Table 5, RON+MON


@pytest.mark.parametrize("name", CANON)
def test_system_id_convention(name):
    df = pd.read_csv(STAGING_DIR / name)
    named = df[~df.system_id.astype(str).str.startswith("fuel:")]
    assert len(named) > 0
    for iks, sid in zip(named.comp_inchikeys.head(400), named.system_id.head(400)):
        want = hashlib.md5("|".join(sorted(json.loads(iks))).encode()).hexdigest()[:12]
        assert sid == want == system_id(json.loads(iks))


def test_ai4fuel_corrections():
    m = pd.read_csv(STAGING_DIR / "canonical_ai4fuel_measurements.csv")
    assert (m.provenance_tier == "compiled").all()
    assert "cetane_number_D613" not in set(m.property) and "cetane_number" in set(m.property)
    assert "melting_point_C" in set(m.property) and "freezing_point_C" not in set(m.property)
    d = m[m.property == "density_kg_m3"]
    assert d.T_uncertain.all() and (d.T_K_lo == 288.15).all() and (d.T_K_hi == 303.15).all() and d.T_K.isna().all()
    assert {"dup_in_unified", "compound_in_unified"} <= set(m.columns)
    c = pd.read_csv(STAGING_DIR / "canonical_ai4fuel_compounds.csv")
    assert int(c.in_unified_compounds.sum()) == 603
    from rdkit.Chem.Descriptors import MolWt
    r = c.iloc[0]
    assert abs(r.mw - MolWt(Chem.MolFromSmiles(r.smiles))) < 1e-6


def test_dropped_staging_and_literals_removed():
    d = REPO / "scripts/external_ingest"
    for f in ["parse_boehm2022_cold_viscosity.py", "ingest_crc_face1_diesel.py", "ingest_ternary_flash.py"]:
        assert not (d / f).exists()
    assert "FACE_BASE_PIONA" not in (d / "resolve_face_gasoline_dha.py").read_text(encoding="utf-8")
