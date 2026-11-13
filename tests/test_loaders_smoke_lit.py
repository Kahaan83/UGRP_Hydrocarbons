"""Tests for the smoke point / TSI / YSI literature loader (data/raw/smoke_lit).  Page numbers are PDF page indices."""
import json
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fuelblend.data.loaders_smoke_lit import FILES, RAW, load_smoke_lit  # noqa: E402


def _one(d, sample, prop, **kw):
    s = d[(d.sample_id == sample) & (d.property == prop)]
    for k, v in kw.items():
        s = s[s[k] == v]
    assert len(s) == 1, (sample, prop, len(s))
    return s.iloc[0]


def test_row_counts():
    d = load_smoke_lit()
    assert d.groupby("property").size().to_dict() == {"smoke_point_mm": 166, "tsi": 128, "ysi": 229}
    assert d.groupby("source_id").size().to_dict() == {"S3_mensch2009": 228, "S4_rainville2025": 37,
                                                      "S5_botero_preprint153": 13, "S6_boehm2022": 16,
                                                      "S7_ysi_unified": 229}


def test_spot_values_against_pdf():
    d = load_smoke_lit()
    # S3 Table 4 (PDF p45): 1-methylnaphthalene SP 5.5 mm, TSI 100; toluene TSI 40
    assert _one(d, "pure_1-methylnaphthalene", "smoke_point_mm").value == 5.5
    assert _one(d, "pure_toluene", "tsi").value == 40
    # S3 Appendix B (PDF p72): iso-octane/toluene row 4 (0.49/0.51), SP 14.5, TSI 24
    r = _one(d, "Mensch_bin_iso-octane__toluene_4", "smoke_point_mm")
    assert r.value == 14.5 and r.fraction_basis == "mole"
    assert _one(d, "Mensch_bin_iso-octane__toluene_4", "tsi").value == 24
    # S3 Table 7 (PDF p62): 3-component A SP 22.5, TSI 29.0
    assert _one(d, "Mensch_MURI_3comp_A", "smoke_point_mm").value == 22.5
    assert _one(d, "Mensch_MURI_3comp_A", "tsi").value == 29.0
    # S3 Appendix A (PDF p70): Standard 1 first replicate 14.3 mm; 19 + 4 replicate dates for Std1 / Std6 (x2 properties)
    r = d[(d.sample_id == "Mensch_Std1_iC8_tol_60_40v") & (d.property == "smoke_point_mm")].iloc[0]
    assert r.value == 14.3 and r.qualifier == "replicate" and r.fraction_basis == "volume"
    assert (d.sample_id == "Mensch_Std1_iC8_tol_60_40v").sum() == 38
    assert (d.sample_id == "Mensch_Std6_iC8_tol_95_5v").sum() == 8
    # S3 Table 8 (PDF p63): JP-8 drum reported mean 23.0 mm
    assert _one(d, "Mensch_JP8_drum", "smoke_point_mm", qualifier="reported_mean").value == 23.0
    # S4 Table 5 (PDF p7): fuel 3 (80/20) VBA 22.37, FURTI 22.04; Table 6 (PDF p9): S8 24.1, Jet-A 22.8
    assert _one(d, "Rainville_iC8_tol_80_20v", "smoke_point_mm", method="VBA").value == 22.37
    assert _one(d, "Rainville_iC8_tol_80_20v", "smoke_point_mm", method="FURTI").value == 22.04
    assert _one(d, "Rainville2022_S8", "smoke_point_mm").value == 24.1
    assert _one(d, "JetA_Rainville2025", "smoke_point_mm").value == 22.8
    # S5 Table 1 (PDF p8): H50T50 12.0 mm, TSI 19.9; gasoline 11.8
    assert _one(d, "Botero_H50T50", "smoke_point_mm").value == 12.0 and _one(d, "Botero_H50T50", "tsi").value == 19.9
    assert _one(d, "Botero_gasoline", "smoke_point_mm").value == 11.8
    # S6 Table 5 (PDF p10): n-butylcyclohexane 45.1; surrogate 1 TSI 18.5; Table 2 (PDF p5) fuel F 14.88
    assert _one(d, "Boehm_n-butylcyclohexane", "smoke_point_mm").value == 45.1
    assert _one(d, "Boehm_surrogate1", "tsi").value == 18.5
    assert _one(d, "Boehm_T2_toluene+iso-octane_40_60v", "smoke_point_mm").value == 14.88
    # S7 YSI Volume 2 (sheet Database): benzene unified YSI 100.3
    y = d[(d.property == "ysi") & (d.components == '["benzene"]')]
    assert len(y) == 1 and abs(y.value.iloc[0] - 100.3) < 0.01


def test_smiles_parse_and_fractions():
    raw = pd.concat([pd.read_csv(RAW / f, keep_default_na=False) for f in FILES])
    for s in raw.smiles.unique():
        for smi in json.loads(s):
            assert Chem.MolFromSmiles(smi) is not None
    for r in raw.itertuples():
        fr = json.loads(r.fractions)
        if r.fraction_basis == "none":
            assert fr == [] and r.fraction_unit == ""
            continue
        assert len(fr) == len(json.loads(r.smiles))
        tot = sum(fr) / (100.0 if r.fraction_unit == "percent" else 1.0)
        # printed values are rounded; the largest gap is 0.001 (Boehm surrogate 1, MURI 4-component)
        assert abs(tot - 1) <= 1.5e-3, (r.sample_id, tot)
    d = load_smoke_lit()
    for r in d.itertuples():
        if r.fraction_basis != "none":
            assert abs(sum(json.loads(r.fractions)) - 1) < 1e-5
        assert json.loads(r.comp_inchikeys) == sorted(json.loads(r.comp_inchikeys))


def test_ids_and_flags():
    d = load_smoke_lit()
    assert d[d.kind == "real_fuel"].system_id.str.startswith("fuel:").all()
    assert d[d.kind != "real_fuel"].system_id.str.len().eq(12).all()
    assert d[d.kind != "real_fuel"].is_hc.all()
    # the same binary from two sources shares a system_id; a different pair does not
    a = d[d.sample_id == "Mensch_bin_iso-octane__toluene_1"].system_id.iloc[0]
    b = d[d.sample_id == "Botero_H50T50"].system_id.iloc[0]
    c = d[d.sample_id == "Rainville_iC8_tol_60_40v"].system_id.iloc[0]
    assert a == c and a != b
    # recipes shared with data/raw/surrogates/fuel_surrogates.csv resolve to the same set of molecules
    sur = pd.read_csv(Path(RAW).parents[0] / "surrogates" / "fuel_surrogates.csv")
    import hashlib
    from rdkit.Chem import inchi
    for code in ("Dooley2012_2ndGen", "Honnet2009_Aachen", "KimVioli2014_UM2"):
        iks = sorted(inchi.MolToInchiKey(Chem.MolFromSmiles(s)) for s in sur[sur.fuel_code == code].component_smiles)
        sid = hashlib.md5("|".join(iks).encode()).hexdigest()[:12]
        assert d[d.sample_id == code].system_id.iloc[0] == sid, code
