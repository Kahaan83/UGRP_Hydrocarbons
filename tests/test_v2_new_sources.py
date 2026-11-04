from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CETANE_CSV = ROOT / "data" / "processed" / "v2" / "cetane_nrel_compendium.csv"
OCTANE_CSV = ROOT / "data" / "processed" / "v2" / "octane_co_optima.csv"
WAX_CSV = ROOT / "data" / "processed" / "v2" / "wax_dauphin_coutinho.csv"
SOURCES_CSV = ROOT / "docs" / "data_sources.csv"


def test_new_sources_files_exist():
    assert CETANE_CSV.exists(), f"Missing {CETANE_CSV}"
    assert OCTANE_CSV.exists(), f"Missing {OCTANE_CSV}"
    assert WAX_CSV.exists(), f"Missing {WAX_CSV}"
    assert SOURCES_CSV.exists(), f"Missing {SOURCES_CSV}"


def test_nrel_cetane_compendium():
    df = pd.read_csv(CETANE_CSV)
    assert len(df) >= 200
    for col in ["compound_name", "smiles", "inchikey", "mw", "n_carbon", "dcn"]:
        assert col in df.columns
    assert df["dcn"].notna().all()
    assert (df["dcn"] > -10.0).all() and (df["dcn"] < 200.0).all()


def test_co_optima_octane():
    df = pd.read_csv(OCTANE_CSV)
    assert len(df) >= 100
    for col in ["compound_name", "ron", "mon", "sensitivity"]:
        assert col in df.columns
    # Check valid numeric values
    has_ron = df["ron"].notna()
    assert has_ron.sum() > 50


def test_coutinho_wax_sle():
    df = pd.read_csv(WAX_CSV)
    assert len(df) >= 10
    for col in ["system_id", "system_name", "n_comp", "wat_measured_K", "wat_measured_C"]:
        assert col in df.columns
    assert (df["wat_measured_K"] > 200.0).all()


def test_data_sources_registry():
    df = pd.read_csv(SOURCES_CSV)
    sids = df["source_id"].tolist()
    assert "nrel_compendium_2017" in sids
    assert "doe_co_optima" in sids
    assert "coutinho_dauphin_wax" in sids
