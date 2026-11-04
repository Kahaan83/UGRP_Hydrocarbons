import os
import pandas as pd
import pytest

DENSITY_CSV = "data/processed/nist_density_compressed_liquid.csv"

def test_nist_density_exists():
    assert os.path.exists(DENSITY_CSV)

def test_nist_density_schema():
    df = pd.read_csv(DENSITY_CSV)
    required = ["fuel_id", "fuel_name", "category", "T_K", "P_MPa", "P_kPa", "density_kg_m3", "uncertainty_kg_m3", "source_table", "reference"]
    for col in required:
        assert col in df.columns

def test_nist_density_fuels():
    df = pd.read_csv(DENSITY_CSV)
    fuels = set(df["fuel_id"].unique())
    expected = {"Jet-A-3602", "Jet-A-3638", "Jet-A-4658", "JP-8-3773", "S-8"}
    assert fuels == expected
    assert len(df) == 634

def test_nist_density_physical_ranges():
    df = pd.read_csv(DENSITY_CSV)
    assert (df["T_K"] >= 270.0).all() and (df["T_K"] <= 470.0).all()
    assert (df["P_MPa"] >= 0.4).all() and (df["P_MPa"] <= 41.0).all()
    assert (df["density_kg_m3"] >= 600.0).all() and (df["density_kg_m3"] <= 860.0).all()
