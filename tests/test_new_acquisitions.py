"""Unit tests for newly generated datasets:
- data/processed/octane_ron_mon_hc.csv (T16)
- data/processed/freezing_point_pure_hc.csv (T17)
- data/processed/surface_tension_hc.csv (T20)
- data/processed/sooting_real_fuels.csv (T23)
"""
import os, re
import pandas as pd

def is_hc(formula):
    return bool(isinstance(formula, str) and re.fullmatch(r'C(\d*)H(\d*)', formula))

def test_octane_dataset():
    p = 'data/processed/octane_ron_mon_hc.csv'
    assert os.path.exists(p), f"Missing {p}"
    df = pd.read_csv(p)
    assert len(df) >= 400, f"Expected >= 400 rows, got {len(df)}"
    assert all(df['formula'].apply(is_hc)), "All compounds must be pure hydrocarbons"
    assert set(df['tier'].unique()).issubset({'A', 'B', 'C'}), "Tiers must be A, B, or C"
    # Isooctane check: CAS 540-84-1 RON=100, MON=100
    iso = df[df['cas'] == '540-84-1']
    assert len(iso) == 1, "Isooctane not found"
    assert float(iso['ron'].iloc[0]) == 100.0, f"Isooctane RON should be 100, got {iso['ron'].iloc[0]}"
    assert float(iso['mon'].iloc[0]) == 100.0, f"Isooctane MON should be 100, got {iso['mon'].iloc[0]}"
    # n-heptane check: CAS 142-82-5 RON=0, MON=0
    hep = df[df['cas'] == '142-82-5']
    assert len(hep) == 1, "n-Heptane not found"
    assert float(hep['ron'].iloc[0]) == 0.0, f"n-Heptane RON should be 0, got {hep['ron'].iloc[0]}"
    assert float(hep['mon'].iloc[0]) == 0.0, f"n-Heptane MON should be 0, got {hep['mon'].iloc[0]}"
    print("test_octane_dataset PASSED")

def test_freezing_point_dataset():
    p = 'data/processed/freezing_point_pure_hc.csv'
    assert os.path.exists(p), f"Missing {p}"
    df = pd.read_csv(p)
    assert len(df) >= 600, f"Expected >= 600 rows, got {len(df)}"
    assert all(df['formula'].apply(is_hc)), "All compounds must be pure hydrocarbons"
    assert df['unit'].iloc[0] == 'degC', "Canonical unit must be degC"
    # Physical range check
    assert df['freezing_point_C'].min() >= -200.0, f"Freezing point minimum too low: {df['freezing_point_C'].min()}"
    assert df['freezing_point_C'].max() <= 500.0, f"Freezing point maximum too high: {df['freezing_point_C'].max()}"
    # Benzene check: CAS 71-43-2 Tm ~ 5.5 C (278.65 K)
    bz = df[df['cas'] == '71-43-2']
    assert len(bz) == 1, "Benzene not found"
    assert abs(float(bz['freezing_point_C'].iloc[0]) - 5.5) < 1.0, f"Benzene Tm should be ~5.5 C, got {bz['freezing_point_C'].iloc[0]}"
    print("test_freezing_point_dataset PASSED")

def test_surface_tension_dataset():
    p = 'data/processed/surface_tension_hc.csv'
    assert os.path.exists(p), f"Missing {p}"
    df = pd.read_csv(p)
    assert len(df) >= 100, f"Expected >= 100 rows, got {len(df)}"
    assert all(df['formula'].apply(is_hc)), "All compounds must be pure hydrocarbons"
    assert df['unit'].iloc[0] == 'mN/m', "Canonical unit must be mN/m"
    # Range check on liquid phase hydrocarbons at 20 C
    liq = df[df['is_liquid_phase_at_20C'] & (df['eval_method'] != 'supercritical_T_ge_Tc')]
    assert len(liq) >= 100, "Expected >= 100 liquid-phase hydrocarbons at 20 C"
    # Butane and heavier liquids have sigma between ~10 and ~45 mN/m at 20 C
    assert liq['surface_tension_20C_mN_m'].min() >= 0.0, "Surface tension must be >= 0"
    assert liq['surface_tension_20C_mN_m'].max() <= 100.0, "Surface tension too high"
    # n-Octane check: CAS 111-65-9 sigma ~ 21.6 mN/m
    octane = df[df['cas'] == '111-65-9']
    assert len(octane) == 1, "Octane not found"
    assert abs(float(octane['surface_tension_20C_mN_m'].iloc[0]) - 21.6) < 1.5, f"Octane sigma should be ~21.6 mN/m, got {octane['surface_tension_20C_mN_m'].iloc[0]}"
    print("test_surface_tension_dataset PASSED")

def test_sooting_dataset():
    p = 'data/processed/sooting_real_fuels.csv'
    assert os.path.exists(p), f"Missing {p}"
    df = pd.read_csv(p)
    assert len(df) == 100, f"Expected 100 fuels, got {len(df)}"
    assert df['unit'].iloc[0] == 'mm', "Canonical unit must be mm"
    assert df['method'].iloc[0] in ['ASTM D1322', 'D1322'], "Method should be ASTM D1322"
    # Jet fuel smoke point range ~ 15-35 mm
    assert df['smoke_point_mm'].min() >= 15.0, "Smoke point too low"
    assert df['smoke_point_mm'].max() <= 35.0, "Smoke point too high"
    print("test_sooting_dataset PASSED")

if __name__ == '__main__':
    test_octane_dataset()
    test_freezing_point_dataset()
    test_surface_tension_dataset()
    test_sooting_dataset()
    print("ALL TESTS PASSED SUCCESSFULLY!")
