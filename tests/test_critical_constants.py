"""Unit test for T39 and T40: critical constants and derived parameters."""
import os
import pandas as pd
import numpy as np


def test_critical_constants_ranges():
    csv_path = "data/processed/critical_constants_hc.csv"
    assert os.path.exists(csv_path), f"File {csv_path} not found"

    df = pd.read_csv(csv_path)
    assert len(df) == 1038, f"Expected 1038 rows, got {len(df)}"

    required_cols = [
        'inchikey', 'smiles', 'formula', 'mw', 'tc_K', 'pc_kPa', 'vc_m3_mol',
        'omega', 'v_star_m3_mol', 'source_tc', 'source_pc', 'source_omega', 'tier', 'derived'
    ]
    for c in required_cols:
        assert c in df.columns, f"Missing required column {c}"

    # Check n-octane (SMILES: CCCCCCCC)
    octane = df[df['smiles'] == 'CCCCCCCC'].iloc[0]
    assert np.isclose(octane['tc_K'], 568.9, atol=2.0), f"Octane Tc out of range: {octane['tc_K']}"
    assert np.isclose(octane['pc_kPa'], 2490.0, atol=20.0), f"Octane Pc out of range: {octane['pc_kPa']}"
    assert np.isclose(octane['omega'], 0.398, atol=0.03), f"Octane omega out of range: {octane['omega']}"
    assert octane['tier'] == 'A', f"Octane should be tier A, got {octane['tier']}"
    assert octane['derived'] == True, "Octane derived should be True"

    # Check benzene (c1ccccc1)
    benzene = df[df['smiles'] == 'c1ccccc1'].iloc[0]
    assert np.isclose(benzene['tc_K'], 562.0, atol=3.0), f"Benzene Tc out of range: {benzene['tc_K']}"
    assert np.isclose(benzene['pc_kPa'], 4890.0, atol=50.0), f"Benzene Pc out of range: {benzene['pc_kPa']}"
    assert np.isclose(benzene['omega'], 0.211, atol=0.03), f"Benzene omega out of range: {benzene['omega']}"

    # Check toluene (Cc1ccccc1)
    toluene = df[df['smiles'] == 'Cc1ccccc1'].iloc[0]
    assert np.isclose(toluene['tc_K'], 593.0, atol=3.0), f"Toluene Tc out of range: {toluene['tc_K']}"
    assert np.isclose(toluene['pc_kPa'], 4100.0, atol=50.0), f"Toluene Pc out of range: {toluene['pc_kPa']}"
    assert np.isclose(toluene['omega'], 0.266, atol=0.03), f"Toluene omega out of range: {toluene['omega']}"

    # Check reasonable general physical ranges for populated entries
    valid_tc = df['tc_K'].dropna()
    assert (valid_tc > 150.0).all() and (valid_tc < 1300.0).all(), "Tc values outside physical bounds"

    valid_pc = df['pc_kPa'].dropna()
    assert (valid_pc > 100.0).all() and (valid_pc < 10000.0).all(), "Pc values outside physical bounds"

    valid_om = df['omega'].dropna()
    assert (valid_om > -0.2).all() and (valid_om < 2.5).all(), "Omega values outside physical bounds"

    # Check COSTALD V* formula holds: V* = Vc * (0.29056 - 0.08775 * omega)
    valid_vstar = df[df['v_star_m3_mol'].notna()]
    expected_vstar = valid_vstar['vc_m3_mol'] * (0.29056 - 0.08775 * valid_vstar['omega'])
    assert np.allclose(valid_vstar['v_star_m3_mol'], expected_vstar, rtol=1e-4), "COSTALD V* mismatch"

    # Check high coverage
    assert df['tc_K'].notna().sum() >= 930
    assert df['pc_kPa'].notna().sum() >= 930
    assert df['vc_m3_mol'].notna().sum() >= 930
    assert df['omega'].notna().sum() >= 900
    assert df['v_star_m3_mol'].notna().sum() >= 900

    print("test_critical_constants passed successfully!")


if __name__ == '__main__':
    test_critical_constants_ranges()
