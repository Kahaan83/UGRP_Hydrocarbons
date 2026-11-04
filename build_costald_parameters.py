"""
T40: Compute and verify COSTALD characteristic volume V* and SRK-tuned omega.
Derived=True, tier=C.
Saves to data/processed/costald_parameters.csv.
"""

import json
import math
import os
import pandas as pd
import chemicals
import chemicals.volume as cv

COMPOUNDS_CSV = 'data/unified/compounds.csv'
CRITICAL_CSV = 'data/processed/critical_constants.csv'
MEASUREMENTS_CSV = 'data/unified/measurements_dedup.csv'
OUTPUT_CSV = 'data/processed/costald_parameters.csv'

df_compounds = pd.read_csv(COMPOUNDS_CSV)
df_crit = pd.read_csv(CRITICAL_CSV).set_index('inchikey')

# Load density data to compute characteristic volume V*
df_meas = pd.read_csv(MEASUREMENTS_CSV, low_memory=False)
pure_density = df_meas[(df_meas['property'] == 'density_kg_m3') & (df_meas['kind'] == 'pure')]

# Map inchikey to list of (T_K, density_kg_m3)
density_by_inchikey = {}
for idx, r in pure_density.iterrows():
    s = r.get('comp_inchikeys')
    if pd.isna(s): continue
    try:
        keys = json.loads(s)
        if keys:
            ik = keys[0]
            val = r.get('value')
            T = r.get('T_K')
            if pd.notna(val) and pd.notna(T) and val > 0 and T > 0:
                density_by_inchikey.setdefault(ik, []).append((float(T), float(val)))
    except Exception:
        pass

rows = []

for idx, r in df_compounds.iterrows():
    ik = r['inchikey']
    smiles = r['smiles']
    formula = r.get('formula')
    mw = r.get('mw')
    
    crit_row = df_crit.loc[ik] if ik in df_crit.index else None
    Tc = crit_row['Tc_K'] if crit_row is not None and pd.notna(crit_row['Tc_K']) else None
    Pc_kPa = crit_row['Pc_kPa'] if crit_row is not None and pd.notna(crit_row['Pc_kPa']) else None
    Vc_standard = crit_row['Vc_m3_mol'] if crit_row is not None and pd.notna(crit_row['Vc_m3_mol']) else None
    omega_standard = crit_row['omega'] if crit_row is not None and pd.notna(crit_row['omega']) else None
    
    # 1. Compute SRK-tuned omega from Antoine / Psat at Tr = 0.7 if available
    omega_srk = None
    omega_srk_method = None
    
    A, B, C = r.get('antoine_A'), r.get('antoine_B'), r.get('antoine_C')
    if Tc is not None and Pc_kPa is not None and pd.notna(A) and pd.notna(B) and pd.notna(C):
        try:
            T_07 = 0.7 * float(Tc)
            # log10(P/kPa) = A - B / (T + C)
            log10_Psat = float(A) - float(B) / (T_07 + float(C))
            Psat_kPa = 10.0 ** log10_Psat
            if Psat_kPa > 0 and float(Pc_kPa) > 0:
                calc_om = -math.log10(Psat_kPa / float(Pc_kPa)) - 1.0
                if -0.5 < calc_om < 2.0:
                    omega_srk = round(calc_om, 5)
                    omega_srk_method = 'antoine_at_0.7Tc'
        except Exception:
            pass
            
    if omega_srk is None and omega_standard is not None:
        omega_srk = float(omega_standard)
        omega_srk_method = 'critical_omega_fallback'
        
    # 2. Compute Characteristic Volume V* from density measurements
    V_star = None
    V_star_method = None
    n_density_pts = 0
    
    effective_omega = omega_srk if omega_srk is not None else (omega_standard if omega_standard is not None else 0.2)
    
    if ik in density_by_inchikey and Tc is not None and pd.notna(mw) and mw > 0:
        dens_pts = density_by_inchikey[ik]
        v_star_estimates = []
        for T_meas, rho_meas in dens_pts:
            Tr = T_meas / float(Tc)
            if 0.25 <= Tr <= 0.98:
                try:
                    # COSTALD with characteristic volume = 1.0 gives the dimensionless expansion factor V_factor
                    v_factor = cv.COSTALD(T_meas, float(Tc), 1.0, float(effective_omega))
                    if v_factor > 0:
                        v_mol_meas = (float(mw) / (1000.0 * rho_meas))  # m3/mol
                        v_star_est = v_mol_meas / v_factor
                        if 1e-5 < v_star_est < 0.01:
                            v_star_estimates.append(v_star_est)
                except Exception:
                    pass
        if v_star_estimates:
            V_star = float(pd.Series(v_star_estimates).median())
            V_star_method = 'fit_to_experimental_density'
            n_density_pts = len(v_star_estimates)
            
    if V_star is None and Vc_standard is not None:
        V_star = float(Vc_standard)
        V_star_method = 'critical_Vc_fallback'
        
    rows.append({
        'inchikey': ik,
        'smiles': smiles,
        'formula': formula,
        'mw': mw,
        'Tc_K': Tc,
        'Pc_kPa': Pc_kPa,
        'omega_srk': omega_srk,
        'omega_srk_method': omega_srk_method,
        'V_star_m3_mol': V_star,
        'V_star_method': V_star_method,
        'n_density_pts': n_density_pts,
        'derived': True,
        'tier': 'C'
    })

res_df = pd.DataFrame(rows)
res_df.to_csv(OUTPUT_CSV, index=False)
print(f"T40 Exported {len(res_df)} COSTALD parameter rows to {OUTPUT_CSV}")
print("Coverage summary:")
print(res_df[['omega_srk', 'V_star_m3_mol']].notna().sum())
print("V* methods breakdown:")
print(res_df['V_star_method'].value_counts())
print("omega_srk methods breakdown:")
print(res_df['omega_srk_method'].value_counts())
