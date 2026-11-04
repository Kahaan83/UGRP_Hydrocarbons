"""
Extract critical constants for all compounds in data/unified/compounds.csv (and registry)
using `chemicals` library with priority:
NIST/WebBook (HEOS, WEBBOOK, IUPAC, CRC) > Poling/Perry (PSRK, PD, PINAMARTINES) > Yaws (YAWS) > Joback/Wilson.
"""

import math
import os
import pandas as pd
import chemicals

COMPOUNDS_CSV = 'data/unified/compounds.csv'
REGISTRY_CSV = 'compounds.csv'
OUTPUT_RAW = 'data/raw/chemicals_lib/critical.csv'
OUTPUT_PROCESSED = 'data/processed/critical_constants.csv'

os.makedirs('data/raw/chemicals_lib', exist_ok=True)
os.makedirs('data/processed', exist_ok=True)

df_compounds = pd.read_csv(COMPOUNDS_CSV)
if os.path.exists(REGISTRY_CSV):
    df_reg = pd.read_csv(REGISTRY_CSV)
else:
    df_reg = pd.DataFrame()

# Methods priority lists per property
tc_pref = ['HEOS', 'WEBBOOK', 'IUPAC', 'CRC', 'PSRK', 'PD', 'PINAMARTINES', 'YAWS', 'WILSON_JASPERSON', 'JOBACK']
pc_pref = ['HEOS', 'WEBBOOK', 'IUPAC', 'CRC', 'PSRK', 'PD', 'PINAMARTINES', 'YAWS', 'WILSON_JASPERSON', 'JOBACK']
vc_pref = ['HEOS', 'WEBBOOK', 'IUPAC', 'CRC', 'PSRK', 'PINAMARTINES', 'YAWS', 'JOBACK', 'FEDORS']
omega_pref = ['HEOS', 'PSRK', 'PD', 'YAWS', 'ACENTRIC_DEFINITION']

def resolve_cas(row):
    for id_val in ['InChIKey=' + str(row.get('inchikey', '')), str(row.get('inchi', '')), str(row.get('smiles', ''))]:
        if not id_val or id_val == 'nan' or id_val == 'InChIKey=nan':
            continue
        try:
            cas = chemicals.CAS_from_any(id_val)
            if cas:
                return cas
        except Exception:
            pass
    return None

records = []
seen_inchikeys = set()

for idx, r in df_compounds.iterrows():
    inchikey = r.get('inchikey')
    seen_inchikeys.add(inchikey)
    cas = resolve_cas(r)
    
    rec = {
        'inchikey': inchikey,
        'smiles': r.get('smiles'),
        'formula': r.get('formula'),
        'mw': r.get('mw'),
        'cas': cas,
        'Tc_K': None,
        'Tc_source': None,
        'Pc_kPa': None,
        'Pc_source': None,
        'Vc_m3_mol': None,
        'Vc_source': None,
        'omega': None,
        'omega_source': None,
    }
    
    if cas:
        tc_avail = chemicals.critical.Tc_methods(cas)
        tc_m = next((m for m in tc_pref if m in tc_avail), None)
        if tc_m:
            rec['Tc_K'] = chemicals.critical.Tc(cas, method=tc_m)
            rec['Tc_source'] = tc_m
            
        pc_avail = chemicals.critical.Pc_methods(cas)
        pc_m = next((m for m in pc_pref if m in pc_avail), None)
        if pc_m:
            p = chemicals.critical.Pc(cas, method=pc_m)
            rec['Pc_kPa'] = p / 1000.0 if p is not None else None
            rec['Pc_source'] = pc_m
            
        vc_avail = chemicals.critical.Vc_methods(cas)
        vc_m = next((m for m in vc_pref if m in vc_avail), None)
        if vc_m:
            rec['Vc_m3_mol'] = chemicals.critical.Vc(cas, method=vc_m)
            rec['Vc_source'] = vc_m
            
        om_avail = chemicals.acentric.omega_methods(cas)
        om_m = next((m for m in omega_pref if m in om_avail), None)
        if om_m:
            rec['omega'] = chemicals.acentric.omega(cas, method=om_m)
            rec['omega_source'] = om_m
            
    records.append(rec)

res_df = pd.DataFrame(records)
res_df.to_csv(OUTPUT_RAW, index=False)
res_df.to_csv(OUTPUT_PROCESSED, index=False)

print(f"T39 Exported {len(res_df)} critical property rows to {OUTPUT_PROCESSED}")
print("Coverage counts:")
print(res_df[['cas', 'Tc_K', 'Pc_kPa', 'Vc_m3_mol', 'omega']].notna().sum())
