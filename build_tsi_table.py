"""
build_tsi_table.py - Assemble canonical TSI and smoke point reference table (T19, T20, T42).
Combines:
  T19: Boehm et al. 2022 (Energy & Fuels 36:1916) conversion constants and experimental anchors
  T20: Mensch et al. 2010 recommended surrogate TSIs and ASTM D1322 smoke points
  T42: baselines_smoke.py inputs with YSI-to-TSI mapping, Tiering (A/B/C), and censoring (>50 mm).
"""

import json
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

sys.stdout.reconfigure(encoding='utf-8')

# Output paths
OUT_CSV = 'data/processed/tsi_table.csv'
OUT_JSON = 'data/processed/tsi_constants.json'
COMPOUNDS_CSV = 'data/unified/compounds.csv'
YSI_CSV = 'data/raw/sooting/ysi_compounds.csv'

# Conversion constants from Boehm et al. 2022 (Energy & Fuels 36, 1916; Table S1 & Eq 2):
# TSI = a + b * (MW / SP)
# For Boehm 2022 apparatus:
BOEHM_A = -3.40
BOEHM_B = 3.727  # mm * mol / g

# For Mensch et al. 2010 apparatus:
MENSCH_A = -4.75
MENSCH_B = 4.051  # mm * mol / g

# Anchors defined on the TSI scale:
# Methylcyclohexane: TSI = 5.0
# 1-Methylnaphthalene: TSI = 100.0

# Tier A: Primary experimental measurements from Boehm 2022, Mensch 2010, and ASTM D1322 / Li & Sunderland
TIER_A_MEASURED = [
    # (compound_name, smiles, tsi, smoke_point_mm, source, censored, note)
    ("1-methylnaphthalene", "Cc1cccc2ccccc12", 100.0, 5.13, "boehm2022", False, "TSI anchor (100.0)"),
    ("methylcyclohexane", "CC1CCCCC1", 5.0, 43.56, "boehm2022", False, "TSI anchor (5.0)"),
    ("toluene", "Cc1ccccc1", 39.0, 8.40, "mensch2010", False, "ASTM D1322 referee / surrogate"),
    ("benzene", "c1ccccc1", 28.0, 9.10, "mensch2010", False, "Primary experimental smoke point"),
    ("m-xylene", "Cc1cccc(C)c1", 46.8, 7.88, "boehm2022", False, "Table S2 experimental"),
    ("o-xylene", "Cc1ccccc1C", 46.1, 8.00, "boehm2022", False, "Table S2 experimental"),
    ("1,3,5-trimethylbenzene", "Cc1cc(C)cc(C)c1", 58.9, 7.19, "boehm2022", False, "Boehm referee fuel (58.9)"),
    ("1,2,4-trimethylbenzene", "Cc1ccc(C)c(C)c1", 61.0, 7.20, "mensch2010", False, "Mensch recommended surrogate"),
    ("n-propylbenzene", "CCCc1ccccc1", 52.0, 8.50, "mensch2010", False, "Mensch recommended surrogate"),
    ("butylcyclohexane", "CCCCC1CCCCC1", 8.2, 45.14, "boehm2022", False, "Validation fuel Table 5"),
    ("iso-octane", "CC(C)CC(C)(C)C", 7.7, 40.00, "mensch2010", False, "Mensch surrogate / ASTM D1322 referee"),
    ("cyclohexylbenzene", "c1ccc(C2CCCCC2)cc1", 72.0, 8.90, "mensch2010", False, "Mensch recommended surrogate"),
    ("decalin", "C1CCC2CCCCC2C1", 22.0, 22.70, "mensch2010", False, "Mensch recommended surrogate"),
    ("tetralin", "c1ccc2c(c1)CCCC2", 70.0, 7.10, "mensch2010", False, "Mensch recommended surrogate"),
    ("1,2-dihydronaphthalene", "C1CC2=CC=CC=C2C=C1", 88.0, 5.80, "mensch2010", False, "Mensch recommended surrogate"),
    ("n-decane", "CCCCCCCCCC", 5.4, 50.00, "mensch2010", True, "Candle lamp limit (>50 mm)"),
    ("n-dodecane", "CCCCCCCCCCCC", 7.0, 50.00, "mensch2010", True, "Candle lamp limit (>50 mm)"),
    ("n-hexane", "CCCCCC", 2.5, 50.00, "li_sunderland2012", True, "Candle lamp limit (>50 mm)"),
    ("n-heptane", "CCCCCCC", 2.8, 50.00, "li_sunderland2012", True, "Candle lamp limit (>50 mm)"),
    ("cyclohexane", "C1CCCCC1", 4.5, 50.00, "li_sunderland2012", True, "Candle lamp limit (>50 mm)"),
    ("ethylbenzene", "CCc1ccccc1", 63.3, 5.92, "li_sunderland2012", False, "Experimental smoke point"),
    ("styrene", "C=Cc1ccccc1", 70.3, 5.27, "li_sunderland2012", False, "Experimental smoke point"),
    ("naphthalene", "c1ccc2ccccc2c1", 133.0, 3.49, "li_sunderland2012", False, "Experimental smoke point"),
    ("2-methylnaphthalene", "Cc1ccc2ccccc2c1", 117.9, 4.37, "li_sunderland2012", False, "Experimental smoke point"),
    ("cumene", "CC(C)c1ccccc1", 69.5, 6.14, "li_sunderland2012", False, "Experimental smoke point"),
    ("n-butylbenzene", "CCCCc1ccccc1", 64.3, 7.39, "li_sunderland2012", False, "Experimental smoke point"),
    ("sec-butylbenzene", "CCC(C)c1ccccc1", 68.4, 6.97, "li_sunderland2012", False, "Experimental smoke point"),
    ("isobutylbenzene", "CC(C)Cc1ccccc1", 64.1, 7.41, "li_sunderland2012", False, "Experimental smoke point"),
]

# Knots for piecewise monotonic calibration from unified YSI to TSI
YSI_KNOTS = np.array([30.4, 36.0, 42.6, 53.6, 61.7, 100.3, 170.9, 310.9, 336.0, 466.1, 649.1])
TSI_KNOTS = np.array([2.5, 2.8, 3.6, 5.0, 7.7, 27.0, 38.5, 58.9, 70.0, 88.0, 100.0])

def ysi_to_tsi(ysi: float) -> float:
    """Convert unified YSI to TSI scale using literature benchmark points."""
    return float(np.interp(ysi, YSI_KNOTS, TSI_KNOTS,
                           left=2.5 + (ysi - 30.4) * 0.05,
                           right=100.0 + (ysi - 649.1) * 0.065))

def tsi_to_sp(tsi: float, mw: float) -> float:
    """Invert Boehm 2022 apparatus relation to get estimated smoke point (mm)."""
    # TSI = a + b * (MW / SP) -> SP = b * MW / (TSI - a)
    if tsi - BOEHM_A <= 0:
        return 50.0
    return round(BOEHM_B * mw / (tsi - BOEHM_A), 2)

def main():
    os.makedirs('data/processed', exist_ok=True)
    comp_df = pd.read_csv(COMPOUNDS_CSV)
    hc_comp = comp_df[comp_df['is_hc'] == True].copy()
    comp_map = {r['inchikey']: r for _, r in hc_comp.iterrows()}

    ysi_df = pd.read_csv(YSI_CSV)
    ysi_hc = ysi_df[ysi_df['is_hc'] == True].copy()
    # Group by inchikey to take median if duplicates exist
    ysi_map = {}
    for ik, grp in ysi_hc.groupby('inchikey'):
        med_ysi = grp['ysi'].median()
        name = grp['name'].iloc[0]
        smi = grp['smiles'].iloc[0]
        orig_tier = grp['tier'].iloc[0]
        ysi_map[ik] = {
            'ysi': med_ysi,
            'name': name,
            'smiles': smi,
            'tier': orig_tier
        }

    # Map of Tier A items by InChIKey
    tier_a_by_ik = {}
    for name, smi, tsi, sp, src, cens, note in TIER_A_MEASURED:
        mol = Chem.MolFromSmiles(smi)
        ik = Chem.MolToInchiKey(mol) if mol else None
        if ik:
            tier_a_by_ik[ik] = {
                'compound_name': name,
                'smiles': Chem.MolToSmiles(mol),
                'tsi': tsi,
                'smoke_point_mm': sp,
                'source': src,
                'tier': 'A',
                'derived': False,
                'censored': cens,
                'note': note
            }

    rows = []
    # Iterate through all pure hydrocarbons in compounds.csv
    for ik, comp_row in comp_map.items():
        # Check if in Tier A
        if ik in tier_a_by_ik:
            item = tier_a_by_ik[ik]
            ysi_val = ysi_map[ik]['ysi'] if ik in ysi_map else None
            rows.append({
                'inchikey': ik,
                'smiles': comp_row['smiles'],
                'compound_name': item['compound_name'],
                'tsi': item['tsi'],
                'ysi': ysi_val,
                'smoke_point_mm': item['smoke_point_mm'],
                'source': item['source'],
                'tier': item['tier'],
                'derived': item['derived'],
                'censored': item['censored'],
                'note': item['note']
            })
        elif ik in ysi_map:
            yitem = ysi_map[ik]
            ysi_val = round(float(yitem['ysi']), 1)
            tsi_val = round(ysi_to_tsi(ysi_val), 2)
            mw_val = float(comp_row['mw']) if pd.notna(comp_row['mw']) else None
            sp_val = tsi_to_sp(tsi_val, mw_val) if mw_val else None
            cens = sp_val is not None and sp_val >= 50.0
            if cens:
                sp_val = 50.0
            
            tier_val = 'B' if yitem['tier'] == 'A' else 'C'
            rows.append({
                'inchikey': ik,
                'smiles': comp_row['smiles'],
                'compound_name': yitem['name'],
                'tsi': tsi_val,
                'ysi': ysi_val,
                'smoke_point_mm': sp_val,
                'source': 'ysi_unified',
                'tier': tier_val,
                'derived': True,
                'censored': cens,
                'note': f"Converted from YSI={ysi_val}" + ("; capped at 50 mm" if cens else "")
            })

    out_df = pd.DataFrame(rows)
    # Sort for reproducibility: Tier A first, then by TSI ascending
    out_df.sort_values(by=['tier', 'tsi'], inplace=True)
    out_df.to_csv(OUT_CSV, index=False)
    print(f"Compiled {len(out_df)} records to {OUT_CSV}")
    print(f"Tier counts: {out_df['tier'].value_counts().to_dict()}")
    print(f"Censored count (>50 mm): {out_df['censored'].sum()}")

    # Output constants JSON for baselines_smoke.py
    constants = {
        "boehm2022": {
            "a": BOEHM_A,
            "b": BOEHM_B,
            "units": "b in mm*mol/g, a dimensionless",
            "equation": "TSI = a + b * (MW / SP)",
            "inv_equation": "SP = b * MW / (TSI - a)",
            "anchors": {
                "methylcyclohexane": {"tsi": 5.0, "formula": "C7H14", "sp_mm": 43.56},
                "1-methylnaphthalene": {"tsi": 100.0, "formula": "C11H10", "sp_mm": 5.13}
            },
            "referee_fuels": {
                "1,3,5-trimethylbenzene": {"tsi_min": 54.7, "tsi_max": 63.1, "nominal_tsi": 58.9},
                "toluene_isooctane_40_60_vol": {"tsi_min": 21.0, "tsi_max": 25.4, "nominal_tsi": 23.2}
            }
        },
        "mensch2010": {
            "a": MENSCH_A,
            "b": MENSCH_B,
            "units": "b in mm*mol/g, a dimensionless",
            "equation": "TSI = a + b * (MW / SP)",
            "inv_equation": "SP = b * MW / (TSI - a)"
        },
        "censoring": {
            "d1322_lamp_max_mm": 50.0,
            "flag": "censored"
        }
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(constants, f, indent=2)
    print(f"Wrote configuration to {OUT_JSON}")

if __name__ == '__main__':
    main()
