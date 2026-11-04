"""T26: Ensure 100% GCV coverage for compounds in data/unified/compounds.csv.
Compounds missing experimental / stoichiometric GCV are filled via modified Dulong formula:
GCV (MJ/kg) = 33.8*w_C + 144.2*(w_H - w_O/8) + 9.4*w_S
Flag tier: C, derived: True.
"""
import re
import pandas as pd

aw = {'C': 12.011, 'H': 1.008, 'O': 15.999, 'S': 32.06, 'N': 14.007}

comp_path = 'data/unified/compounds.csv'
gcv_path = 'data/processed/gcv_MJ_kg.csv'

comp = pd.read_csv(comp_path)
gcv = pd.read_csv(gcv_path)

existing_keys = set(gcv['inchikey'])
print(f'Starting with {len(gcv)} GCV entries out of {len(comp)} compounds.')

new_rows = []
for _, r in comp.iterrows():
    ik = r['inchikey']
    if ik in existing_keys:
        continue
    smiles = r['smiles']
    f_str = r['formula']
    elems = re.findall(r'([A-Z][a-z]?)(\d*)', str(f_str))
    f_dict = {el: int(num) if num else 1 for el, num in elems}
    
    mw = sum(f_dict[el] * aw.get(el, 0.0) for el in f_dict)
    if mw <= 0:
        continue
    w_C = f_dict.get('C', 0) * aw['C'] / mw
    w_H = f_dict.get('H', 0) * aw['H'] / mw
    w_O = f_dict.get('O', 0) * aw['O'] / mw
    w_S = f_dict.get('S', 0) * aw['S'] / mw
    
    # Modified Dulong formula
    val = 33.8 * w_C + 144.2 * (w_H - w_O / 8.0) + 9.4 * w_S
    
    new_rows.append({
        'inchikey': ik,
        'smiles': smiles,
        'property': 'gcv_MJ_kg',
        'value': f'{val:.4f}',
        'unit': 'MJ/kg',
        'method': 'modified_Dulong',
        'tier': 'C',
        'derived': True,
        'T_ref_K': 298.15,
        'ref': 'elemental_composition_modified_Dulong'
    })

if new_rows:
    df_new = pd.DataFrame(new_rows)
    df_combined = pd.concat([gcv, df_new], ignore_index=True)
    df_combined.to_csv(gcv_path, index=False)
    print(f'Added {len(new_rows)} tier C fallback entries. Total now: {len(df_combined)}.')
else:
    print('All compounds already had GCV.')
