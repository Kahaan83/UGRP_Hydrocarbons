"""
Script to fetch Harvard Dataverse Yield Sooting Index (YSI) Volume 2 dataset
and build data/raw/sooting/ysi_compounds.csv cross-referenced with compounds.csv.
"""
import io
import os
import sys
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from common import is_hydrocarbon, hc_class

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
OUT_DIR = os.path.join(BASE_DIR, 'data', 'raw', 'sooting')
COMP_PATH = os.path.join(BASE_DIR, 'data', 'unified', 'compounds.csv')
OUT_CSV = os.path.join(OUT_DIR, 'ysi_compounds.csv')

os.makedirs(OUT_DIR, exist_ok=True)

# 1. Fetch raw files if not already present
s = requests.Session()
s.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})

files_to_download = [
    ('3048960', 'YSI_Database_Volume_2.tab', False),
    ('3048959', 'Detailed_YSI_Database_Volume_2.tab', False),
    ('3048961', 'Readme.pdf', True),
    ('3048974', 'ysi-primer.pdf', True),
]

for file_id, filename, is_binary in files_to_download:
    dest_path = os.path.join(OUT_DIR, filename)
    if not os.path.exists(dest_path):
        url = f'https://dataverse.harvard.edu/api/access/datafile/{file_id}'
        r = s.get(url, timeout=30, verify=False)
        r.raise_for_status()
        mode = 'wb' if is_binary else 'w'
        kwargs = {} if is_binary else {'encoding': 'utf-8'}
        with open(dest_path, mode, **kwargs) as f:
            if is_binary:
                f.write(r.content)
            else:
                f.write(r.text)

# 2. Parse Detailed YSI Database Volume 2
det_path = os.path.join(OUT_DIR, 'Detailed_YSI_Database_Volume_2.tab')
with open(det_path, 'r', encoding='utf-8') as f:
    text = f.read()

# Fix encoding mojibake for primes: \xc3\xa2\xc2\x80\xc2\x99 -> '
text = text.replace('\u00e2\u0080\u0099', "'").replace('\u2019', "'").replace('\u2032', "'")
df_det = pd.read_csv(io.StringIO(text), sep='\t')

# 3. Load compounds.csv registry
comp = pd.read_csv(COMP_PATH)
comp_inchikeys = set(comp['inchikey'].dropna())

# 4. Process molecules and cross-reference
def parse_mol(smi):
    try:
        return Chem.MolFromSmiles(smi)
    except Exception:
        return None

df_det['mol'] = df_det['SMILES'].apply(parse_mol)
df_det['inchikey'] = df_det['mol'].apply(lambda m: inchi.MolToInchiKey(m) if m else None)
df_det['smiles_can'] = df_det['mol'].apply(lambda m: Chem.MolToSmiles(m) if m else None)
df_det['is_hc'] = df_det['mol'].apply(lambda m: is_hydrocarbon(m) if m else False)
df_det['hc_class'] = df_det['mol'].apply(lambda m: hc_class(m) if (m and is_hydrocarbon(m)) else 'other')
df_det['in_unified_compounds'] = df_det['inchikey'].isin(comp_inchikeys)

ref_map = {
    1.0: 'Kuhn et al. 2011 (Proc. Combust. Inst. 33:743)',
    2.0: 'Das et al. 2017 (Proc. Combust. Inst. 36:871)',
    3.0: 'Das et al. 2018 (Combust. Flame 197:445 / McEnally & Pfefferle)',
    4.0: 'McEnally & Pfefferle (Proc. Combust. Inst. 32:673)',
}

def get_source(r):
    ref_num = r['Ref#']
    base = ref_map.get(ref_num, 'Harvard Dataverse doi:10.7910/DVN/7HGFT8')
    return f"{base}; Das et al. 2018 / McEnally & Pfefferle"

df_det['source'] = df_det.apply(get_source, axis=1)
df_det['tier'] = 'A'
df_det['scale'] = 'unified (benzene=100.3, toluene=170.9, n-heptane=36.0, YSI-molar)'

out_df = pd.DataFrame({
    'inchikey': df_det['inchikey'],
    'smiles': df_det['smiles_can'],
    'smiles_raw': df_det['SMILES'],
    'name': df_det['Species'],
    'formula': df_det['Formula'],
    'cas': df_det['CAS#'],
    'ysi': df_det['UnifiedYSI'],
    'ysi_unc': df_det['UnifiedYSIError'],
    'scale': df_det['scale'],
    'tier': df_det['tier'],
    'source': df_det['source'],
    'is_hc': df_det['is_hc'],
    'cls': df_det['hc_class'],
    'in_unified_compounds': df_det['in_unified_compounds'],
    'original_ysi': df_det['OriginalYSI'],
    'original_scale': df_det['OriginalScale'],
})

out_df.to_csv(OUT_CSV, index=False, encoding='utf-8')
print(f"Written {len(out_df)} rows to {OUT_CSV}")
print(f"Total compounds: {len(out_df)}")
print(f"Pure hydrocarbons: {out_df['is_hc'].sum()}")
print(f"Matched with unified compounds.csv: {out_df['in_unified_compounds'].sum()}")
print(f"Matched pure hydrocarbons: {(out_df['in_unified_compounds'] & out_df['is_hc']).sum()}")
