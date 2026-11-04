"""Extract Autoignition Temperature (AIT) from chemicals.safety for compounds.csv."""
import pandas as pd
from chemicals import identifiers, safety

def main():
    safety._load_safety_data()
    db = identifiers.get_pubchem_db()
    
    # Root compounds for exact lookup
    root_comp = pd.read_csv('compounds.csv')
    root_cas_by_smi = dict(zip(root_comp['smiles'], root_comp['cas']))
    
    df_u = pd.read_csv('data/unified/compounds.csv')
    
    def get_cas(row):
        ik = row.get('inchikey')
        smi = row.get('smiles')
        inchi = row.get('inchi')
        if pd.notna(ik):
            m = db.search_InChI_key(ik)
            if m and m.CASs:
                return m.CASs
        if pd.notna(smi):
            m = db.search_smiles(smi)
            if m and m.CASs:
                return m.CASs
        if pd.notna(inchi):
            m = db.search_InChI(inchi)
            if m and m.CASs:
                return m.CASs
        if smi in root_cas_by_smi:
            return root_cas_by_smi[smi]
        return None

    hits = []
    for idx, row in df_u.iterrows():
        ik = row['inchikey']
        smi = row['smiles']
        cas = get_cas(row)
        if cas:
            t_k = safety.T_autoignition(cas)
            if t_k is not None:
                methods = safety.T_autoignition_methods(cas)
                hits.append({
                    'inchikey': ik,
                    'smiles': smi,
                    'cas': cas,
                    'T_autoignition_C': round(t_k - 273.15, 2),
                    'T_autoignition_K': round(t_k, 2),
                    'method': methods[0] if methods else 'chemicals.safety',
                    'tier': 'A',
                    'derived': False
                })

    df_hits = pd.DataFrame(hits)
    out_path = 'data/processed/ait_compounds.csv'
    df_hits.to_csv(out_path, index=False)
    print(f"Total checked: {len(df_u)}")
    print(f"Hits saved: {len(df_hits)} to {out_path}")
    print(f"Range: {df_hits['T_autoignition_C'].min()} to {df_hits['T_autoignition_C'].max()} deg C")

if __name__ == '__main__':
    main()
