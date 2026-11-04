"""Generate data/processed/freezing_point_pure_hc.csv (Task T17).
Pure hydrocarbon freezing / melting points compiled from authoritative sources:
1. CRC Handbook of Chemistry and Physics (Organics table) - Tier A
2. NIST WebBook (Triple point / Melting point) - Tier A
3. Bradley Open Notebook Melting Point Dataset (Open_NTBKM) - Tier B
"""
import math, re, csv
import pandas as pd
import chemicals.phase_change as cpc
import chemicals.identifiers as ci
from chemicals import search_chemical

def is_hc(f):
    return bool(isinstance(f, str) and re.fullmatch(r'C(\d*)H(\d*)', f))

def main():
    cpc._load_phase_change_constants()
    
    pref_order = [
        ('CRC_ORG', 'CRC_Handbook_Organics', 'A'),
        ('WEBBOOK', 'NIST_WebBook', 'A'),
        ('OPEN_NTBKM', 'Open_Notebook_Melting_Points', 'B'),
        ('COMMON_CHEMISTRY', 'CAS_Common_Chemistry', 'B'),
        ('WIKIDATA', 'Wikidata', 'B')
    ]

    all_cas = set()
    for sname, _, _ in pref_order:
        all_cas |= set(str(c) for c in cpc.Tm_sources[sname].index)

    rows = []
    for cas in sorted(all_cas):
        try:
            s = search_chemical(cas)
            if not is_hc(s.formula):
                continue
            
            val_K = None
            chosen_source = None
            chosen_tier = None
            
            for sname, label, tier in pref_order:
                df = cpc.Tm_sources[sname]
                if cas in df.index:
                    val = df.loc[cas, 'Tm'] if 'Tm' in df.columns else None
                    if (val is None or (isinstance(val, float) and math.isnan(val))) and sname == 'WEBBOOK' and 'Tt' in df.columns:
                        val = df.loc[cas, 'Tt']
                    if isinstance(val, (int, float)) and not math.isnan(val):
                        val_K = float(val)
                        chosen_source = label
                        chosen_tier = tier
                        break
            
            if val_K is not None:
                val_C = round(val_K - 273.15, 2)
                rows.append({
                    'cas': cas,
                    'compound_name': s.common_name,
                    'formula': s.formula,
                    'smiles': s.smiles,
                    'inchikey': s.InChI_key,
                    'freezing_point_C': val_C,
                    'freezing_point_K': round(val_K, 2),
                    'unit': 'degC',
                    'source': chosen_source,
                    'tier': chosen_tier,
                    'derived': False
                })
        except Exception:
            pass

    out_path = 'data/processed/freezing_point_pure_hc.csv'
    fieldnames = ['cas', 'compound_name', 'formula', 'smiles', 'inchikey', 'freezing_point_C', 'freezing_point_K', 'unit', 'source', 'tier', 'derived']
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f'Wrote {len(rows)} pure hydrocarbon freezing point rows to {out_path}')

if __name__ == '__main__':
    main()
