"""Generate data/processed/octane_ron_mon_hc.csv (Task T16).
Hydrocarbon RON/MON compiled from vom Lehn et al. 2021, CombustDB, and Leenhouts 2026 / Kuzhagaliyeva 2022.
Includes experimental data (Tier A) and ANN model predictions (Tier C).
"""
import math, re, csv
import pandas as pd
import chemicals.combustion as cc
import chemicals.identifiers as ci
from chemicals import search_chemical

def is_hc(f):
    return bool(isinstance(f, str) and re.fullmatch(r'C(\d*)H(\d*)', f))

def main():
    # 1. Load Leenhouts pure MON (Kuzhagaliyeva 2022 experimental)
    lh_path = 'data/interim/leenhouts_benchmark/datasets/datasets/fuel_ignition_numbers/processed_data/mon_pure_exp.csv'
    df_lh = pd.read_csv(lh_path)
    lh_map = {}
    for idx, row in df_lh.iterrows():
        inchi = row['component_inchi_0']
        val = float(row['value'])
        try:
            cas = ci.CAS_from_any(inchi)
            s = search_chemical(cas)
            if is_hc(s.formula):
                lh_map[cas] = (val, s)
        except Exception:
            pass

    # 2. Load chemicals combustion databases (vom Lehn 2021 & CombustDB)
    cc._load_combustion_data()
    fl_ron = cc.RON_sources['FLORIAN_LIMING']
    fl_mon = cc.MON_sources['FLORIAN_LIMING']
    fl_ann = cc.RON_sources['FLORIAN_LIMING_ANN']
    comb_ron = cc.RON_sources['COMBUSTDB']
    comb_mon = cc.MON_sources['COMBUSTDB']

    all_cas = set(str(c) for c in (set(fl_ron.index) | set(fl_mon.index) | set(comb_ron.index) | set(comb_mon.index) | set(fl_ann.index) | set(lh_map.keys())))

    rows = []
    for cas in sorted(all_cas):
        try:
            if cas in lh_map:
                s = lh_map[cas][1]
            else:
                s = search_chemical(cas)
            if not is_hc(s.formula):
                continue
                
            ron_exp = None
            mon_exp = None
            source_parts = []
            
            # Priority for RON: FLORIAN_LIMING (vom Lehn) > COMBUSTDB
            if cas in fl_ron.index and not math.isnan(fl_ron.loc[cas, 'RON']):
                ron_exp = float(fl_ron.loc[cas, 'RON'])
                source_parts.append('vom_Lehn_2021_exp')
            elif cas in comb_ron.index and not math.isnan(comb_ron.loc[cas, 'RON']):
                ron_exp = float(comb_ron.loc[cas, 'RON'])
                source_parts.append('CombustDB_exp')
                
            # Priority for MON: FLORIAN_LIMING > COMBUSTDB > Leenhouts/Kuzhagaliyeva
            if cas in fl_mon.index and not math.isnan(fl_mon.loc[cas, 'MON']):
                mon_exp = float(fl_mon.loc[cas, 'MON'])
                if 'vom_Lehn_2021_exp' not in source_parts:
                    source_parts.append('vom_Lehn_2021_exp')
            elif cas in comb_mon.index and not math.isnan(comb_mon.loc[cas, 'MON']):
                mon_exp = float(comb_mon.loc[cas, 'MON'])
                if 'CombustDB_exp' not in source_parts:
                    source_parts.append('CombustDB_exp')
            elif cas in lh_map:
                mon_exp = lh_map[cas][0]
                source_parts.append('Leenhouts_2026_Kuzhagaliyeva_2022')
                
            # ANN model fallback (vom Lehn 2021 ANN predictions)
            ron_ann = None
            mon_ann = None
            if cas in fl_ann.index:
                r = fl_ann.loc[cas, 'RON']
                m = fl_ann.loc[cas, 'MON']
                if not math.isnan(r): ron_ann = float(r)
                if not math.isnan(m): mon_ann = float(m)
                
            has_exp = (ron_exp is not None or mon_exp is not None)
            tier = 'A' if has_exp else 'C'
            derived = not has_exp
            
            ron_val = ron_exp if ron_exp is not None else ron_ann
            mon_val = mon_exp if mon_exp is not None else mon_ann
            
            source_final = '; '.join(source_parts) if source_parts else 'vom_Lehn_2021_ANN'
                
            sens = None
            if ron_val is not None and mon_val is not None:
                sens = round(ron_val - mon_val, 2)
                
            rows.append({
                'cas': cas,
                'compound_name': s.common_name,
                'formula': s.formula,
                'smiles': s.smiles,
                'inchikey': s.InChI_key,
                'ron': round(ron_val, 2) if ron_val is not None else '',
                'mon': round(mon_val, 2) if mon_val is not None else '',
                'octane_sensitivity': sens if sens is not None else '',
                'source': source_final,
                'tier': tier,
                'derived': derived,
                'is_experimental': has_exp
            })
        except Exception:
            pass

    out_path = 'data/processed/octane_ron_mon_hc.csv'
    fieldnames = ['cas', 'compound_name', 'formula', 'smiles', 'inchikey', 'ron', 'mon', 'octane_sensitivity', 'source', 'tier', 'derived', 'is_experimental']
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f'Wrote {len(rows)} pure hydrocarbon octane rows to {out_path}')

if __name__ == '__main__':
    main()
