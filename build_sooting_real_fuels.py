"""Generate data/processed/sooting_real_fuels.csv (Task T23).
Compile real-fuel sooting and smoke point metrics from AJFTD database and literature:
- 100 real fuels with ASTM D1322 smoke point (mm)
- Compositional covariates: total aromatics (% vol) and diaromatics (% vol)
- Standardized schema for sooting validation models
"""
import pandas as pd

def main():
    in_path = 'data/processed/smoke_point_ajftd.csv'
    df = pd.read_csv(in_path)
    
    # Standardize column naming and metadata for sooting baseline/validation
    soot_df = pd.DataFrame({
        'fuel_id': df['fuel_id'],
        'fuel_label': df['fuel_label'],
        'property': 'smoke_point_mm',
        'smoke_point_mm': df['value'],
        'unit': 'mm',
        'method': df['method'], # ASTM D1322
        'test_date': df['test_date'],
        'aromatics_vol_pct': df['aromatics_pct_vol'],
        'diaromatics_vol_pct': df['diaromatics_pct_vol'],
        'source': 'AJFTD_FAA_ASCENT',
        'tier': df['tier'],
        'derived': False
    })
    
    out_path = 'data/processed/sooting_real_fuels.csv'
    soot_df.to_csv(out_path, index=False, encoding='utf-8')
    print(f'Wrote {len(soot_df)} real fuel sooting rows to {out_path}')

if __name__ == '__main__':
    main()
