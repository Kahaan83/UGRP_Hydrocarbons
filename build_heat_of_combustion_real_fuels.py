"""T27: Extract heat of combustion measurements and hydrogen content from AJFTD fuels.
Converts Net Heat of Combustion (NCV) to Gross Calorific Value (GCV):
GCV = NCV + 0.2122 * w_H (with w_H in wt% H).
Uses measured w_H where available (ASTM D3343, D3701, D7171, GCxGC);
otherwise falls back to default jet fuel w_H = 13.8 wt% (derived_w_H=True).

Output saved to data/processed/heat_of_combustion_real_fuels.csv.
"""
from __future__ import annotations
import os, re
import pandas as pd
from bs4 import BeautifulSoup

DEFAULT_W_H = 13.8

def clean_num(s: str | None) -> float | None:
    if s is None:
        return None
    s = str(s).strip()
    m = re.match(r'^([0-9]+(?:\.[0-9]+)?)', s)
    if not m:
        return None
    val = float(m.group(1))
    # Correct scale artifacts if present:
    # 1) AJFTD bug where Btu/lb was divided by 429.9226 instead of 0.002326 (yielding ~0.10 MJ/kg)
    if 0.08 <= val <= 0.12:
        val = val * 429.9226
    # 2) Raw kJ/kg or unscaled ~43000
    elif 35000 <= val <= 55000:
        val = val / 1000.0
    return val

def extract_ajftd_combustion(debug_dir: str = 'debug') -> pd.DataFrame:
    files = [f for f in os.listdir(debug_dir) if f.startswith('detail_') and f.endswith('.html')]
    re_btn = re.compile(rb'<button[^>]*onclick=["\']myFunction\(\s*["\']([^"\']+)["\']\s*\)["\'][^>]*>(.*?)</button>', re.DOTALL | re.IGNORECASE)
    re_label = re.compile(rb'<h[1-4][^>]*>\s*Fuel Details for\s*([^<]+)</h[1-4]>', re.IGNORECASE)

    raw_heat = []
    raw_h = []

    for fn in files:
        fid = fn.replace('detail_', '').replace('.html', '')
        path = os.path.join(debug_dir, fn)
        with open(path, 'rb') as fp:
            content = fp.read()
        
        has_heat = b'net_heat_of_combustion' in content
        has_h = b'hydrogen_content' in content
        if not (has_heat or has_h):
            continue
            
        m_lbl = re_label.search(content)
        label = m_lbl.group(1).decode('utf-8', errors='ignore').strip() if m_lbl else f'fuel_{fid}'
        
        for m in re_btn.finditer(content):
            tid = m.group(1).decode('utf-8', errors='ignore')
            btn_text = re.sub(rb'<[^>]+>', b'', m.group(2)).decode('utf-8', errors='ignore').strip().lower()
            if 'net_heat_of_combustion' in btn_text:
                tbl_pat = re.compile(rf'<table[^>]*id=["\']{re.escape(tid)}["\'][^>]*>(.*?)</table>'.encode('utf-8'), re.DOTALL | re.IGNORECASE)
                tbl_m = tbl_pat.search(content)
                if tbl_m:
                    soup = BeautifulSoup(b'<table>' + tbl_m.group(1) + b'</table>', 'html.parser')
                    trs = soup.find_all('tr')
                    if len(trs) > 1:
                        headers = [c.get_text(strip=True) for c in trs[0].find_all(['th', 'td'])]
                        for tr in trs[1:]:
                            vals = [c.get_text(strip=True) for c in tr.find_all(['th', 'td'])]
                            d = dict(zip(headers, vals))
                            d['fuel_id'] = fid
                            d['fuel_label'] = label
                            raw_heat.append(d)
            elif 'hydrogen_content' in btn_text:
                tbl_pat = re.compile(rf'<table[^>]*id=["\']{re.escape(tid)}["\'][^>]*>(.*?)</table>'.encode('utf-8'), re.DOTALL | re.IGNORECASE)
                tbl_m = tbl_pat.search(content)
                if tbl_m:
                    soup = BeautifulSoup(b'<table>' + tbl_m.group(1) + b'</table>', 'html.parser')
                    trs = soup.find_all('tr')
                    if len(trs) > 1:
                        headers = [c.get_text(strip=True) for c in trs[0].find_all(['th', 'td'])]
                        for tr in trs[1:]:
                            vals = [c.get_text(strip=True) for c in tr.find_all(['th', 'td'])]
                            d = dict(zip(headers, vals))
                            d['fuel_id'] = fid
                            d['fuel_label'] = label
                            raw_h.append(d)

    # Process hydrogen content by fuel_id
    h_by_fuel: dict[str, dict] = {}
    for r in raw_h:
        fid = r['fuel_id']
        val = clean_num(r.get('Value'))
        if val is not None and 8.0 <= val <= 25.0:
            meth = str(r.get('Test Method', '')).strip()
            entry = h_by_fuel.setdefault(fid, {'vals': [], 'methods': []})
            entry['vals'].append(val)
            if meth and meth.lower() not in ['nan', 'none', '']:
                entry['methods'].append(meth)

    # Process net heat of combustion by fuel_id
    heat_by_fuel: dict[str, dict] = {}
    for r in raw_heat:
        fid = r['fuel_id']
        val = clean_num(r.get('Value'))
        # Fuels have net heat of combustion typically between 35 and 55 MJ/kg
        if val is not None and 35.0 <= val <= 55.0:
            meth = str(r.get('Test Method', '')).strip()
            dt = str(r.get('Test Date', '')).strip()
            entry = heat_by_fuel.setdefault(fid, {
                'fuel_label': r['fuel_label'],
                'vals': [],
                'methods': [],
                'dates': []
            })
            entry['vals'].append(val)
            if meth and meth.lower() not in ['nan', 'none', '']:
                entry['methods'].append(meth)
            if dt and dt.lower() not in ['nan', 'none', 'not found', '']:
                entry['dates'].append(dt)

    records = []
    # Sort fuels naturally by integer id if possible
    sorted_fids = sorted(heat_by_fuel.keys(), key=lambda x: int(x) if x.isdigit() else x)
    for fid in sorted_fids:
        hdata = heat_by_fuel[fid]
        ncv = float(pd.Series(hdata['vals']).median())
        meth_str = '; '.join(sorted(set(hdata['methods']))) if hdata['methods'] else 'unspecified'
        dt_str = hdata['dates'][0] if hdata['dates'] else ''
        
        if fid in h_by_fuel:
            h_info = h_by_fuel[fid]
            w_H = float(pd.Series(h_info['vals']).median())
            derived_w_H = False
            h_meth_str = '; '.join(sorted(set(h_info['methods']))) if h_info['methods'] else 'unspecified'
        else:
            w_H = DEFAULT_W_H
            derived_w_H = True
            h_meth_str = 'default_jet_fuel_fallback'

        # Conversion formula from prompt: GCV = NCV + 0.2122 * w_H
        gcv = ncv + 0.2122 * w_H

        records.append({
            'fuel_id': fid,
            'fuel_label': hdata['fuel_label'],
            'ncv_MJ_kg': round(ncv, 4),
            'gcv_MJ_kg': round(gcv, 4),
            'w_H_wt_pct': round(w_H, 3),
            'derived_w_H': derived_w_H,
            'derived_gcv': True,  # GCV is calculated from NCV + w_H
            'n_ncv_obs': len(hdata['vals']),
            'ncv_method': meth_str,
            'w_H_method': h_meth_str,
            'test_date': dt_str,
            'source': 'AJFTD'
        })

    df_out = pd.DataFrame(records)
    return df_out

if __name__ == '__main__':
    out_csv = 'data/processed/heat_of_combustion_real_fuels.csv'
    df = extract_ajftd_combustion()
    df.to_csv(out_csv, index=False)
    print(f'Wrote {len(df)} real fuel calorific records to {out_csv}')
    print('Sample rows:')
    print(df.head())
    print('Summary:')
    print(df[['ncv_MJ_kg', 'gcv_MJ_kg', 'w_H_wt_pct']].describe())
    print('Hydrogen measurement breakdown:')
    print(df['derived_w_H'].value_counts())
