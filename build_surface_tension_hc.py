"""Generate data/processed/surface_tension_hc.csv (Task T20).
Pure hydrocarbon liquid surface tension parameters and 20 degC values compiled from:
1. Jasper (1972) / Lange Handbook linear temperature dependence parameters
2. Somayajulu (1988)
3. VDI PPDS-11 critical temperature scaling
4. Mulero & Cachadina (2012)
Evaluated at standard reference temperature T = 20 degC (293.15 K) in canonical units of mN/m.
For subcritical / liquid-phase hydrocarbons at 20 degC (293.15 K < Tc).
"""
import math, re, csv
import pandas as pd
import chemicals.interface as ci
from chemicals import search_chemical
from thermo import Chemical

def is_hc(f):
    return bool(isinstance(f, str) and re.fullmatch(r'C(\d*)H(\d*)', f))

def main():
    ci.load_interface_dfs()
    j = ci.sigma_data_Jasper_Lange
    s1 = ci.sigma_data_Somayajulu
    s2 = ci.sigma_data_Somayajulu2
    ppds = ci.sigma_data_VDI_PPDS_11
    mc = ci.sigma_data_Mulero_Cachadina

    all_cas = set()
    for d in [j, s1, s2, ppds, mc]:
        all_cas |= set(str(c) for c in d.index)

    rows = []
    for cas in sorted(all_cas):
        try:
            s = search_chemical(cas)
            if not is_hc(s.formula):
                continue
            
            c = Chemical(cas)
            st = c.SurfaceTension
            valid = st.valid_methods()
            
            val_20C_mN_m = None
            chosen_m = None
            
            # If 293.15 K is above critical temperature Tc, liquid surface tension is zero / supercritical
            if c.Tc and 293.15 >= c.Tc:
                val_20C_mN_m = 0.0
                chosen_m = 'supercritical_T_ge_Tc'
            else:
                # Prioritize standard experimental correlations
                for m in ['JASPER', 'SOMAYAJULU', 'SOMAYAJULU2', 'VDI_PPDS', 'REFPROP', 'REFPROP_FIT']:
                    if m in valid:
                        try:
                            v = st.calculate(293.15, m)
                            if v is not None and not math.isnan(v):
                                val_20C_mN_m = round(v * 1000.0, 3) # N/m to mN/m (dyn/cm)
                                chosen_m = m
                                break
                        except Exception:
                            pass
            
            # Record Jasper parameters if available
            jasper_a = None
            jasper_b = None
            Tmin_K = None
            Tmax_K = None
            if cas in j.index:
                row_j = j.loc[cas]
                jasper_a = float(row_j['a'])
                jasper_b = float(row_j['b'])
                Tmin_K = float(row_j['Tmin'])
                Tmax_K = float(row_j['Tmax'])
                
            available_sources = []
            if cas in j.index: available_sources.append('Jasper_Lange')
            if cas in s1.index: available_sources.append('Somayajulu')
            if cas in ppds.index: available_sources.append('VDI_PPDS_11')
            if cas in mc.index: available_sources.append('Mulero_Cachadina')
            
            is_liquid_20C = bool(c.Tc and 293.15 < c.Tc)
            
            rows.append({
                'cas': cas,
                'compound_name': s.common_name,
                'formula': s.formula,
                'smiles': s.smiles,
                'inchikey': s.InChI_key,
                'surface_tension_20C_mN_m': val_20C_mN_m if val_20C_mN_m is not None else '',
                'unit': 'mN/m',
                'T_ref_C': 20.0,
                'eval_method': chosen_m if chosen_m is not None else '',
                'is_liquid_phase_at_20C': is_liquid_20C,
                'Tc_K': round(float(c.Tc), 2) if c.Tc else '',
                'jasper_a': jasper_a if jasper_a is not None else '',
                'jasper_b': jasper_b if jasper_b is not None else '',
                'Tmin_K': Tmin_K if Tmin_K is not None else '',
                'Tmax_K': Tmax_K if Tmax_K is not None else '',
                'sources': '; '.join(available_sources),
                'tier': 'A',
                'derived': False
            })
        except Exception:
            pass

    out_path = 'data/processed/surface_tension_hc.csv'
    fieldnames = ['cas', 'compound_name', 'formula', 'smiles', 'inchikey', 'surface_tension_20C_mN_m', 'unit', 'T_ref_C', 'eval_method', 'is_liquid_phase_at_20C', 'Tc_K', 'jasper_a', 'jasper_b', 'Tmin_K', 'Tmax_K', 'sources', 'tier', 'derived']
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f'Wrote {len(rows)} pure hydrocarbon surface tension rows to {out_path}')

if __name__ == '__main__':
    main()
