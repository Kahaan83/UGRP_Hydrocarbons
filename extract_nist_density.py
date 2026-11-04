"""
T34: Extract Outcalt, Laesecke & Freund compressed liquid density measurements
for Jet-A (3 samples: Jet-A-3602, Jet-A-3638, Jet-A-4658), JP-8-3773, and S-8.
Extracted from NISTIR 6659 (Tables 20-24).
Saves to data/processed/nist_density_compressed_liquid.csv.
"""

import os
import pandas as pd
import numpy as np

NIST_TXT = 'data/raw/gc_composition/nistir6659.txt'
OUTPUT_CSV = 'data/processed/nist_density_compressed_liquid.csv'

with open(NIST_TXT, 'r', encoding='utf-8', errors='ignore') as f:
    text = f.read()

sub = text[text.find('Table 20:'):text.find('Table 25:')]
pages = sub.split('\x0c')

all_records = []

# --- Table 20: Jet-A-3602 ---
p0_lines = [l.strip() for l in pages[0].splitlines() if l.strip() and l.strip() not in ('58', '59')]
T20_1 = [float(x) for x in p0_lines[3:34]]
P20_1 = [float(x) for x in p0_lines[36:67]]
D20_1 = [float(x) for x in p0_lines[69:100]]

T20_2 = [float(x) for x in p0_lines[102:133]]
P20_2 = [float(x) for x in p0_lines[135:166]]
D20_2 = [float(x) for x in p0_lines[168:199]]

p1_lines = [l.strip() for l in pages[1].splitlines() if l.strip()][:219]
T20_3 = [float(x) for x in p1_lines[0:37]]
P20_3 = [float(x) for x in p1_lines[37:74]]
D20_3 = [float(x) for x in p1_lines[74:111]]

T20_4 = [float(x) for x in p1_lines[111:147]]
P20_4 = [float(x) for x in p1_lines[147:183]]
D20_4 = [float(x) for x in p1_lines[183:219]]

T20 = T20_1 + T20_2 + T20_3 + T20_4
P20 = P20_1 + P20_2 + P20_3 + P20_4
D20 = D20_1 + D20_2 + D20_3 + D20_4

for t, p, d in zip(T20, P20, D20):
    all_records.append({
        'fuel_id': 'Jet-A-3602',
        'fuel_name': 'Jet A (POSF 3602)',
        'category': 'conventional_jet',
        'T_K': t,
        'P_MPa': p,
        'P_kPa': p * 1000.0,
        'density_kg_m3': d,
        'uncertainty_kg_m3': 1.0,
        'source_table': 'Table 20',
        'reference': 'Outcalt, Laesecke & Freund, Energy & Fuels 2009 / NISTIR 6659'
    })

# --- Table 21: Jet-A-3638 ---
t21_p1_T = [270.00, 270.00, 270.00, 370.00, 370.00, 370.00]
t21_p1_P = [30.045, 25.104, 20.148, 5.088, 4.066, 3.048]
t21_p1_D = [827.22, 824.61, 821.94, 740.92, 739.86, 738.73]

p2_lines = [l.strip() for l in pages[2].splitlines() if l.strip() and l != '61']
left_T21 = [float(x) for x in p2_lines[0:43]]
left_P21 = [float(x) for x in p2_lines[43:86]]
left_D21 = [float(x) for x in p2_lines[129:172]]

right_T21 = [float(x) for x in p2_lines[86:129]]
right_P21 = [float(x) for x in p2_lines[172:215]]
right_D21 = [float(x) for x in p2_lines[215:258]]

p3_lines = [l.strip() for l in pages[3].splitlines() if l.strip()][:36]
p3_T21 = [float(x) for x in p3_lines[0:9]]
p3_P21 = [float(x) for x in p3_lines[9:18]]
p3_D21 = [float(x) for x in p3_lines[18:27]]

t_470_21 = [float(x) for x in p3_lines[27:36]]
p_470_21 = [float(x) for x in ['25.002', '20.051', '15.078', '10.062', '4.997', '3.981', '2.957', '1.941', '0.922']]
d_470_21 = [float(x) for x in ['694.28', '687.75', '680.53', '672.39', '662.99', '660.92', '658.76', '656.54', '654.23']]

T21 = [270.0, 270.0, 270.0] + left_T21 + [370.0, 370.0, 370.0] + p3_T21 + right_T21 + t_470_21
P21 = [30.045, 25.104, 20.148] + left_P21 + [5.088, 4.066, 3.048] + p3_P21 + right_P21 + p_470_21
D21 = [827.22, 824.61, 821.94] + left_D21 + [740.92, 739.86, 738.73] + p3_D21 + right_D21 + d_470_21

for t, p, d in zip(T21, P21, D21):
    all_records.append({
        'fuel_id': 'Jet-A-3638',
        'fuel_name': 'Jet A (POSF 3638)',
        'category': 'conventional_jet',
        'T_K': t,
        'P_MPa': p,
        'P_kPa': p * 1000.0,
        'density_kg_m3': d,
        'uncertainty_kg_m3': 1.0,
        'source_table': 'Table 21',
        'reference': 'Outcalt, Laesecke & Freund, Energy & Fuels 2009 / NISTIR 6659'
    })

# --- Table 22: Jet-A-4658 ---
p3_all = [l.strip() for l in pages[3].splitlines() if l.strip() and l != '62']
T_p3_1 = [float(x) for x in p3_all[39:70]]
P_p3_1 = [float(x) for x in p3_all[72:103]]
D_p3_1 = [float(x) for x in p3_all[105:136]]

p4_lines = [l.strip() for l in pages[4].splitlines() if l.strip()][:119]
T_block1 = [float(x) for x in p4_lines[0:30]]
P_block1 = [float(x) for x in p4_lines[30:60]]
D_block1 = [float(x) for x in p4_lines[60:90]]

T_p3_2 = [float(x) for x in p3_all[138:169]]
P_p3_2 = [float(x) for x in p3_all[180:211]]
D_p3_2 = [float(x) for x in p3_all[222:253]]

T_29 = [float(x) for x in p4_lines[90:119]]
lines4 = [l.strip() for l in pages[4].splitlines() if l.strip()][119:]
idx_63 = lines4.index('63')
P_29 = [float(x) for x in lines4[idx_63+1:idx_63+30]]
D_29 = [float(x) for x in lines4[idx_63+42:idx_63+71]]

T22 = T_p3_1 + T_block1 + T_p3_2 + T_29
P22 = P_p3_1 + P_block1 + P_p3_2 + P_29
D22 = D_p3_1 + D_block1 + D_p3_2 + D_29

for t, p, d in zip(T22, P22, D22):
    all_records.append({
        'fuel_id': 'Jet-A-4658',
        'fuel_name': 'Jet A (POSF 4658)',
        'category': 'conventional_jet',
        'T_K': t,
        'P_MPa': p,
        'P_kPa': p * 1000.0,
        'density_kg_m3': d,
        'uncertainty_kg_m3': 1.0,
        'source_table': 'Table 22',
        'reference': 'Outcalt, Laesecke & Freund, Energy & Fuels 2009 / NISTIR 6659'
    })

# --- Table 23: JP-8-3773 ---
T_t23_1 = [float(x) for x in lines4[3:13]]
P_t23_1 = [float(x) for x in lines4[15:25]]
D_t23_1 = [float(x) for x in lines4[27:37]]
idx_P2 = lines4.index('Pressure', 49)
idx_D2 = lines4.index('Density', idx_P2)
T_t23_2 = [float(x) for x in lines4[39:49]]
P_t23_2 = [float(x) for x in lines4[idx_P2+2:idx_P2+12]]
D_t23_2 = [float(x) for x in lines4[idx_D2+2:idx_D2+12]]

lines5 = [l.strip() for l in pages[5].splitlines() if l.strip() and l != '64']
left_T23 = [float(x) for x in lines5[0:46]]
left_P23 = [float(x) for x in lines5[46:92]]
left_D23 = [float(x) for x in lines5[138:184]]
right_T23 = [float(x) for x in lines5[92:138]]
right_P23 = [float(x) for x in lines5[184:230]]
right_D23 = [float(x) for x in lines5[230:276]]

lines6_t23 = [l.strip() for l in pages[6].splitlines() if l.strip()][:63]
T_16 = [float(x) for x in lines6_t23[0:16]]
P_16 = [float(x) for x in lines6_t23[16:32]]
D_16 = [float(x) for x in lines6_t23[32:48]]
T_15 = [float(x) for x in lines6_t23[48:63]]
lines6_t24 = [l.strip() for l in pages[6].splitlines() if l.strip()][63:]
idx_65 = lines6_t24.index('65')
P_15 = [float(x) for x in lines6_t24[idx_65+1:idx_65+16]]
D_15 = [float(x) for x in lines6_t24[-41:-26]]

T23_all = T_t23_1 + left_T23 + T_t23_2 + right_T23 + T_16 + T_15
P23_all = P_t23_1 + left_P23 + P_t23_2 + right_P23 + P_16 + P_15
D23_all = D_t23_1 + left_D23 + D_t23_2 + right_D23 + D_16 + D_15

for t, p, d in zip(T23_all, P23_all, D23_all):
    all_records.append({
        'fuel_id': 'JP-8-3773',
        'fuel_name': 'JP-8 (POSF 3773)',
        'category': 'military_jet',
        'T_K': t,
        'P_MPa': p,
        'P_kPa': p * 1000.0,
        'density_kg_m3': d,
        'uncertainty_kg_m3': 1.0,
        'source_table': 'Table 23',
        'reference': 'Outcalt, Laesecke & Freund, Energy & Fuels 2009 / NISTIR 6659'
    })

# --- Table 24: S-8 ---
pt_pairs = lines6_t24[5:53]
P_t24_1 = [float(pt_pairs[i]) for i in range(0, len(pt_pairs), 2)]
T_t24_1 = [float(pt_pairs[i+1]) for i in range(0, len(pt_pairs), 2)]
D_t24_1 = [float(x) for x in lines6_t24[55:55+24]]

T_t24_2 = [float(x) for x in lines6_t24[81:105]]
idx_p = lines6_t24.index('Pressure', 105)
P_t24_2 = [float(x) for x in lines6_t24[idx_p+2:idx_p+2+24]]
D_t24_2 = [float(x) for x in lines6_t24[-24:]]

lines7 = [l.strip() for l in pages[7].splitlines() if l.strip() and l != '66']
P_l = [float(x) for x in lines7[0:38]]
D_l = [float(x) for x in lines7[38:76]]
T_l = [float(x) for x in lines7[76:114]]

rem = lines7[114:]
T_r = [float(x) for x in rem if 250 <= float(x) <= 500]
P_r = [float(x) for x in rem if float(x) < 250]
D_r = [float(x) for x in rem if float(x) > 500]

T_s8 = T_t24_1 + T_t24_2 + T_l + T_r
P_s8 = P_t24_1 + P_t24_2 + P_l + P_r
D_s8 = D_t24_1 + D_t24_2 + D_l + D_r

for t, p, d in zip(T_s8, P_s8, D_s8):
    all_records.append({
        'fuel_id': 'S-8',
        'fuel_name': 'Synthetic Kerosene S-8 (POSF 4734)',
        'category': 'synthetic_kerosene',
        'T_K': t,
        'P_MPa': p,
        'P_kPa': p * 1000.0,
        'density_kg_m3': d,
        'uncertainty_kg_m3': 1.0,
        'source_table': 'Table 24',
        'reference': 'Outcalt, Laesecke & Freund, Energy & Fuels 2009 / NISTIR 6659'
    })

df_out = pd.DataFrame(all_records)
os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
df_out.to_csv(OUTPUT_CSV, index=False)

print(f"T34 Extracted {len(df_out)} compressed liquid density records across 5 fuels into {OUTPUT_CSV}")
print("Per-fuel row counts:")
print(df_out['fuel_id'].value_counts())
