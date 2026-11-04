# T0 — Data coverage audit

Data: `data/processed/model_table.csv` (42,662 rows, hydrocarbon-only), `data/unified/measurements.csv`
(101,710 raw rows incl. real fuels), `data/processed/excluded_by_R1.csv` (33,513 rows).

## 1. Model table: property × composition kind

| property | pure | binary | multi (≥3) |
|---|---|---|---|
| viscosity_dyn_mPas | 2,270 | 16,037 | 625 |
| density_kg_m3 | 1,752 | 12,747 | 897 |
| vapor_pressure_kPa | 4,415 | 2,021 | 64 |
| viscosity_kin_cSt | 67 | 28 | 0 |
| flash_point_C | 789 | 123 | 3 |
| mon | 267 | 30 | 76 |
| dcn | 80 | 213 | 133 |
| freezing_point_C | 25 | 0 | 0 |

## 2. Temperature coverage (@204 see measurements source split)

| property | T range (K) | n with T |
|---|---|---|
| viscosity_dyn_mPas | 230–453 | 18,929 |
| density_kg_m3 | 230–469 | 15,393 |
| vapor_pressure_kPa | 230–520 | 6,500 |
| viscosity_kin_cSt | 283–448 | 95 |
| flash point | single T per compound (method standard) | — |

## 3. Sources (top, per property)

- viscosity_dyn_mPas: nist_tde 15,775 · thermoml 12,051 · leenhouts 3,438 · chew2024 845
- density_kg_m3: thermoml 20,600 · ajftd 2,670 (real fuel, L0) · crc647 398 · av2718 95
- vapor_pressure_kPa: thermoml 8,713
- viscosity_kin_cSt: ajftd 2,481 (L0 real fuels) · thermoml 106
- flash_point_C: ajftd 2,413 (L0) · flashpoint_li 2,159 (pure) · leenhouts 366 · saldana11 327
- freezing_point_C: ajftd 2,672 (L0) · crc647 57 · misc 27 measured pure
- pour_point_C: 61 (CRC 647 aviation fuels, L0 real)

## 4. Real-fuel vs characterized composition (Rule R1)

| property | L3 (characterized, usable) | L0 (real fuel, no comp) |
|---|---|---|
| viscosity_dyn_mPas | 32,109 | 0 |
| flash_point_C | 2,852 | 2,487 |
| density_kg_m3 | 20,600 | 3,163 |
| viscosity_kin_cSt | 106 | 2,747 |
| freezing_point_C | 27 | 2,732 |
| pour_point_C | 0 | 61 |

## 5. GAP RANKING (rows < ~200 in model table = thin signal)

| rank | gap | rows now | impact on program |
|---|---|---|---|
| 1 | **Mixture flash point** (binary 123, multi 3) | 126 | CRITICAL for T5; safety metric is one-sided. Leenhouts binary 123 + 3 ternary is all we have. |
| 2 | **freezing_point_C pure** | 25 | T6 blocked. AJFTD has 2,672 freezing but L0 (no composition). |
| 3 | **viscosity_kin_cSt mixtures** | 95 | Low importance: kin = dyn/rho, derivable. |
| 4 | **pour_point_C** | 61 (L0 only) | T6 blocked; no characterized data exists in open sources. |
| 5 | **cloud_point_C** | 3 | T6 excluded. |
| 6 | **GCV/heat of combustion** | 0 pure HC | T6-1 blocked; only net_heat_content_btu_lb 57 real-fuel rows (crc647). |
| 7 | **octane (mon/dcn) mixture** | 76 (mon multi) / 133 (dcn multi) | marginal; not a core target property. |

## 6. excluded_by_R1.csv recoverability

**All 33,513 excluded rows have n_comp = 0** — none carry a component list. Recovery requires
finding composition for the labeled fuels (AJFTD / CRC aviation fuels are bulk spec fuels whose
composition was never published as molecule lists). **Recoverable only for:**
- AJFTD: if GC×GC or DHA (detailed hydrocarbon analysis) papers for specific fuels can be located.
- RP-1/JP surrogates: NISTIR 6646 has one 33-component RP-1 — already in L3.

**Do NOT expect bulk recovery from excluded_by_R1. Treat L0 real-fuel rows as a separate
bridge/stacking evaluation set (Rule R1 §9), not recoverable training data.**

## 7. Verdict for the task list

- T5 (mixture flash): data exists but only 123 binary + 3 ternary L3 rows — thin. Score on
  these plus keep locked real-fuel set separate.
- T6: GCV blocked (0 pure HC heats — find ISO/DIPPR-style tables or compute from group
  contribution), octane+DCN marginal, pour/cloud/freezing largely out of scope per rule.
- T7: highest-yield collection target = **more mixture flash point binaries with known
  composition** (e.g., Li et al. theoretical flash point mixtures, USNA mixture dataset).
- T2 wrapper: well-supported for viscosity + density at 230–470 K; flash limited to the 126
  mixture rows for residual fitting.

Done: coverage table exists; gaps ranked by impact; excluded rows assessed (n_comp=0, not
recoverable in bulk).
