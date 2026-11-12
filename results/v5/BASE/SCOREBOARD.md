# v5 BASE scoreboard

781 cells; 44 `no_evaluation`, 191 `n_too_small` (see `data_gaps.csv`).  Every T2 fuel has been read, so T2 cells are `descriptive`, not held-out.  Error is predicted minus measured; viscosity error is ln(pred/meas).  Bias, MAE and unsafe-side rate are in the property's unit (viscosity: ln).

## Readings

**R-BASE-1** (48 required (property, tier) cells): all exist with a status.

**R-BASE-2** (recomputed cells with a committed counterpart):

| cell | observed | committed | diff | tol | within | committed_source | cause |
|---|---|---|---|---|---|---|---|
| flash T2 AJFTD GCxGC, all 37 eligible rows: bias | -0.426 | -0.426 | -3.21e-08 | 0.05 | True | results/v5/G1/summary.csv |  |
| flash T2 AJFTD GCxGC, all 37 eligible rows: MAE | 3.54 | 3.54 | 1.08e-07 | 0.05 | True | results/v5/G1/summary.csv |  |
| flash T2 AJFTD GCxGC, headline one-per-group (31): bias | -0.634 | -0.426 | -0.208 | 0.05 | False | results/v5/G1/summary.csv | BASE scores one fuel per double-count group (31); G1 scored all 37 eligible rows (6 duplicate-group fuels); the 37-row rerun above matches G1 |
| ln viscosity AV-33-22 233.15 K: bias | -0.637 | -0.637 | -0.000286 | 0.005 | True | results/v5/C1 (R-C1-5, rounded to 3 decimals) |  |
| ln viscosity AV-33-22 313.15 K: bias | 0.0656 | 0.066 | -0.00036 | 0.005 | True | results/v5/C1 (R-C1-5, rounded to 3 decimals) |  |
| flash CRC 647 r3 dev (45): bias | 0.436 | 0.436 | 1.27e-09 | 0.05 | True | results/v5/L1/comparison.csv (before) |  |
| freezing T2 pooled F1: MAE | 3.46 | 3.46 | 0.00429 | 0.05 | True | results/v5/S5 (README, rounded) |  |

**R-BASE-3** (designated model vs the leave-one-group-out constant, aggregate rows without temperature split, by tier; a tier appears in several columns when its sets or models disagree, see the per-property tables):

| property | beats | no_skill | worse | too_small_only |
|---|---|---|---|---|
| density_kg_m3 | T2_measured | CRC647_r3 | - | diesel |
| flash_point_C | CRC647_r3, T2_measured | T2_measured, T3_ajftd | T2_measured | diesel |
| freezing_point_C | T3_ajftd | CRC647_r3, T2_measured, T3_ajftd | - | - |
| heat_of_combustion_net_MJ_kg | T3_ajftd | T2_measured, T3_ajftd | CRC647_r3, T3_ajftd | - |
| mon | - | - | T1_lab | - |
| pour_point_C | - | CRC647_r3 | - | - |
| ron | - | - | T1_lab | - |
| smoke_point_mm | T3_ajftd | T3_ajftd | - | - |
| viscosity_kin_cSt | CRC647_r3, T2_measured | T3_ajftd | T2_measured | diesel |

**R-BASE-4** (largest unsafe-side rate per property, designated rows only; flash = over-prediction > 5 C, cold viscosity = under-prediction beyond the S1 margin):

| property | max_unsafe_rate | tier | set | model | variant | label | n_groups |
|---|---|---|---|---|---|---|---|
| flash_point_C | 0.25 | CRC647_r3 | crc647_r3 | ideal_le_chatelier | designated | fold_b | 12 |
| viscosity_kin_cSt | 1 | T2_measured | ajftd_gcxgc | refutas_raw | designated | descriptive | 10 |

**R-BASE-5** data gaps: `data_gaps.csv`, 233 rows (44 no_evaluation, 191 n_too_small).

## Decisions waiting on the user (U1)

Designated (v4 routing) stays unchanged until the user decides.  Each row: the same cell scored with and without the change.

| decision | property | tier | set | T_K | model | variant | label | n_groups | bias_designated | bias_variant | mae_designated | mae_variant | unsafe_designated | unsafe_variant |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LFL overlay | flash_point_C | CRC647_r3 | crc647_r3 |  | ideal_le_chatelier | lfl_overlay | dev | 45 | 0.436 | -0.726 | 4.03 | 4.09 | 0.133 | 0.0444 |
| LFL overlay | flash_point_C | CRC647_r3 | crc647_r3 |  | ideal_le_chatelier | lfl_overlay | fold_b | 12 | 2.46 | 1.41 | 3.72 | 3.26 | 0.25 | 0.0833 |
| LFL overlay | flash_point_C | T1_lab | lab_mixtures_31_systems |  | ideal_le_chatelier | lfl_overlay | dev | 31 | -3.3 | -3.25 | 3.74 | 3.69 | 0.00794 | 0.00794 |
| LFL overlay | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | lfl_overlay | descriptive | 31 | -0.634 | -0.478 | 3.38 | 3.31 | 0.0968 | 0.0968 |
| LFL overlay | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 31 | -0.634 | -0.841 | 3.38 | 3.65 | 0.0968 | 0.0968 |
| LFL overlay | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | lfl_overlay | descriptive | 35 | -1.06 | -1.14 | 1.5 | 1.55 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 35 | -1.06 | -1.46 | 1.5 | 1.69 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | lfl_overlay | descriptive | 28 | -1.07 | -1.13 | 1.56 | 1.61 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 28 | -1.07 | -1.42 | 1.56 | 1.7 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | lfl_overlay | descriptive | 32 | -0.387 | -0.271 | 3.1 | 3.05 | 0.0938 | 0.0938 |
| LFL overlay | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 32 | -0.387 | -0.661 | 3.1 | 3.38 | 0.0938 | 0.0938 |
| LFL overlay | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | lfl_overlay | descriptive | 10 | -0.921 | -0.799 | 3.79 | 3.77 | 0.1 | 0.1 |
| LFL overlay | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 10 | -0.921 | -1.12 | 3.79 | 4.09 | 0.1 | 0.1 |
| LFL overlay | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | lfl_overlay | descriptive | 4 | 1.36 | 1.57 | 5.82 | 5.86 | 0.25 | 0.25 |
| LFL overlay | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 4 | 1.36 | 1.31 | 5.82 | 6.12 | 0.25 | 0.25 |
| LFL overlay | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | lfl_overlay | descriptive | 70 | -0.734 | -0.692 | 2.58 | 2.58 | 0.0571 | 0.0571 |
| LFL overlay | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 70 | -0.734 | -1.03 | 2.58 | 2.81 | 0.0571 | 0.0571 |
| Corrected mesitylene Antoine set | flash_point_C | T1_lab | lab_mixtures_31_systems |  | ideal_le_chatelier | mesitylene_fixed | dev | 31 | -3.3 | -3.59 | 3.74 | 3.89 | 0.00794 | 0.00794 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 31 | -0.634 | -0.841 | 3.38 | 3.65 | 0.0968 | 0.0968 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 31 | -0.634 | -0.994 | 3.38 | 3.72 | 0.0968 | 0.0968 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 35 | -1.06 | -1.46 | 1.5 | 1.69 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 35 | -1.06 | -1.39 | 1.5 | 1.62 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 28 | -1.07 | -1.42 | 1.56 | 1.7 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 28 | -1.07 | -1.35 | 1.56 | 1.65 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 32 | -0.387 | -0.661 | 3.1 | 3.38 | 0.0938 | 0.0938 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 32 | -0.387 | -0.775 | 3.1 | 3.41 | 0.0938 | 0.0938 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 10 | -0.921 | -1.12 | 3.79 | 4.09 | 0.1 | 0.1 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 10 | -0.921 | -1.24 | 3.79 | 4.11 | 0.1 | 0.1 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 4 | 1.36 | 1.31 | 5.82 | 6.12 | 0.25 | 0.25 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 4 | 1.36 | 1.1 | 5.82 | 6.07 | 0.25 | 0.25 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 70 | -0.734 | -1.03 | 2.58 | 2.81 | 0.0571 | 0.0571 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 70 | -0.734 | -1.07 | 2.58 | 2.81 | 0.0571 | 0.0571 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 10 | -0.698 | -0.0607 | 0.698 | 0.108 | 1 | 0.1 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 28 | -0.291 | -0.0269 | 0.344 | 0.124 | 0.179 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 273 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.0753 | -0.0199 | 0.0753 | 0.0199 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 288 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.00338 | -0.0225 | 0.00909 | 0.0225 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 293 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | 0.00233 | -0.0332 | 0.0206 | 0.0332 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 16 | -0.00696 | -0.0726 | 0.0541 | 0.0744 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc_T_263_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.0753 | -0.0199 | 0.0753 | 0.0199 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 17 | -0.018 | -0.0746 | 0.0441 | 0.0763 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 28 | -0.405 | -0.0209 | 0.458 | 0.134 | 0.464 | 0.0357 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.637 | 0.0607 | 0.637 | 0.116 | 0.914 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.264 | 0.0269 | 0.264 | 0.0624 | 0.0286 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 273 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.0554 | 0.0199 | 0.062 | 0.042 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 288 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.0191 | 0.0225 | 0.0309 | 0.0326 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 293 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.0355 | 0.0332 | 0.0401 | 0.0384 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.0656 | 0.0726 | 0.0656 | 0.0726 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322_T_263_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.0554 | 0.0199 | 0.062 | 0.042 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.0191 | 0.0225 | 0.0309 | 0.0326 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.637 | 0.0607 | 0.637 | 0.116 | 0.914 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | -0.664 | -0.0272 | 0.664 | 0.125 | 1 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | -0.312 | -0.048 | 0.312 | 0.0844 | 0.25 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | 0.0325 | -0.0331 | 0.0325 | 0.0331 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 373 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.0197 | -0.0127 | 0.0197 | 0.0127 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | 0.0319 | -0.0332 | 0.0319 | 0.0332 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | -0.664 | -0.0272 | 0.664 | 0.125 | 1 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 48 | -0.65 | 0.0312 | 0.65 | 0.115 | 0.938 | 0.0208 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 65 | -0.275 | 0.00332 | 0.298 | 0.0888 | 0.0923 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 273 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 38 | -0.057 | 0.0167 | 0.063 | 0.0403 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 288 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 38 | 0.0173 | 0.019 | 0.0292 | 0.0318 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 293 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 39 | 0.0321 | 0.0264 | 0.0381 | 0.0379 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 54 | 0.042 | 0.0234 | 0.0601 | 0.0712 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 373 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.008 | -0.0127 | 0.0162 | 0.0127 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled_T_263_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 38 | -0.057 | 0.0167 | 0.063 | 0.0403 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 53 | 0.00976 | -0.00893 | 0.0356 | 0.0459 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 65 | -0.547 | 0.0254 | 0.57 | 0.125 | 0.738 | 0.0154 |
| SLE2 two-compilations tier | freezing_point_C | T1_lab | sle_alkane_blends |  | ideal_svl | sle2_single_compilation | dev | 38 | -1.1 | -4.02 | 2.02 | 7.03 |  |  |
| SLE2 two-compilations tier | freezing_point_C | T1_lab | sle_alkane_blends |  | ideal_svl | sle2_two_compilations | dev | 11 | -1.1 | -1.41 | 2.02 | 2.13 |  |  |

## viscosity_kin_cSt

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_r3 | Recalibrated Refutas | designated | dev |  | 45 | -6.55e-16 | 0.0873 | 0 | beats_trivial | 0.347 | 0.978 |
| CRC647_r3 | crc647_r3 | Recalibrated Refutas | designated | fold_b |  | 12 | 0.0248 | 0.109 | 0 | beats_trivial | 0.347 | 0.917 |
| T1_lab | T_263_283K | refutas_raw | designated | already_read |  | 4 | 0.0382 | 0.0718 | 0 | not_applicable | 0.608 | 1 |
| T1_lab | T_263_283K | refutas_raw | designated | dev |  | 16 | -0.0233 | 0.0677 | 0 | not_applicable | 0.608 | 1 |
| T1_lab | T_ge_283K | refutas_raw | designated | already_read |  | 70 | 0.0415 | 0.0533 |  | not_applicable | 0.186 | 0.943 |
| T1_lab | T_ge_283K | refutas_raw | designated | dev |  | 264 | 0.0245 | 0.068 |  | not_applicable | 0.186 | 0.909 |
| T1_lab | T_lt_263K | refutas_raw | designated | already_read |  | 3 | -0.079 | 0.0973 | 0 | not_applicable | 0.415 | 1 |
| T1_lab | T_lt_263K | refutas_raw | designated | dev |  | 20 | -0.113 | 0.118 | 0.0354 | not_applicable | 0.415 | 0.9 |
| T2_measured | ajftd_gcxgc_T_263_283K | refutas_raw | designated | descriptive |  | 3 | -0.0753 | 0.0753 | 0 | n_too_small | inf |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | refutas_raw | designated | descriptive |  | 17 | -0.018 | 0.0441 |  | beats_trivial | 0.186 | 0.941 |
| T2_measured | ajftd_gcxgc_T_lt_263K | refutas_raw | designated | descriptive |  | 28 | -0.405 | 0.458 | 0.464 | worse_than_trivial | inf |  |
| T2_measured | av3322_T_263_283K | refutas_raw | designated | descriptive |  | 35 | -0.0554 | 0.062 | 0 | beats_trivial | inf |  |
| T2_measured | av3322_T_ge_283K | refutas_raw | designated | descriptive |  | 35 | 0.0191 | 0.0309 |  | beats_trivial | 0.186 | 1 |
| T2_measured | av3322_T_lt_263K | refutas_raw | designated | descriptive |  | 35 | -0.637 | 0.637 | 0.914 | worse_than_trivial | inf |  |
| T2_measured | fuellib_T_ge_283K | refutas_raw | designated | descriptive |  | 4 | 0.0319 | 0.0319 |  | n_too_small | 0.186 | 1 |
| T2_measured | fuellib_T_lt_263K | refutas_raw | designated | descriptive |  | 4 | -0.664 | 0.664 | 1 | n_too_small | inf |  |
| T2_measured | pooled_T_263_283K | refutas_raw | designated | descriptive |  | 38 | -0.057 | 0.063 | 0 | beats_trivial | inf |  |
| T2_measured | pooled_T_ge_283K | refutas_raw | designated | descriptive |  | 54 | 0.00976 | 0.0356 |  | beats_trivial | 0.186 | 0.981 |
| T2_measured | pooled_T_lt_263K | refutas_raw | designated | descriptive |  | 65 | -0.547 | 0.57 | 0.738 | worse_than_trivial | inf |  |
| T3_ajftd | ajftd_bulk_locked_all | const | designated | already_read |  | 7 | 0.201 | 0.654 |  | not_applicable | 1.14 | 0.946 |
| T3_ajftd | ajftd_bulk_locked_all | hgb | designated | already_read |  | 7 | -0.118 | 0.323 |  | n_too_small | 1.66 | 0.965 |
| T3_ajftd | ajftd_bulk_locked_all | ridge | designated | already_read |  | 7 | -0.0177 | 0.514 |  | n_too_small | 2.26 | 0.96 |
| T3_ajftd | ajftd_bulk_logo_all | const | designated | dev |  | 31 | 0.00171 | 0.608 |  | not_applicable | 1.14 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | hgb | designated | dev |  | 31 | -0.214 | 0.5 |  | no_skill_shown | 1.66 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | ridge | designated | dev |  | 31 | -0.204 | 0.854 |  | no_skill_shown | 2.26 | 0.901 |
| diesel | face_diesel | refutas | designated | descriptive |  | 4 | -0.412 | 0.412 |  | n_too_small |  |  |
| diesel | face_surrogates | refutas | designated | descriptive |  | 3 | 0.12 | 0.12 |  | not_applicable |  |  |
| CRC647_r3 | crc647_r3 | Bilinear pairwise | carried | dev |  | 45 | 0.024 | 0.0803 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Bilinear pairwise | carried | fold_b |  | 12 | -0.0228 | 0.0923 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Lab-trained structured residual (zero-shot) | carried | dev |  | 45 | -0.395 | 0.395 | 0.407 | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | Lab-trained structured residual (zero-shot) | carried | fold_b |  | 12 | -0.37 | 0.372 | 0.361 | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | Raw Refutas | carried | dev |  | 45 | -0.4 | 0.401 | 0.407 | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | Raw Refutas | carried | fold_b |  | 12 | -0.375 | 0.378 | 0.361 | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | Structured GP | carried | dev |  | 45 | -6.61e-05 | 0.0652 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Structured GP | carried | fold_b |  | 12 | 0.00252 | 0.0944 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Structured ridge | carried | dev |  | 45 | 0.000255 | 0.0622 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Structured ridge | carried | fold_b |  | 12 | -0.00152 | 0.0935 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Walther (D341) | carried | dev |  | 45 | -0.374 | 0.375 | 0.4 | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | Walther (D341) | carried | fold_b |  | 12 | -0.349 | 0.353 | 0.361 | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | Walther + recalibration | carried | dev |  | 45 | 0.069 | 0.109 | 0 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | Walther + recalibration | carried | fold_b |  | 12 | 0.0934 | 0.156 | 0 | beats_trivial |  |  |
| T2_measured | ajftd_gcxgc_T_263_283K | api_as_shipped | known | descriptive |  | 3 | -0.0157 | 0.0157 | 0 | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_263_283K | api_as_shipped | reconstructed_route | descriptive |  | 3 | -0.0157 | 0.0157 | 0 | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_263_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 3 | -0.0199 | 0.0199 | 0 | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | api_as_shipped | known | descriptive |  | 17 | 0.0243 | 0.068 |  | beats_trivial |  | 0.75 |
| T2_measured | ajftd_gcxgc_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 17 | 0.0243 | 0.068 |  | beats_trivial |  |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 17 | -0.0746 | 0.0763 |  | beats_trivial |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | api_as_shipped | known | descriptive |  | 28 | -0.263 | 0.359 | 0.321 | no_skill_shown |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 28 | -0.263 | 0.359 | 0.321 | no_skill_shown |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 28 | -0.0209 | 0.134 | 0.0357 | beats_trivial |  |  |
| T2_measured | av3322_T_263_283K | api_as_shipped | known | descriptive |  | 35 | -0.0128 | 0.0422 | 0 | beats_trivial |  |  |
| T2_measured | av3322_T_263_283K | api_as_shipped | reconstructed_route | descriptive |  | 35 | -0.0128 | 0.0422 | 0 | beats_trivial |  |  |
| T2_measured | av3322_T_263_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 35 | 0.0199 | 0.042 | 0 | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | api_as_shipped | known | descriptive |  | 35 | 0.0477 | 0.0529 |  | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 35 | 0.0477 | 0.0529 |  | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 35 | 0.0225 | 0.0326 |  | beats_trivial |  |  |
| T2_measured | av3322_T_lt_263K | api_as_shipped | known | descriptive |  | 35 | -0.604 | 0.604 | 0.829 | worse_than_trivial |  |  |
| T2_measured | av3322_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 35 | -0.604 | 0.604 | 0.829 | worse_than_trivial |  |  |
| T2_measured | av3322_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 35 | 0.0607 | 0.116 | 0 | beats_trivial |  |  |
| T2_measured | fuellib_T_ge_283K | api_as_shipped | known | descriptive |  | 4 | 0.0709 | 0.0709 |  | n_too_small |  | 1 |
| T2_measured | fuellib_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 4 | 0.0709 | 0.0709 |  | n_too_small |  |  |
| T2_measured | fuellib_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 3 | -0.0332 | 0.0332 |  | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | api_as_shipped | known | descriptive |  | 4 | -0.651 | 0.651 | 1 | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 4 | -0.651 | 0.651 | 1 | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 4 | -0.0272 | 0.125 | 0 | n_too_small |  |  |
| T2_measured | pooled_T_263_283K | api_as_shipped | known | descriptive |  | 38 | -0.013 | 0.0401 | 0 | beats_trivial |  |  |
| T2_measured | pooled_T_263_283K | api_as_shipped | reconstructed_route | descriptive |  | 38 | -0.013 | 0.0401 | 0 | beats_trivial |  |  |
| T2_measured | pooled_T_263_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 38 | 0.0167 | 0.0403 | 0 | beats_trivial |  |  |
| T2_measured | pooled_T_ge_283K | api_as_shipped | known | descriptive |  | 54 | 0.0413 | 0.0584 |  | beats_trivial |  | 0.8 |
| T2_measured | pooled_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 54 | 0.0413 | 0.0584 |  | beats_trivial |  |  |
| T2_measured | pooled_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 53 | -0.00893 | 0.0459 |  | beats_trivial |  |  |
| T2_measured | pooled_T_lt_263K | api_as_shipped | known | descriptive |  | 65 | -0.472 | 0.513 | 0.646 | worse_than_trivial |  |  |
| T2_measured | pooled_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 65 | -0.472 | 0.513 | 0.646 | worse_than_trivial |  |  |
| T2_measured | pooled_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 65 | 0.0254 | 0.125 | 0.0154 | beats_trivial |  |  |

## density_kg_m3

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_r3 | ideal_volume | designated | dev |  | 45 | -14.2 | 14.3 |  | no_skill_shown | 23.6 | 0.822 |
| CRC647_r3 | crc647_r3 | ideal_volume | designated | fold_b |  | 12 | -15.1 | 15.1 |  | no_skill_shown | 23.6 | 0.833 |
| T1_lab | T_263_283K | ideal_volume | designated | already_read |  | 1 | -0.525 | 0.525 |  | not_applicable |  |  |
| T1_lab | T_263_283K | ideal_volume | designated | dev |  | 10 | 0.256 | 0.717 |  | not_applicable |  |  |
| T1_lab | T_ge_283K | ideal_volume | designated | already_read |  | 41 | 0.185 | 1.18 |  | not_applicable |  |  |
| T1_lab | T_ge_283K | ideal_volume | designated | dev |  | 167 | 0.258 | 1.35 |  | not_applicable |  |  |
| T1_lab | T_lt_263K | ideal_volume | designated | dev |  | 8 | 0.55 | 0.837 |  | not_applicable |  |  |
| T1_lab | thermoml2020_new_systems | ideal_volume | designated | dev |  | 4 | 1.28 | 1.34 |  | not_applicable | 4.59 | 0.75 |
| T2_measured | ajftd_gcxgc_T_263_283K | ideal_volume | designated | descriptive |  | 3 | -8.13 | 8.13 |  | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | ideal_volume | designated | descriptive |  | 32 | -5.3 | 7.55 |  | beats_trivial |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | ideal_volume | designated | descriptive |  | 3 | -7.49 | 7.49 |  | n_too_small |  |  |
| T2_measured | av3322_T_263_283K | ideal_volume | designated | descriptive |  | 35 | -1.87 | 2.22 |  | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | ideal_volume | designated | descriptive |  | 35 | -2.22 | 2.47 |  | beats_trivial |  |  |
| T2_measured | av3322_T_lt_263K | ideal_volume | designated | descriptive |  | 35 | -0.874 | 1.85 |  | beats_trivial |  |  |
| T2_measured | fuellib_T_263_283K | ideal_volume | designated | descriptive |  | 4 | -7.21 | 7.21 |  | n_too_small |  |  |
| T2_measured | fuellib_T_ge_283K | ideal_volume | designated | descriptive |  | 4 | -8.18 | 8.18 |  | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | ideal_volume | designated | descriptive |  | 2 | -7.74 | 7.74 |  | n_too_small |  |  |
| T2_measured | pooled_T_263_283K | ideal_volume | designated | descriptive |  | 41 | -2.67 | 2.97 |  | beats_trivial |  |  |
| T2_measured | pooled_T_ge_283K | ideal_volume | designated | descriptive |  | 69 | -3.8 | 4.97 |  | beats_trivial |  |  |
| T2_measured | pooled_T_lt_263K | ideal_volume | designated | descriptive |  | 39 | -1.54 | 2.41 |  | beats_trivial |  |  |
| diesel | face_diesel | ideal_volume | designated | descriptive |  | 4 | 5.84 | 15 |  | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_263_283K | api_as_shipped | known | descriptive |  | 3 | -4.82 | 4.82 |  | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | api_as_shipped | known | descriptive |  | 32 | -1.55 | 4.96 |  | beats_trivial |  | 0.5 |
| T2_measured | ajftd_gcxgc_T_lt_263K | api_as_shipped | known | descriptive |  | 3 | -14.9 | 14.9 |  | n_too_small |  |  |
| T2_measured | av3322_T_263_283K | api_as_shipped | known | descriptive |  | 35 | 3.63 | 3.63 |  | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | api_as_shipped | known | descriptive |  | 35 | 4.67 | 4.67 |  | beats_trivial |  |  |
| T2_measured | av3322_T_lt_263K | api_as_shipped | known | descriptive |  | 35 | -6.56 | 6.56 |  | no_skill_shown |  |  |
| T2_measured | fuellib_T_263_283K | api_as_shipped | known | descriptive |  | 4 | -4.39 | 4.39 |  | n_too_small |  |  |
| T2_measured | fuellib_T_ge_283K | api_as_shipped | known | descriptive |  | 4 | -4.11 | 4.11 |  | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | api_as_shipped | known | descriptive |  | 2 | -15.8 | 15.8 |  | n_too_small |  |  |
| T2_measured | pooled_T_263_283K | api_as_shipped | known | descriptive |  | 41 | 2.44 | 3.75 |  | beats_trivial |  |  |
| T2_measured | pooled_T_ge_283K | api_as_shipped | known | descriptive |  | 69 | 1.49 | 4.83 |  | beats_trivial |  | 0.5 |
| T2_measured | pooled_T_lt_263K | api_as_shipped | known | descriptive |  | 39 | -7.44 | 7.44 |  | no_skill_shown |  |  |

## flash_point_C

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_r3 | ideal_le_chatelier | designated | dev |  | 45 | 0.436 | 4.03 | 0.133 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | ideal_le_chatelier | designated | fold_b |  | 12 | 2.46 | 3.72 | 0.25 | beats_trivial |  |  |
| T1_lab | lab_mixtures_31_systems | ideal_le_chatelier | designated | dev |  | 31 | -3.3 | 3.74 | 0.00794 | not_applicable | 3.7 | 0.935 |
| T2_measured | ajftd_gcxgc | api_as_shipped | designated | descriptive |  | 31 | -23.4 | 24 | 0.0323 | worse_than_trivial |  | 0 |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | designated | descriptive |  | 31 | -0.634 | 3.38 | 0.0968 | no_skill_shown | 3.7 | 0.903 |
| T2_measured | av3322 | api_as_shipped | designated | descriptive |  | 35 | -19.7 | 19.7 | 0 | worse_than_trivial |  | 0 |
| T2_measured | av3322 | ideal_le_chatelier | designated | descriptive |  | 35 | -1.06 | 1.5 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | designated | descriptive |  | 28 | -1.07 | 1.56 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | designated | descriptive |  | 32 | -0.387 | 3.1 | 0.0938 | beats_trivial | 3.7 | 0.906 |
| T2_measured | family_other | ideal_le_chatelier | designated | descriptive |  | 10 | -0.921 | 3.79 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | fuellib | api_as_shipped | designated | descriptive |  | 4 | -24.8 | 24.8 | 0 | n_too_small |  | 0 |
| T2_measured | fuellib | ideal_le_chatelier | designated | descriptive |  | 4 | 1.36 | 5.82 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | pooled | api_as_shipped | designated | descriptive |  | 70 | -21.7 | 21.9 | 0.0143 | worse_than_trivial |  | 0 |
| T2_measured | pooled | ideal_le_chatelier | designated | descriptive |  | 70 | -0.734 | 2.58 | 0.0571 | beats_trivial | 3.7 | 0.943 |
| T3_ajftd | ajftd_bulk_locked_all | const | designated | already_read |  | 8 | 1.28 | 4.69 |  | not_applicable | 7.03 | 0.864 |
| T3_ajftd | ajftd_bulk_locked_all | hgb | designated | already_read |  | 8 | 0.322 | 2.55 |  | no_skill_shown | 13.4 | 0.992 |
| T3_ajftd | ajftd_bulk_locked_all | ridge | designated | already_read |  | 8 | 1.33 | 3.76 |  | no_skill_shown | 23.8 | 1 |
| T3_ajftd | ajftd_bulk_logo_all | const | designated | dev |  | 32 | 0.0348 | 3.94 |  | not_applicable | 7.03 | 0.905 |
| T3_ajftd | ajftd_bulk_logo_all | hgb | designated | dev |  | 32 | 0.0501 | 4.64 |  | no_skill_shown | 13.4 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | hgb | designated | dev |  | 32 | 0.0501 | 4.64 | 0.101 | no_skill_shown |  |  |
| T3_ajftd | ajftd_bulk_logo_all | ridge | designated | dev |  | 32 | -1.2 | 6.72 |  | no_skill_shown | 23.8 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | ridge | designated | dev |  | 32 | -1.2 | 6.72 | 0.149 | no_skill_shown |  |  |
| diesel | face_diesel | ideal_le_chatelier | designated | descriptive |  | 4 | -12.8 | 12.8 | 0 | n_too_small |  |  |
| CRC647_r3 | crc647_r3 | ideal_le_chatelier | lfl_overlay | dev |  | 45 | -0.726 | 4.09 | 0.0444 | not_applicable |  |  |
| CRC647_r3 | crc647_r3 | ideal_le_chatelier | lfl_overlay | fold_b |  | 12 | 1.41 | 3.26 | 0.0833 | not_applicable |  |  |
| CRC647_r3 | crc647_r3 | unifac_le_chatelier | unifac | dev |  | 45 | -0.485 | 3.99 | 0.0667 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | unifac_le_chatelier | unifac | fold_b |  | 12 | 1.6 | 3.35 | 0.0833 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | unifac_le_chatelier | unifac_lfl_overlay | dev |  | 45 | -1.68 | 4.13 | 0.0444 | not_applicable |  |  |
| CRC647_r3 | crc647_r3 | unifac_le_chatelier | unifac_lfl_overlay | fold_b |  | 12 | 0.514 | 3.15 | 0.0833 | not_applicable |  |  |
| T1_lab | lab_mixtures_31_systems | ideal_le_chatelier | lfl_overlay | dev |  | 31 | -3.25 | 3.69 | 0.00794 | not_applicable |  |  |
| T1_lab | lab_mixtures_31_systems | ideal_le_chatelier | mesitylene_fixed | dev |  | 31 | -3.59 | 3.89 | 0.00794 | not_applicable |  |  |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | lfl_overlay | descriptive |  | 31 | -0.478 | 3.31 | 0.0968 | no_skill_shown | 3.7 | 0.903 |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 31 | -0.841 | 3.65 | 0.0968 | no_skill_shown | 3.7 | 0.903 |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 31 | -0.994 | 3.72 | 0.0968 | no_skill_shown | 3.7 | 0.903 |
| T2_measured | av3322 | ideal_le_chatelier | lfl_overlay | descriptive |  | 35 | -1.14 | 1.55 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | av3322 | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 35 | -1.46 | 1.69 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | av3322 | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 35 | -1.39 | 1.62 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | lfl_overlay | descriptive |  | 28 | -1.13 | 1.61 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 28 | -1.42 | 1.7 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 28 | -1.35 | 1.65 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | lfl_overlay | descriptive |  | 32 | -0.271 | 3.05 | 0.0938 | beats_trivial | 3.7 | 0.906 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 32 | -0.661 | 3.38 | 0.0938 | beats_trivial | 3.7 | 0.906 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 32 | -0.775 | 3.41 | 0.0938 | beats_trivial | 3.7 | 0.906 |
| T2_measured | family_other | ideal_le_chatelier | lfl_overlay | descriptive |  | 10 | -0.799 | 3.77 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | family_other | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 10 | -1.12 | 4.09 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | family_other | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 10 | -1.24 | 4.11 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | fuellib | ideal_le_chatelier | lfl_overlay | descriptive |  | 4 | 1.57 | 5.86 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | fuellib | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 4 | 1.31 | 6.12 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | fuellib | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 4 | 1.1 | 6.07 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | pooled | ideal_le_chatelier | lfl_overlay | descriptive |  | 70 | -0.692 | 2.58 | 0.0571 | beats_trivial | 3.7 | 0.943 |
| T2_measured | pooled | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 70 | -1.03 | 2.81 | 0.0571 | beats_trivial | 3.7 | 0.943 |
| T2_measured | pooled | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 70 | -1.07 | 2.81 | 0.0571 | beats_trivial | 3.7 | 0.943 |
| T3_ajftd | ajftd_bulk_locked_all | api_2b7_t10 | unverified_coefficients | already_read |  | 8 | 6.27 | 6.35 |  | no_skill_shown | 11.1 | 0.936 |
| T3_ajftd | ajftd_bulk_locked_all | butler_t10 | unverified_coefficients | already_read |  | 8 | 2.65 | 3.91 |  | no_skill_shown | 7.33 | 0.941 |
| T3_ajftd | ajftd_bulk_logo_all | api_2b7_t10 | unverified_coefficients | dev |  | 32 | 1.96 | 9.5 |  | worse_than_trivial | 11.1 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | butler_t10 | unverified_coefficients | dev |  | 32 | -1.21 | 6.41 |  | no_skill_shown | 7.33 | 0.901 |

## freezing_point_C

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_r3 | ideal_svl | designated | dev |  | 45 | 4.14 | 5.26 |  | no_skill_shown |  |  |
| CRC647_r3 | crc647_r3 | ideal_svl | designated | fold_b |  | 12 | 2.59 | 3.14 |  | no_skill_shown |  |  |
| T1_lab | sle_alkane_blends | ideal_svl | designated | dev |  | 8 | -1.1 | 2.02 |  | not_applicable | inf |  |
| T2_measured | AJFTD | ideal_svl | designated | descriptive |  | 12 | 9.26 | 9.26 |  | no_skill_shown |  |  |
| T2_measured | AJFTD_locked | ideal_svl | designated | already_read |  | 8 | 4.64 | 9.56 |  | worse_than_trivial |  |  |
| T2_measured | AV-33-22 | ideal_svl | designated | descriptive |  | 34 | 7.54 | 7.54 |  | no_skill_shown |  |  |
| T2_measured | FuelLib | ideal_svl | designated | descriptive |  | 3 | 7.92 | 7.92 |  | n_too_small |  |  |
| T2_measured | pooled | ideal_svl | designated | descriptive |  | 49 | 7.99 | 7.99 |  | no_skill_shown |  |  |
| T2_measured | pooled_locked | ideal_svl | designated | already_read |  | 8 | 4.64 | 9.56 |  | worse_than_trivial |  |  |
| T3_ajftd | ajftd_bulk_locked_all | const | designated | already_read |  | 7 | 1.25 | 3.62 |  | not_applicable | 8.47 | 0.896 |
| T3_ajftd | ajftd_bulk_locked_all | hgb | designated | already_read |  | 7 | -0.818 | 3.23 |  | n_too_small | 6.83 | 0.931 |
| T3_ajftd | ajftd_bulk_locked_all | ridge | designated | already_read |  | 7 | -0.46 | 3.71 |  | n_too_small | 10.5 | 0.954 |
| T3_ajftd | ajftd_bulk_logo_all | const | designated | dev |  | 33 | -0.144 | 4.96 |  | not_applicable | 8.47 | 0.902 |
| T3_ajftd | ajftd_bulk_logo_all | hgb | designated | dev |  | 33 | -0.318 | 3.52 |  | beats_trivial | 6.83 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | ridge | designated | dev |  | 33 | -0.524 | 5.59 |  | no_skill_shown | 10.5 | 0.901 |
| CRC647_r3 | crc647_r3 | coutinho_library | carried | dev |  | 45 | -10.3 | 10.6 |  | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | coutinho_library | carried | fold_b |  | 12 | -11.7 | 11.7 |  | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | wat_ideal_library | carried | dev |  | 45 | 6.2 | 6.93 |  | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | wat_ideal_library | carried | fold_b |  | 12 | 4.61 | 4.76 |  | worse_than_trivial |  |  |
| CRC647_r3 | crc647_r3 | wat_minus_devfit_offset | carried | dev |  | 45 | -2.41e-15 | 3.79 |  | no_skill_shown |  |  |
| CRC647_r3 | crc647_r3 | wat_minus_devfit_offset | carried | fold_b |  | 12 | -1.59 | 2.26 |  | no_skill_shown |  |  |
| T1_lab | sle_alkane_blends | ideal_svl | sle2_single_compilation | dev |  | 38 | -4.02 | 7.03 |  | not_applicable | 24.2 |  |
| T1_lab | sle_alkane_blends | ideal_svl | sle2_two_compilations | dev |  | 11 | -1.41 | 2.13 |  | not_applicable | 12.3 |  |
| T2_measured | AJFTD | coutinho_f1 | f1_recalibrated | descriptive |  | 12 | -1.48 | 5.62 |  | no_skill_shown |  |  |
| T2_measured | AJFTD | ridge_p | p_ridge | descriptive |  | 13 | 4.99 | 8.41 |  | no_skill_shown |  |  |
| T2_measured | AJFTD_locked | coutinho_f1 | f1_recalibrated | already_read |  | 8 | -4.06 | 8.85 |  | no_skill_shown |  |  |
| T2_measured | AJFTD_locked | ridge_p | p_ridge | already_read |  | 8 | -3.4 | 4.91 |  | no_skill_shown |  |  |
| T2_measured | AV-33-22 | coutinho_f1 | f1_recalibrated | descriptive |  | 34 | 0.192 | 2.92 |  | beats_trivial |  |  |
| T2_measured | AV-33-22 | ridge_p | p_ridge | descriptive |  | 34 | -0.0731 | 3.58 |  | beats_trivial |  |  |
| T2_measured | FuelLib | coutinho_f1 | f1_recalibrated | descriptive |  | 3 | 0.537 | 0.977 |  | n_too_small |  |  |
| T2_measured | FuelLib | ridge_p | p_ridge | descriptive |  | 3 | -0.781 | 2.92 |  | n_too_small |  |  |
| T2_measured | pooled | coutinho_f1 | f1_recalibrated | descriptive |  | 49 | -0.197 | 3.46 |  | beats_trivial |  |  |
| T2_measured | pooled | ridge_p | p_ridge | descriptive |  | 50 | 1.2 | 4.8 |  | beats_trivial |  |  |
| T2_measured | pooled_locked | coutinho_f1 | f1_recalibrated | already_read |  | 8 | -4.06 | 8.85 |  | no_skill_shown |  |  |
| T2_measured | pooled_locked | ridge_p | p_ridge | already_read |  | 8 | -3.4 | 4.91 |  | no_skill_shown |  |  |

## heat_of_combustion_net_MJ_kg

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_r3 | ga_comp_H | designated | dev |  | 45 | 0.376 | 0.376 |  | worse_than_trivial | 0.224 | 0.0222 |
| CRC647_r3 | crc647_r3 | ga_comp_H | designated | fold_b |  | 12 | 0.395 | 0.395 |  | worse_than_trivial | 0.224 | 0.0833 |
| T2_measured | ajftd_gcxgc_measured | ga_comp_H | designated | descriptive |  | 28 | 0.348 | 0.35 |  | no_skill_shown |  |  |
| T2_measured | av3322 | ga_comp_H | designated | descriptive |  | 35 | 0.0938 | 0.0971 |  | no_skill_shown |  |  |
| T2_measured | measured_pooled | ga_comp_H | designated | descriptive |  | 28 | 0.348 | 0.35 |  | no_skill_shown |  |  |
| T3_ajftd | R-S4b-1_measured | C_all | designated | dev |  | 16 |  | 0.407 |  | worse_than_trivial |  |  |
| T3_ajftd | R-S4b-1_measured | hgb | designated | dev |  | 16 |  | 0.203 |  | beats_trivial |  |  |
| T3_ajftd | R-S4b-1_measured | ridge | designated | dev |  | 16 |  | 0.26 |  | beats_trivial |  |  |
| T3_ajftd | R-S4b-2_calculated | C_all | designated | dev |  | 19 |  | 0.158 |  | no_skill_shown |  |  |
| T3_ajftd | R-S4b-2_calculated | hgb | designated | dev |  | 19 |  | 0.0427 |  | beats_trivial |  |  |
| T3_ajftd | R-S4b-2_calculated | ridge | designated | dev |  | 19 |  | 0.0984 |  | no_skill_shown |  |  |
| T3_ajftd | R-S4b-3_mixed | C_all | designated | dev |  | 33 |  | 0.214 |  | worse_than_trivial |  |  |
| T3_ajftd | R-S4b-3_mixed | hgb | designated | dev |  | 33 |  | 0.155 |  | beats_trivial |  |  |
| T3_ajftd | R-S4b-3_mixed | ridge | designated | dev |  | 33 |  | 0.16 |  | beats_trivial |  |  |
| T3_ajftd | R-S4b-4_locked_measured | C_all | designated | already_read |  | 5 |  | 0.582 |  | n_too_small |  |  |
| T3_ajftd | R-S4b-4_locked_measured | hgb | designated | already_read |  | 5 |  | 0.252 |  | n_too_small |  |  |
| T3_ajftd | R-S4b-4_locked_measured | ridge | designated | already_read |  | 5 |  | 0.276 |  | n_too_small |  |  |
| T3_ajftd | ajftd_bulk_locked_all | const | designated | already_read |  | 8 | -0.0568 | 0.183 |  | not_applicable | 0.493 | 0.884 |
| T3_ajftd | ajftd_bulk_locked_all | hgb | designated | already_read |  | 8 | -0.00257 | 0.0942 |  | beats_trivial | 0.36 | 0.953 |
| T3_ajftd | ajftd_bulk_locked_all | ridge | designated | already_read |  | 8 | -0.0288 | 0.111 |  | beats_trivial | 0.377 | 0.944 |
| T3_ajftd | ajftd_bulk_logo_all | const | designated | dev |  | 33 | 0.00121 | 0.214 |  | not_applicable | 0.493 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | hgb | designated | dev |  | 33 | -0.00149 | 0.155 |  | beats_trivial | 0.36 | 0.901 |
| T3_ajftd | ajftd_bulk_logo_all | ridge | designated | dev |  | 33 | -0.00713 | 0.16 |  | beats_trivial | 0.377 | 0.901 |
| CRC647_r3 | crc647_r3 | ga_const_H | const_H | dev |  | 45 | 0.433 | 0.433 |  | worse_than_trivial | 0.194 | 0.0222 |
| CRC647_r3 | crc647_r3 | ga_const_H | const_H | fold_b |  | 12 | 0.47 | 0.47 |  | worse_than_trivial | 0.194 | 0.0833 |
| T2_measured | ajftd_gcxgc_measured | ga_const_H | const_H | descriptive |  | 28 | 0.553 | 0.553 |  | worse_than_trivial |  |  |
| T2_measured | av3322 | ga_const_H | const_H | descriptive |  | 35 | 0.102 | 0.113 |  | no_skill_shown |  |  |
| T2_measured | measured_pooled | ga_const_H | const_H | descriptive |  | 28 | 0.553 | 0.553 |  | worse_than_trivial |  |  |

## smoke_point_mm

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T3_ajftd | ajftd100_bulk | arom+diarom_fixed | designated | dev |  |  | 0.00581 | 0.904 |  | beats_trivial |  |  |
| T3_ajftd | ajftd100_bulk_groupTcode | arom+diarom_fixed | designated | dev |  |  | -0.477 | 1.31 |  | beats_trivial |  |  |
| T3_ajftd | ajftd100_bulk_groupTestDate | arom+diarom_fixed | designated | dev |  |  | -0.00753 | 0.895 |  | beats_trivial |  |  |
| T3_ajftd | ajftd_broad_bulk | arom+diarom_fixed | designated | dev |  |  | -0.000923 | 1.68 |  | no_skill_shown |  |  |
| T3_ajftd | ajftd100_bulk | const_mean | carried | dev |  |  | 2.1e-15 | 1.14 |  | not_applicable |  |  |
| T3_ajftd | ajftd100_bulk | selected_in_fold | carried | dev |  |  | 0.00581 | 0.904 |  | beats_trivial |  |  |
| T3_ajftd | ajftd100_bulk_groupTcode | const_mean | carried | dev |  |  | -0.161 | 1.47 |  | not_applicable |  |  |
| T3_ajftd | ajftd100_bulk_groupTcode | selected_in_fold | carried | dev |  |  | -0.453 | 1.29 |  | beats_trivial |  |  |
| T3_ajftd | ajftd100_bulk_groupTestDate | const_mean | carried | dev |  |  | 6.75e-16 | 1.13 |  | not_applicable |  |  |
| T3_ajftd | ajftd100_bulk_groupTestDate | selected_in_fold | carried | dev |  |  | -0.00753 | 0.895 |  | beats_trivial |  |  |
| T3_ajftd | ajftd_broad_bulk | const_mean | carried | dev |  |  | -5.77e-05 | 1.65 |  | not_applicable |  |  |
| T3_ajftd | ajftd_broad_bulk | selected_in_fold | carried | dev |  |  | -1.61e-05 | 1.68 |  | no_skill_shown |  |  |

## pour_point_C

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_measured_bulk_n57 | const | designated | dev |  | 57 | -4.24e-15 | 1.82 |  | not_applicable |  |  |
| CRC647_r3 | crc647_measured_bulk_n57 | ridge | designated | dev |  | 57 | 0.103 | 1.92 |  | no_skill_shown |  |  |
| CRC647_r3 | crc647_measured_bulk_n57 | freezing_offset | uses_measured_freezing | dev |  | 57 | -1.37e-15 | 3.9 |  | not_applicable |  |  |

## ron

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1_lab | octane_dev_cv5_stored_folds | base_unchanged | designated | dev |  | 18 | -5.07 | 5.1 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | base_unchanged | designated | dev |  | 18 | -5.07 | 5.1 |  | worse_than_trivial |  |  |
| T1_lab | octane_locked_thin_test | base_unchanged | designated | already_read |  | 4 | -5.69 | 5.69 |  | n_too_small |  |  |
| T1_lab | octane_dev_cv5_stored_folds | constant_mean | carried | dev |  | 18 | -0.115 | 3.22 |  | not_applicable |  |  |
| T1_lab | octane_dev_cv5_stored_folds | gbm_resid | carried | dev |  | 18 | -0.264 | 2.55 |  | beats_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | ghosh_lumps | carried | dev |  | 18 | -0.514 | 3.19 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_cv5_stored_folds | gp_resid | carried | dev |  | 18 | -0.537 | 2.29 |  | beats_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | linear_mole_prior | carried | dev |  | 16 | -6.76 | 6.79 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | linear_volume_pooled_bov | carried | dev |  | 18 | -0.0146 | 2.73 |  | beats_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | linear_volume_prior | carried | dev |  | 18 | -6.46 | 6.49 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | constant_mean | carried | dev |  | 18 | 0.00331 | 3.25 |  | not_applicable |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | gbm_resid | carried | dev |  | 18 | 0.0345 | 3.15 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | ghosh_lumps | carried | dev |  | 18 | -0.396 | 2.98 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | gp_resid | carried | dev |  | 18 | -0.396 | 2.87 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | linear_mole_prior | carried | dev |  | 16 | -6.76 | 6.79 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | linear_volume_pooled_bov | carried | dev |  | 18 | 0.16 | 2.69 |  | beats_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | linear_volume_prior | carried | dev |  | 18 | -6.46 | 6.49 |  | worse_than_trivial |  |  |
| T1_lab | octane_locked_thin_test | constant_mean | carried | already_read |  | 4 | -0.766 | 2.75 |  | not_applicable |  |  |
| T1_lab | octane_locked_thin_test | gbm_resid | carried | already_read |  | 4 | -1.02 | 1.61 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | ghosh_lumps | carried | already_read |  | 4 | -0.588 | 1.48 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | gp_resid | carried | already_read |  | 4 | -1.48 | 1.8 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | linear_mole_prior | carried | already_read |  | 4 | -5.5 | 5.5 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | linear_volume_pooled_bov | carried | already_read |  | 4 | -0.492 | 1.84 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | linear_volume_prior | carried | already_read |  | 4 | -6.02 | 6.02 |  | n_too_small |  |  |

## mon

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1_lab | octane_dev_cv5_stored_folds | base_unchanged | designated | dev |  | 18 | -3.26 | 3.26 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | base_unchanged | designated | dev |  | 18 | -3.26 | 3.26 |  | worse_than_trivial |  |  |
| T1_lab | octane_locked_thin_test | base_unchanged | designated | already_read |  | 4 | -4.27 | 4.27 |  | n_too_small |  |  |
| T1_lab | octane_dev_cv5_stored_folds | constant_mean | carried | dev |  | 18 | 0.012 | 1.62 |  | not_applicable |  |  |
| T1_lab | octane_dev_cv5_stored_folds | gbm_resid | carried | dev |  | 18 | -0.119 | 1.11 |  | beats_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | ghosh_lumps | carried | dev |  | 18 | -0.104 | 1.74 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_cv5_stored_folds | gp_resid | carried | dev |  | 18 | -0.0971 | 1.15 |  | beats_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | linear_mole_prior | carried | dev |  | 16 | -4.64 | 4.77 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | linear_volume_pooled_bov | carried | dev |  | 18 | 0.19 | 1.25 |  | beats_trivial |  |  |
| T1_lab | octane_dev_cv5_stored_folds | linear_volume_prior | carried | dev |  | 18 | -4.36 | 4.5 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | constant_mean | carried | dev |  | 18 | 3.16e-15 | 1.65 |  | not_applicable |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | gbm_resid | carried | dev |  | 18 | -0.222 | 1.3 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | ghosh_lumps | carried | dev |  | 18 | -0.0942 | 1.8 |  | no_skill_shown |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | gp_resid | carried | dev |  | 18 | -0.148 | 1.25 |  | beats_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | linear_mole_prior | carried | dev |  | 16 | -4.64 | 4.77 |  | worse_than_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | linear_volume_pooled_bov | carried | dev |  | 18 | 0.211 | 1.28 |  | beats_trivial |  |  |
| T1_lab | octane_dev_leave_one_blendstock_out | linear_volume_prior | carried | dev |  | 18 | -4.36 | 4.5 |  | worse_than_trivial |  |  |
| T1_lab | octane_locked_thin_test | constant_mean | carried | already_read |  | 4 | -1.23 | 1.77 |  | not_applicable |  |  |
| T1_lab | octane_locked_thin_test | gbm_resid | carried | already_read |  | 4 | -1.51 | 1.52 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | ghosh_lumps | carried | already_read |  | 4 | -0.936 | 0.983 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | gp_resid | carried | already_read |  | 4 | -1.02 | 1.22 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | linear_mole_prior | carried | already_read |  | 4 | -4.98 | 4.98 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | linear_volume_pooled_bov | carried | already_read |  | 4 | -0.884 | 1.44 |  | n_too_small |  |  |
| T1_lab | octane_locked_thin_test | linear_volume_prior | carried | already_read |  | 4 | -5.03 | 5.03 |  | n_too_small |  |  |

## designs_certified

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| designs | d4 | v4_designer | designated | descriptive |  |  | 0 |  |  | not_applicable |  |  |
