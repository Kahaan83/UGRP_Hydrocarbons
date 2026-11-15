# v5 BASE scoreboard

781 cells; 44 `no_evaluation`, 191 `n_too_small` (see `data_gaps.csv`).  Every T2 fuel has been read, so T2 cells are `descriptive`, not held-out.  Error is predicted minus measured; viscosity error is ln(pred/meas).  Bias, MAE and unsafe-side rate are in the property's unit (viscosity: ln).

## Readings

**R-BASE-1** (48 required (property, tier) cells): all exist with a status.

**R-BASE-2** (recomputed cells with a committed counterpart):

| cell | observed | committed | diff | tol | within | committed_source | cause |
|---|---|---|---|---|---|---|---|
| flash T2 AJFTD GCxGC, all 37 eligible rows: bias | 0.141 | -0.426 | 0.567 | 0.05 | False | results/v5/G1/summary.csv |  |
| flash T2 AJFTD GCxGC, all 37 eligible rows: MAE | 3.09 | 3.54 | -0.448 | 0.05 | False | results/v5/G1/summary.csv |  |
| flash T2 AJFTD GCxGC, headline one-per-group (31): bias | -0.152 | -0.426 | 0.275 | 0.05 | False | results/v5/G1/summary.csv | BASE scores one fuel per double-count group (31); G1 scored all 37 eligible rows (6 duplicate-group fuels); the 37-row rerun above matches G1 |
| ln viscosity AV-33-22 233.15 K: bias | -0.574 | -0.637 | 0.063 | 0.005 | False | results/v5/C1 (R-C1-5, rounded to 3 decimals) | C1 used S3's mapping route; BASE uses the shared composition route |
| ln viscosity AV-33-22 313.15 K: bias | 0.109 | 0.066 | 0.0426 | 0.005 | False | results/v5/C1 (R-C1-5, rounded to 3 decimals) | C1 used S3's mapping route; BASE uses the shared composition route |
| flash CRC 647 r3 dev (45): bias | 0.436 | 0.436 | 1.27e-09 | 0.05 | True | results/v5/L1/comparison.csv (before) |  |
| freezing T2 pooled F1: MAE | 3.46 | 3.46 | 0.00429 | 0.05 | True | results/v5/S5 (README, rounded) |  |

**R-BASE-3** (designated model vs the leave-one-group-out constant, aggregate rows without temperature split, by tier; a tier appears in several columns when its sets or models disagree, see the per-property tables):

| property | beats | no_skill | worse | too_small_only |
|---|---|---|---|---|
| density_kg_m3 | T2_measured | CRC647_r3 | - | diesel |
| flash_point_C | CRC647_r3, T2_measured | T2_measured, T3_ajftd | - | diesel |
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
| LFL overlay | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | lfl_overlay | descriptive | 31 | -0.152 | 0.00194 | 3 | 2.92 | 0.0968 | 0.0968 |
| LFL overlay | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 31 | -0.152 | -0.371 | 3 | 3.27 | 0.0968 | 0.0968 |
| LFL overlay | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | lfl_overlay | descriptive | 35 | -0.0275 | -0.112 | 1.2 | 1.24 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 35 | -0.0275 | -0.46 | 1.2 | 1.22 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | lfl_overlay | descriptive | 28 | 0.0904 | 0.013 | 1.29 | 1.31 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 28 | 0.0904 | -0.297 | 1.29 | 1.25 | 0 | 0 |
| LFL overlay | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | lfl_overlay | descriptive | 32 | -0.0446 | 0.0688 | 2.88 | 2.84 | 0.0938 | 0.0938 |
| LFL overlay | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 32 | -0.0446 | -0.329 | 2.88 | 3.13 | 0.0938 | 0.0938 |
| LFL overlay | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | lfl_overlay | descriptive | 10 | 0.3 | 0.417 | 2.57 | 2.55 | 0.1 | 0.1 |
| LFL overlay | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 10 | 0.3 | 0.0672 | 2.57 | 2.9 | 0.1 | 0.1 |
| LFL overlay | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | lfl_overlay | descriptive | 4 | 2.44 | 2.64 | 4.73 | 4.78 | 0.25 | 0.25 |
| LFL overlay | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 4 | 2.44 | 2.36 | 4.73 | 5.07 | 0.25 | 0.25 |
| LFL overlay | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | lfl_overlay | descriptive | 70 | 0.0585 | 0.0961 | 2.2 | 2.19 | 0.0571 | 0.0571 |
| LFL overlay | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 70 | 0.0585 | -0.26 | 2.2 | 2.35 | 0.0571 | 0.0571 |
| Corrected mesitylene Antoine set | flash_point_C | T1_lab | lab_mixtures_31_systems |  | ideal_le_chatelier | mesitylene_fixed | dev | 31 | -3.3 | -3.59 | 3.74 | 3.89 | 0.00794 | 0.00794 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 31 | -0.152 | -0.371 | 3 | 3.27 | 0.0968 | 0.0968 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | ajftd_gcxgc |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 31 | -0.152 | -0.523 | 3 | 3.34 | 0.0968 | 0.0968 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 35 | -0.0275 | -0.46 | 1.2 | 1.22 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | av3322 |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 35 | -0.0275 | -0.378 | 1.2 | 1.18 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 28 | 0.0904 | -0.297 | 1.29 | 1.25 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D56 |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 28 | 0.0904 | -0.221 | 1.29 | 1.21 | 0 | 0 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 32 | -0.0446 | -0.329 | 2.88 | 3.13 | 0.0938 | 0.0938 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_D93_or_IP170 |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 32 | -0.0446 | -0.441 | 2.88 | 3.18 | 0.0938 | 0.0938 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 10 | 0.3 | 0.0672 | 2.57 | 2.9 | 0.1 | 0.1 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | family_other |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 10 | 0.3 | -0.0506 | 2.57 | 2.92 | 0.1 | 0.1 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 4 | 2.44 | 2.36 | 4.73 | 5.07 | 0.25 | 0.25 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | fuellib |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 4 | 2.44 | 2.16 | 4.73 | 5.01 | 0.25 | 0.25 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive | 70 | 0.0585 | -0.26 | 2.2 | 2.35 | 0.0571 | 0.0571 |
| Corrected mesitylene Antoine set | flash_point_C | T2_measured | pooled |  | ideal_le_chatelier | mesitylene_fixed | descriptive | 70 | 0.0585 | -0.297 | 2.2 | 2.36 | 0.0571 | 0.0571 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 10 | -0.687 | -0.113 | 0.687 | 0.128 | 1 | 0.1 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 28 | -0.273 | -0.0661 | 0.326 | 0.146 | 0.0714 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 273 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.0421 | -0.0393 | 0.0421 | 0.0393 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 288 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | 0.0276 | -0.0404 | 0.0276 | 0.0404 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 293 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | 0.0412 | -0.042 | 0.0412 | 0.042 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 16 | -0.000297 | -0.109 | 0.0594 | 0.109 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc_T_263_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | -0.0421 | -0.0393 | 0.0421 | 0.0393 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 17 | -0.00696 | -0.108 | 0.0492 | 0.108 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | ajftd_gcxgc_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 28 | -0.386 | -0.0617 | 0.439 | 0.147 | 0.357 | 0.0357 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.574 | 0.113 | 0.574 | 0.135 | 0.886 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.207 | 0.0661 | 0.207 | 0.0769 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 273 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.00283 | 0.0393 | 0.0318 | 0.0461 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 288 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.068 | 0.0404 | 0.068 | 0.0418 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 293 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.0832 | 0.042 | 0.0832 | 0.0423 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322 | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.109 | 0.109 | 0.109 | 0.109 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322_T_263_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.00283 | 0.0393 | 0.0318 | 0.0461 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | 0.068 | 0.0404 | 0.068 | 0.0418 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | av3322_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 35 | -0.574 | 0.113 | 0.574 | 0.135 | 0.886 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | -0.628 | -0.0544 | 0.628 | 0.0965 | 1 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | -0.279 | -0.0727 | 0.279 | 0.0727 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | 0.0573 | -0.0513 | 0.0573 | 0.0513 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib | 373 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | 0.00342 | -0.0186 | 0.0144 | 0.0217 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | 0.0567 | -0.0431 | 0.0567 | 0.0431 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | fuellib_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 4 | -0.628 | -0.0544 | 0.628 | 0.0965 | 1 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 233 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 48 | -0.599 | 0.0567 | 0.599 | 0.131 | 0.917 | 0.0208 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 253 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 65 | -0.236 | 0.00602 | 0.259 | 0.105 | 0.0308 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 273 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 38 | -0.00593 | 0.0331 | 0.0327 | 0.0455 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 288 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 38 | 0.0648 | 0.0341 | 0.0648 | 0.0417 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 293 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 39 | 0.0789 | 0.0333 | 0.0789 | 0.0422 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 313 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 54 | 0.0733 | 0.0353 | 0.091 | 0.106 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled | 373 | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 3 | 0.0205 | -0.0186 | 0.0216 | 0.0217 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled_T_263_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 38 | -0.00593 | 0.0331 | 0.0327 | 0.0455 | 0 | 0 |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled_T_ge_283K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 53 | 0.045 | -0.00834 | 0.0626 | 0.0626 |  |  |
| K1 cold offset | viscosity_kin_cSt | T2_measured | pooled_T_lt_263K |  | refutas_raw+k1 | k1_cold_offset_loso | descriptive | 65 | -0.504 | 0.0341 | 0.527 | 0.138 | 0.692 | 0.0154 |
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
| T2_measured | ajftd_gcxgc_T_263_283K | refutas_raw | designated | descriptive |  | 3 | -0.0421 | 0.0421 | 0 | n_too_small | inf |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | refutas_raw | designated | descriptive |  | 17 | -0.00696 | 0.0492 |  | beats_trivial | 0.186 | 0.941 |
| T2_measured | ajftd_gcxgc_T_lt_263K | refutas_raw | designated | descriptive |  | 28 | -0.386 | 0.439 | 0.357 | worse_than_trivial | inf |  |
| T2_measured | av3322_T_263_283K | refutas_raw | designated | descriptive |  | 35 | -0.00283 | 0.0318 | 0 | beats_trivial | inf |  |
| T2_measured | av3322_T_ge_283K | refutas_raw | designated | descriptive |  | 35 | 0.068 | 0.068 |  | beats_trivial | 0.186 | 1 |
| T2_measured | av3322_T_lt_263K | refutas_raw | designated | descriptive |  | 35 | -0.574 | 0.574 | 0.886 | worse_than_trivial | inf |  |
| T2_measured | fuellib_T_ge_283K | refutas_raw | designated | descriptive |  | 4 | 0.0567 | 0.0567 |  | n_too_small | 0.186 | 1 |
| T2_measured | fuellib_T_lt_263K | refutas_raw | designated | descriptive |  | 4 | -0.628 | 0.628 | 1 | n_too_small | inf |  |
| T2_measured | pooled_T_263_283K | refutas_raw | designated | descriptive |  | 38 | -0.00593 | 0.0327 | 0 | beats_trivial | inf |  |
| T2_measured | pooled_T_ge_283K | refutas_raw | designated | descriptive |  | 54 | 0.045 | 0.0626 |  | beats_trivial | 0.186 | 0.981 |
| T2_measured | pooled_T_lt_263K | refutas_raw | designated | descriptive |  | 65 | -0.504 | 0.527 | 0.692 | worse_than_trivial | inf |  |
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
| T2_measured | ajftd_gcxgc_T_263_283K | api_as_shipped | known | descriptive |  | 3 | -0.0421 | 0.0421 | 0 | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_263_283K | api_as_shipped | reconstructed_route | descriptive |  | 3 | 0.0918 | 0.0918 | 0 | n_too_small |  | 1 |
| T2_measured | ajftd_gcxgc_T_263_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 3 | -0.0393 | 0.0393 | 0 | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | api_as_shipped | known | descriptive |  | 17 | -0.00696 | 0.0492 |  | beats_trivial |  | 0.941 |
| T2_measured | ajftd_gcxgc_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 17 | -0.0144 | 0.0472 |  | beats_trivial |  | 1 |
| T2_measured | ajftd_gcxgc_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 17 | -0.108 | 0.108 |  | no_skill_shown |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | api_as_shipped | known | descriptive |  | 28 | -0.386 | 0.439 | 0.357 | worse_than_trivial |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 28 | 0.0467 | 0.125 | 0.0357 | beats_trivial |  | 0.929 |
| T2_measured | ajftd_gcxgc_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 28 | -0.0617 | 0.147 | 0.0357 | beats_trivial |  |  |
| T2_measured | av3322_T_263_283K | api_as_shipped | known | descriptive |  | 35 | -0.00283 | 0.0318 | 0 | beats_trivial |  |  |
| T2_measured | av3322_T_263_283K | api_as_shipped | reconstructed_route | descriptive |  | 35 | 0.131 | 0.131 | 0 | no_skill_shown |  | 1 |
| T2_measured | av3322_T_263_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 35 | 0.0393 | 0.0461 | 0 | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | api_as_shipped | known | descriptive |  | 35 | 0.068 | 0.068 |  | beats_trivial |  | 1 |
| T2_measured | av3322_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 35 | 0.0588 | 0.0588 |  | beats_trivial |  | 1 |
| T2_measured | av3322_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 35 | 0.0404 | 0.0418 |  | beats_trivial |  |  |
| T2_measured | av3322_T_lt_263K | api_as_shipped | known | descriptive |  | 35 | -0.574 | 0.574 | 0.886 | worse_than_trivial |  |  |
| T2_measured | av3322_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 35 | 0.0319 | 0.0995 | 0 | beats_trivial |  | 1 |
| T2_measured | av3322_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 35 | 0.113 | 0.135 | 0 | beats_trivial |  |  |
| T2_measured | fuellib_T_ge_283K | api_as_shipped | known | descriptive |  | 4 | 0.0567 | 0.0567 |  | n_too_small |  | 1 |
| T2_measured | fuellib_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 4 | 0.0567 | 0.0567 |  | n_too_small |  |  |
| T2_measured | fuellib_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 3 | -0.0431 | 0.0431 |  | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | api_as_shipped | known | descriptive |  | 4 | -0.628 | 0.628 | 1 | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 4 | -0.0226 | 0.0965 | 0 | n_too_small |  | 1 |
| T2_measured | fuellib_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 4 | -0.0544 | 0.0965 | 0 | n_too_small |  |  |
| T2_measured | pooled_T_263_283K | api_as_shipped | known | descriptive |  | 38 | -0.00593 | 0.0327 | 0 | beats_trivial |  |  |
| T2_measured | pooled_T_263_283K | api_as_shipped | reconstructed_route | descriptive |  | 38 | 0.128 | 0.128 | 0 | no_skill_shown |  | 1 |
| T2_measured | pooled_T_263_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 38 | 0.0331 | 0.0455 | 0 | beats_trivial |  |  |
| T2_measured | pooled_T_ge_283K | api_as_shipped | known | descriptive |  | 54 | 0.045 | 0.0626 |  | beats_trivial |  | 0.981 |
| T2_measured | pooled_T_ge_283K | api_as_shipped | reconstructed_route | descriptive |  | 54 | 0.0386 | 0.0563 |  | beats_trivial |  | 1 |
| T2_measured | pooled_T_ge_283K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 53 | -0.00834 | 0.0626 |  | beats_trivial |  |  |
| T2_measured | pooled_T_lt_263K | api_as_shipped | known | descriptive |  | 65 | -0.504 | 0.527 | 0.692 | worse_than_trivial |  |  |
| T2_measured | pooled_T_lt_263K | api_as_shipped | reconstructed_route | descriptive |  | 65 | 0.0348 | 0.111 | 0.0154 | beats_trivial |  | 0.969 |
| T2_measured | pooled_T_lt_263K | refutas_raw+k1 | k1_cold_offset_loso | descriptive |  | 65 | 0.0341 | 0.138 | 0.0154 | beats_trivial |  |  |

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
| T2_measured | ajftd_gcxgc_T_263_283K | ideal_volume | designated | descriptive |  | 3 | -6.94 | 6.94 |  | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_ge_283K | ideal_volume | designated | descriptive |  | 32 | -4.57 | 6.81 |  | beats_trivial |  |  |
| T2_measured | ajftd_gcxgc_T_lt_263K | ideal_volume | designated | descriptive |  | 3 | -6.35 | 6.35 |  | n_too_small |  |  |
| T2_measured | av3322_T_263_283K | ideal_volume | designated | descriptive |  | 35 | 0.23 | 1.67 |  | beats_trivial |  |  |
| T2_measured | av3322_T_ge_283K | ideal_volume | designated | descriptive |  | 35 | -0.103 | 1.55 |  | beats_trivial |  |  |
| T2_measured | av3322_T_lt_263K | ideal_volume | designated | descriptive |  | 35 | 1.18 | 1.99 |  | beats_trivial |  |  |
| T2_measured | fuellib_T_263_283K | ideal_volume | designated | descriptive |  | 4 | -6.12 | 6.12 |  | n_too_small |  |  |
| T2_measured | fuellib_T_ge_283K | ideal_volume | designated | descriptive |  | 4 | -7.06 | 7.06 |  | n_too_small |  |  |
| T2_measured | fuellib_T_lt_263K | ideal_volume | designated | descriptive |  | 2 | -6.93 | 6.93 |  | n_too_small |  |  |
| T2_measured | pooled_T_263_283K | ideal_volume | designated | descriptive |  | 41 | -0.716 | 2.34 |  | beats_trivial |  |  |
| T2_measured | pooled_T_ge_283K | ideal_volume | designated | descriptive |  | 69 | -2.38 | 4.15 |  | beats_trivial |  |  |
| T2_measured | pooled_T_lt_263K | ideal_volume | designated | descriptive |  | 39 | 0.402 | 2.44 |  | beats_trivial |  |  |
| diesel | face_diesel | ideal_volume | designated | descriptive |  | 4 | 5.84 | 15 |  | n_too_small |  |  |
| T2_measured | ajftd_gcxgc_T_263_283K | api_as_shipped | known | descriptive |  | 3 | -6.94 | 6.94 |  | n_too_small |  | 0 |
| T2_measured | ajftd_gcxgc_T_ge_283K | api_as_shipped | known | descriptive |  | 32 | -4.57 | 6.81 |  | beats_trivial |  | 0.188 |
| T2_measured | ajftd_gcxgc_T_lt_263K | api_as_shipped | known | descriptive |  | 3 | -6.35 | 6.35 |  | n_too_small |  | 0 |
| T2_measured | av3322_T_263_283K | api_as_shipped | known | descriptive |  | 35 | 0.23 | 1.67 |  | beats_trivial |  | 0.971 |
| T2_measured | av3322_T_ge_283K | api_as_shipped | known | descriptive |  | 35 | -0.103 | 1.55 |  | beats_trivial |  | 0.943 |
| T2_measured | av3322_T_lt_263K | api_as_shipped | known | descriptive |  | 35 | 1.18 | 1.99 |  | beats_trivial |  | 0.943 |
| T2_measured | fuellib_T_263_283K | api_as_shipped | known | descriptive |  | 4 | -6.12 | 6.12 |  | n_too_small |  | 0.25 |
| T2_measured | fuellib_T_ge_283K | api_as_shipped | known | descriptive |  | 4 | -7.06 | 7.06 |  | n_too_small |  | 0.25 |
| T2_measured | fuellib_T_lt_263K | api_as_shipped | known | descriptive |  | 2 | -6.93 | 6.93 |  | n_too_small |  | 0 |
| T2_measured | pooled_T_263_283K | api_as_shipped | known | descriptive |  | 41 | -0.716 | 2.34 |  | beats_trivial |  | 0.854 |
| T2_measured | pooled_T_ge_283K | api_as_shipped | known | descriptive |  | 69 | -2.38 | 4.15 |  | beats_trivial |  | 0.58 |
| T2_measured | pooled_T_lt_263K | api_as_shipped | known | descriptive |  | 39 | 0.402 | 2.44 |  | beats_trivial |  | 0.846 |

## flash_point_C

Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.

| tier | set | model | variant | label | T_K | n_groups | bias | mae | unsafe_rate | verdict | margin90 | coverage90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CRC647_r3 | crc647_r3 | ideal_le_chatelier | designated | dev |  | 45 | 0.436 | 4.03 | 0.133 | beats_trivial |  |  |
| CRC647_r3 | crc647_r3 | ideal_le_chatelier | designated | fold_b |  | 12 | 2.46 | 3.72 | 0.25 | beats_trivial |  |  |
| T1_lab | lab_mixtures_31_systems | ideal_le_chatelier | designated | dev |  | 31 | -3.3 | 3.74 | 0.00794 | not_applicable | 3.7 | 0.935 |
| T2_measured | ajftd_gcxgc | api_as_shipped | designated | descriptive |  | 31 | 0.00194 | 2.92 | 0.0968 | beats_trivial | 2.21 | 0.839 |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | designated | descriptive |  | 31 | -0.152 | 3 | 0.0968 | beats_trivial | 3.7 | 0.871 |
| T2_measured | av3322 | api_as_shipped | designated | descriptive |  | 35 | -0.112 | 1.24 | 0 | beats_trivial | 2.21 | 0.971 |
| T2_measured | av3322 | ideal_le_chatelier | designated | descriptive |  | 35 | -0.0275 | 1.2 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | designated | descriptive |  | 28 | 0.0904 | 1.29 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | designated | descriptive |  | 32 | -0.0446 | 2.88 | 0.0938 | beats_trivial | 3.7 | 0.875 |
| T2_measured | family_other | ideal_le_chatelier | designated | descriptive |  | 10 | 0.3 | 2.57 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | fuellib | api_as_shipped | designated | descriptive |  | 4 | 2.64 | 4.78 | 0.25 | n_too_small | 2.21 | 0.75 |
| T2_measured | fuellib | ideal_le_chatelier | designated | descriptive |  | 4 | 2.44 | 4.73 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | pooled | api_as_shipped | designated | descriptive |  | 70 | 0.0961 | 2.19 | 0.0571 | beats_trivial | 2.21 | 0.9 |
| T2_measured | pooled | ideal_le_chatelier | designated | descriptive |  | 70 | 0.0585 | 2.2 | 0.0571 | beats_trivial | 3.7 | 0.929 |
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
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | lfl_overlay | descriptive |  | 31 | 0.00194 | 2.92 | 0.0968 | beats_trivial | 3.7 | 0.871 |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 31 | -0.371 | 3.27 | 0.0968 | no_skill_shown | 3.7 | 0.903 |
| T2_measured | ajftd_gcxgc | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 31 | -0.523 | 3.34 | 0.0968 | no_skill_shown | 3.7 | 0.903 |
| T2_measured | av3322 | ideal_le_chatelier | lfl_overlay | descriptive |  | 35 | -0.112 | 1.24 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | av3322 | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 35 | -0.46 | 1.22 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | av3322 | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 35 | -0.378 | 1.18 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | lfl_overlay | descriptive |  | 28 | 0.013 | 1.31 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 28 | -0.297 | 1.25 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D56 | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 28 | -0.221 | 1.21 | 0 | beats_trivial | 3.7 | 1 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | lfl_overlay | descriptive |  | 32 | 0.0688 | 2.84 | 0.0938 | beats_trivial | 3.7 | 0.875 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 32 | -0.329 | 3.13 | 0.0938 | beats_trivial | 3.7 | 0.906 |
| T2_measured | family_D93_or_IP170 | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 32 | -0.441 | 3.18 | 0.0938 | beats_trivial | 3.7 | 0.906 |
| T2_measured | family_other | ideal_le_chatelier | lfl_overlay | descriptive |  | 10 | 0.417 | 2.55 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | family_other | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 10 | 0.0672 | 2.9 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | family_other | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 10 | -0.0506 | 2.92 | 0.1 | no_skill_shown | 3.7 | 0.9 |
| T2_measured | fuellib | ideal_le_chatelier | lfl_overlay | descriptive |  | 4 | 2.64 | 4.78 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | fuellib | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 4 | 2.36 | 5.07 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | fuellib | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 4 | 2.16 | 5.01 | 0.25 | n_too_small | 3.7 | 0.75 |
| T2_measured | pooled | ideal_le_chatelier | lfl_overlay | descriptive |  | 70 | 0.0961 | 2.19 | 0.0571 | beats_trivial | 3.7 | 0.929 |
| T2_measured | pooled | ideal_le_chatelier | lfl_overlay+mesitylene_fixed | descriptive |  | 70 | -0.26 | 2.35 | 0.0571 | beats_trivial | 3.7 | 0.943 |
| T2_measured | pooled | ideal_le_chatelier | mesitylene_fixed | descriptive |  | 70 | -0.297 | 2.36 | 0.0571 | beats_trivial | 3.7 | 0.943 |
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
| T2_measured | ajftd_gcxgc_measured | ga_comp_H | designated | descriptive |  | 28 | 0.358 | 0.36 |  | no_skill_shown |  |  |
| T2_measured | av3322 | ga_comp_H | designated | descriptive |  | 35 | 0.131 | 0.131 |  | no_skill_shown |  |  |
| T2_measured | measured_pooled | ga_comp_H | designated | descriptive |  | 28 | 0.358 | 0.36 |  | no_skill_shown |  |  |
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
| T2_measured | ajftd_gcxgc_measured | ga_const_H | const_H | descriptive |  | 28 | 0.564 | 0.564 |  | worse_than_trivial |  |  |
| T2_measured | av3322 | ga_const_H | const_H | descriptive |  | 35 | 0.143 | 0.143 |  | worse_than_trivial |  |  |
| T2_measured | measured_pooled | ga_const_H | const_H | descriptive |  | 28 | 0.564 | 0.564 |  | worse_than_trivial |  |  |

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
