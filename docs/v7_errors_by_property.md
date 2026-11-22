# v7 errors by property and method (dev rows, out-of-fold; 2026-11-22)

MAE in the property's unit; lower is better. Best = lowest point MAE among the methods listed, not a significance claim: where the interval of the gain includes 0 the winner is only nominal (stated in the last column). Dev rows were seen in v5 or v6 and tuned on, so these are exploratory. Sources: `results/v7/track_b/{retrain,descriptor_gbm,curves_v1_long,chemprop_parity,param_heads}/`, `retrain_v5/`, `set_models/`, `mixture/`, `flash_oxy/`, `freezing/`, `flash/hanna/`, `results/v5/BASE/rows_t2_flash.csv`.

## A. Track B pure-compound GNN set (scaffold GroupKFold)
| Property (unit, n) | Constant | v1 D-MPNN | v7 D-MPNN (nested Optuna) | Long-cap D-MPNN | LightGBM | Random forest | Best | Is the win real? |
|---|---|---|---|---|---|---|---|---|
| Flash point, flash_C (degC, 837) | 55.22 | 19.11 | 17.22 | 18.09 | 12.25 | 11.61 | RF | Yes: trees beat the GNN (Holm p 0.010); Chemprop 16.83 (Task 7 parity run) |
| Cetane (259) | 18.56 | 9.04 | 8.64 | 8.66 | 9.31 | 8.57 | RF (v7 GNN 8.64) | No: GNN vs trees CI includes 0 |
| Derived cetane, dcn (89) | 18.37 | 5.56 | 5.40 | 5.34 | 8.09 | 7.09 | Long-cap GNN (v7 5.40) | Yes: GNN beats trees (+1.70, Holm p 0.044) |
| Density (kg/m3, 558) | 107.31 | 21.90 | 20.90 | 22.96 | 24.77 | 33.60 | v7 GNN | No: CI includes 0 |
| Heat of fusion, dhfus (kJ/mol, 81) | 19.14 | 12.80 | 10.98 | 11.35 | 13.17 | 7.23 | RF | Yes: RF beats GNN (Holm p 0.032) |
| Log viscosity at 293 K, logvisc (180) | 0.457 | 0.225 | 0.223 | 0.223 | 0.175 | 0.148 | RF | Yes: trees beat GNN (Holm p 0.056) |
| MON (363) | 14.53 | 8.24 | 7.95 | 7.76 | 9.02 | 9.11 | Long-cap GNN (v7 7.95) | No |
| RON (378) | 16.58 | 8.99 | 8.81 | 8.38 | 10.73 | 9.88 | Long-cap GNN | No (GNN vs trees CI includes 0 after Holm) |
| Melting point, tm_K (K, 547) | 71.82 | 30.42 | 27.22 | 28.01 | 28.33 | 27.38 | v7 GNN | Vs v1 yes (+3.20 K, Holm p 0.015); vs trees no |
| YSI (168) | 262.91 | 74.81 | 69.17 | 73.18 | 91.79 | 73.28 | v7 GNN | No (Holm p 1.0); Chemprop parity run 121.6 vs in-house 125.4 under its own protocol |

## B. Flash point of blends and fuels
| Row set | Methods and MAE (degC) | Best | Is the win real? |
|---|---|---|---|
| Real fuels, T2, 68 dev rows | v5 designated 2.24; LFL overlay 2.26; mesitylene fixed 2.33; both 2.34 (over-prediction: 1.34, 1.36, 1.23, 1.25) | v5 designated | Differences are tiny and not tested in v7 |
| Hydrocarbon binaries, Leenhouts (92 rows, 16 systems) | ideal 3.20; HANNA 3.66; UNIFAC 3.81 | Ideal (production) | Yes: HANNA worse than ideal (-0.48 [-0.75, -0.24]) |
| Oxygenate mixtures (641 rows, 57 systems) | desc_gbm 4.08; v6 MolPool v1 3.44; v7 MolPool 3.61 | v6 MolPool v1 | No: v7 vs v1 -0.17 [-0.65, 0.31]; v7 vs desc_gbm +0.47, Holm p 0.162 |
| Pure compounds from predicted Antoine constants (53 rows) | pooled curve 102.8; D-MPNN head 22.1; registry refit constants 9.3 (reference) | Registry constants | Head over-predicts by 16.4 degC on average: unusable |

## C. Freezing point (degC, 52 dev fuels)
| Method | LOGO MAE | LOFO MAE | Note |
|---|---|---|---|
| Constant | 6.64 | 6.99 | |
| Ridge | 5.83 | 5.23 | |
| UNIFAC-SvL baseline | 8.83 | 8.83 | worse than constant |
| Limiting-species baseline | 9.23 | 9.23 | worse than constant |
| Heavy-tail GP | 6.01 | 6.70 | fails the skill rule |
| v6 LGBM residual on isomer tiles | 3.31 | n/a | |
| v7 monotone LGBM refit (nested Optuna) | 3.07 | n/a | vs refit original +0.25 [-0.11, 0.61] |
| MolPool fuel-only | 3.07 | 2.97 | |
| MolPool joint (pure Tm plus fuel) | **3.00** | **2.80** | vs constant +3.63 [2.08, 5.26]; vs fuel-only +0.07 [-0.21, 0.34] |

Best: joint MolPool, nominally; it is statistically tied with fuel-only MolPool and the LGBM residual. Dev only, so real-fuel freezing stays "status: no skill shown".

## D. Mixture properties (v5 and v6 Track C dev rows)
| Property (unit) | Physics baseline | v5 pair LGBM (refit) | v7 pair LGBM | Set models v7 (pds / gnn / st) | v7 MolPool | Best | Is the win real? |
|---|---|---|---|---|---|---|---|
| Viscosity, ln, 264 groups | 0.0824 | 0.0599 | 0.0604 | 0.0608 / 0.0623 / 0.0577 | n/a | Set transformer 0.0577 | No: gain +0.0022 [-0.0002, 0.0047], gate fails |
| Density (kg/m3), 167 groups | 1.334 | 0.958 | 0.987 | 0.943 / 1.118 / 0.994 | n/a | Deep-sets pds 0.943 | No (CI includes 0); GNN significantly worse |
| Boehm density (kg/m3) | 2.873 | n/a | n/a | n/a | 2.947 | Physics baseline | MolPool vs physics -0.07, CI includes 0 |
| Boehm ln viscosity | 0.0551 (desc_gbm) | n/a | n/a | n/a | 0.0593 | desc_gbm | MolPool beats v6 v1 (+0.019) but not desc_gbm |
| KAUST MON | 2.912 (lin_molar) | n/a | n/a | n/a | 2.697 | v7 MolPool | No: +0.22 [-0.46, 0.37] |
| KAUST RON | 1.486 (desc_gbm) | n/a | n/a | n/a | 2.022 | desc_gbm | Descriptor GBM better |
| Leenhouts MON | 3.084 (desc_gbm) | n/a | n/a | n/a | 3.695 | desc_gbm | Descriptor GBM better; v6 v1 2.294 is better still |

## E. Equation parameters and other
| Quantity | Methods and MAE | Best | Is the win real? |
|---|---|---|---|
| log10 vapour pressure (189 molecules) | pooled curve 1.793; D-MPNN Antoine head 0.474; descriptor ridge 0.220; per-molecule refit 0.005 (uses own labels) | Ridge | Yes: head worse than ridge (-0.34 [-0.40, -0.28]) |
| ln viscosity from Andrade (136 molecules) | pooled 0.726; head 0.280; ridge 0.284; refit 0.033 | Head, nominally | No: -0.038 [-0.079, 0.002] |
| Pure melting point, fusion_estimate (K, LOFO, 109 compounds) | physics 60.12; refit original 54.69; v7 50.48 | v7 | No: -4.9 [-16.4, 5.2] |
| YSI on 169 dev compounds (ysi_lgbm) | physics 234.1; refit original 236.8; v7 tree 220.1 | v7 tree | Yes: +35.9 [11.9, 64.2] (different row set from table A) |

## F. Coverage of the one-sided 90 % margin (Task 19)
Global margin vs family-conditional (Mondrian) coverage on dev rows: YSI poly_arom 0.28 to 0.96, flash D93/IP170 0.84 to 0.87 or 0.90, freezing crc_jeta 0.83 to 0.96 (`results/v7/uq/READING.md`). Mondrian margins are wide where errors are heavy tailed (YSI global margins about 300 to 700 units).

## G. Oxygenate results (added 2026-11-22; same OOF predictions, split by oxygenate vs hydrocarbon rows)
Only these properties have oxygenate labels: the production unified data (`data/unified/`) is hydrocarbon-only (`is_hc` True on every row), so density, viscosity, freezing, melting point, YSI, cetane, dcn, dhfus and vapour pressure have no oxygenate rows to score. Oxygenates enter through the Track B pool (MON, RON), the Track C mixtures (KAUST and Leenhouts MON/RON) and the separate flash-oxygenate set. Numbers are MAE from existing OOF files, dev rows only; "n/a" means that comparator was not computed for the head.

**Track B pure compounds (MON, RON)**
| Property | Rows | Constant | v7 D-MPNN | Long-cap D-MPNN | LightGBM | Random forest | Best |
|---|---|---|---|---|---|---|---|
| MON, oxygenates | 65 | 10.78 | 6.68 | 6.51 | 6.70 | 6.23 | RF |
| MON, hydrocarbons | 298 | 15.04 | 8.22 | 8.03 | 9.53 | 9.74 | Long-cap D-MPNN |
| RON, oxygenates | 79 | 9.04 | 6.81 | 6.92 | 7.94 | 6.96 | v7 D-MPNN |
| RON, hydrocarbons | 299 | 18.52 | 9.34 | 8.76 | 11.47 | 10.65 | Long-cap D-MPNN |

**Track C mixtures (rows containing at least one non-C/H atom vs hydrocarbon-only rows)**
| Head, property | Subset | Rows (groups) | Linear molar | Desc. GBM | v6 MolPool v1 | v7 MolPool | Best |
|---|---|---|---|---|---|---|---|
| KAUST MON | oxygenate | 39 (4) | 5.34 | n/a | 6.37 | 4.78 | v7 MolPool |
| KAUST MON | hydrocarbon | 51 (3) | 1.06 | n/a | 1.03 | 1.10 | v6 MolPool v1 |
| KAUST RON | oxygenate | 39 (4) | n/a | 2.29 | 4.44 | 2.90 | Desc. GBM |
| KAUST RON | hydrocarbon | 51 (3) | n/a | 0.87 | 1.20 | 1.36 | Desc. GBM |
| Leenhouts MON | oxygenate | 106 (21) | n/a | 3.23 | 2.71 | 3.99 | v6 MolPool v1 |
| Leenhouts MON | hydrocarbon | 75 (13) | n/a | 2.88 | 1.71 | 3.28 | v6 MolPool v1 |

**Flash point, oxygenate set (641 rows, 57 systems; Task 11)**: desc_gbm 4.08, v6 MolPool v1 3.44, v7 MolPool 3.61; best v6 v1. On the 95 hydrocarbon rows of the same set, a model trained on hydrocarbons only (arm A) has MAE 4.21 and the model trained with oxygenates (arm B) 6.40, so adding oxygenates hurts hydrocarbon rows (non-inferiority fails). Le Chatelier is not evaluable on oxygenates (no Antoine or LFL constants).

Reading: oxygenate MON and RON are predicted about as well as hydrocarbon ones for pure compounds (smaller MAE, but their constant baseline is also smaller, 10.8 and 9.0 vs 15.0 and 18.5). In mixtures the oxygenate rows have the larger errors, and the v7 MolPool is worse than v6 v1 on Leenhouts MON in both subsets. Counts per cell are small (3 to 21 groups), so these are descriptive, not significance tests. Adding oxygenate labels for the other properties needs new data acquisition, which is not part of v7.
