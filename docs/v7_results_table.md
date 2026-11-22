# v7 complete results table (2026-11-22)

All rows are dev rows already seen in v5 or v6 (the v7 searches were tuned on them), so every gain is exploratory. No locked set was scored except as noted, nothing is adopted into `predict_fuel`, and no v5 or v6 verdict is reopened. Gain = comparator minus model (positive favours the v7 model); intervals are 90 % paired group-bootstrap CIs. Every number is read from the committed file in the last column (paths under `results/v7/`).

## Phase 0 and 1: guardrails, instrumentation
| Task | What | Outcome | Source |
|---|---|---|---|
| 0.5 | Data inventory and sources survey | documents committed; acquisition list for the user | `docs/` data-source docs |
| 1 | Lock guard | refuses v6-scored rows; `guard.LOCKS` covers six locks | `src/fuelblend/v7/guard.py` |
| 2 | Curve logger and hooks | per-epoch CSV logging in torch and LightGBM loops; golden v6 outputs unchanged | `src/fuelblend/ml/curves.py` |
| 3 | GPU venv, device, monitor | RTX 3050 4 GB, torch cu128, strict device resolution | `scripts/v7/` setup |
| 4 | Runner skeleton | scaffold folds on full pool, YSI lock masked after folds | `src/fuelblend/v7/runner.py` |

## Phase 2: cheap diagnostics
| Task | Question | Result | Source |
|---|---|---|---|
| 5 | Is the v1 D-MPNN stopped too early? | Not flat: median best epoch 136 (long cap) vs 36.5 (v1 cap); outer MAE lower with long cap on 8 of 10 targets | `PHASE2_READINGS.md`, `track_b/curves_v1*` |
| 6 | Descriptor GBM and RF comparators | trees equal or better on 9 of 10 targets; GNN earns its place on dcn only (v1 cap), none after long cap | `track_b/descriptor_gbm/` |
| 7 | Chemprop parity | flash_C: in-house 17.19 vs Chemprop 16.83, CI [-0.75, -0.06] (stop rule triggered); ysi 125.4 vs 121.6, CI includes 0 | `track_b/chemprop_parity/READING.md` |
| 7b | flash_C debugging | gap closes with epochs 1000 and patience 100; stop rule lifted for Task 10 | `track_b/debug_flash_chemprop/` |
| 8 | Sample-size slopes | D-MPNN: CI includes 0 on 9 of 10 targets (ron excludes 0); trees: CI excludes 0 (negative) on all 10; freezing ridge and LGBM CI includes 0 | `PHASE2_READINGS.md`, `sample_size/` |

## Phase 3: instrumented retrain with nested HPO
| Task | Model | Comparator | Result | Source |
|---|---|---|---|---|
| 9 | Optuna harness | n/a | nested, resumable studies; budgets 60 trials per fold (20 for freezing) | `PHASE2_READINGS.md` |
| 10 | Track B D-MPNN, 10 targets | v1 D-MPNN; trees | beats v1 on tm_K only (+3.20 K, Holm p 0.015); beats trees on dcn only (+1.70, p 0.044); trees beat GNN on flash_C (RF 11.61 vs 17.22, p 0.010), dhfus (7.23 vs 10.98, p 0.032), logvisc293 (0.148 vs 0.223, p 0.056) | `track_b/retrain/READING.md` |
| 11a | MolPool mixture, 5 heads | best of physics, t_bias, desc_gbm, lin_molar; v1 | none beats its best comparator; beats v1 on boehm ln_visc (+0.019) and kaust ron (+0.58); worse than v1 on leenhouts mon (-1.40) | `mixture/READING.md` |
| 11b | Flash oxygenate arm | desc_gbm 4.08; v1 3.44 | +0.47 [0.03, 0.91] vs desc_gbm, Holm p 0.162; -0.17 vs v1; hydrocarbon non-inferiority fails | `flash_oxy/READING.md` |
| 12 | Set models (3 architectures, 2 properties) | LightGBM pair baselines | no cell passes the gate; best viscosity set transformer +0.0022 [-0.0002, 0.0047]; density GNN significantly worse | `set_models/READING.md` |
| 13a | viscosity_pair | refit original | -0.0005 [-0.0021, 0.0012] | `retrain_v5/READING.md` |
| 13b | density_pair | refit original | -0.0295 [-0.0839, 0.0238] | `retrain_v5/READING.md` |
| 13c | ysi_lgbm (169 dev compounds) | refit original | +35.9 [11.9, 64.2], the only tree gain | `retrain_v5/READING.md` |
| 13d | fusion_estimate | refit original | -4.9 [-16.4, 5.2] | `retrain_v5/READING.md` |
| 13e | freezing_step3 | refit original | +0.247 [-0.111, 0.607] | `retrain_v5/READING.md` |
| 13f | mixture_gbr, 5 units | refit original | 0 of 5 units above 0; boehm ln_visc significantly worse | `retrain_v5/READING.md` |

## Phase 4: literature methods (freezing LOGO or LOFO MAE in degC on 52 dev fuels unless stated)
| Task | Model | Result | Source |
|---|---|---|---|
| 14 | Limiting-species baseline | LOGO MAE 9.23 vs constant 6.64: worse than constant | `freezing/baseline_d/READING.md` |
| 15 | Heavy-tail ARD GP | LOGO 6.01 (constant 6.64, ridge 5.83); LOFO 6.70 (constant 6.99, ridge 5.23); fails the predeclared skill rule; novelty flag rests on 4 to 5 fuels | `freezing/gp/READING.md` |
| 16 | Joint pure plus fuel MolPool | LOGO 3.00, LOFO 2.80; beats constant (+3.63 [2.08, 5.26]; +4.19 [2.65, 5.82]) and ridge (+2.83; +2.43); vs fuel-only +0.07 [-0.21, 0.34] and +0.17 [-0.03, 0.38]; v6 LGBM residual 3.32 (not bootstrapped) | `freezing/molpool_joint/READING.md` |
| 17 | HANNA gamma, 92 flash rows, 16 systems | MAE ideal 3.20, UNIFAC 3.81, HANNA 3.66; HANNA vs ideal -0.48 [-0.75, -0.24]; vs UNIFAC -0.21 [-0.54, 0.10]; gate fails | `flash/hanna/READING.md` |
| 18 | Antoine and Andrade heads | log P MAE 0.474 vs pooled 1.793 and ridge 0.220 (vs ridge -0.34 [-0.40, -0.28]); ln mu 0.280 vs ridge 0.284 (-0.038 [-0.079, 0.002]); flash recompute MAE 22.1 vs registry 9.3, mean over-prediction 16.4 | `track_b/param_heads/READING.md` |
| 19 | Conformal by family | global margin under-covers YSI poly_arom (0.28 ridge, 0.24 LGBM), flash D93/IP170 (0.84), freezing crc_jeta (0.83); Mondrian restores 0.87 to 0.96 at large width; weighted variant over-covers (0.92 to 1.00) | `uq/READING.md` |
| 20 | Fresh lock | 61 unseen systems (658 rows) locked before any fit; primary comparison HANNA vs UNIFAC; unspent (Task 17 gives no dev signal) | `lock/LOCK.md`, `lock/overlap.json` |

## Phase 5: documentation
| Task | What | Source |
|---|---|---|
| 21 | v7 protocol and summary; v6 negative-results rows; literature year note | `docs/v7_protocol.md`, `docs/v7_summary.md`, `docs/v6_negative_results.md` |
| | Full suite: 844 passed, 2 skipped, 1 failed (the CRLF `test_committed_splits_match_rebuild`, unchanged by v7) | `docs/v7_summary.md` |

## Verdict summary
| Question | v7 answer |
|---|---|
| Does any v7 model beat a comparator on real fuels or locked rows? | Not tested; no lock scored. |
| Dev-only gains that survive the CI rule | tm_K (Task 10), dcn (GNN vs trees), boehm ln_visc and kaust ron (MolPool vs v1), ysi_lgbm tree, freezing joint MolPool vs constant and ridge |
| Real-fuel freezing point | status: no skill shown (dev MolPool gain is exploratory and needs a fresh lock) |
| Cold viscosity, flash from HANNA or predicted Antoine | no skill shown |
