# Phase 3 readings (Tasks 10 to 13): v7 versus comparator, and the curve reading
All rows are dev rows already seen in v5 or v6, all searches were tuned on them, so every gain here is exploratory and none reopens a v5 or v6 verdict. Locked sets were never scored. Nothing was adopted into `predict_fuel()`. Details and tables: the READING files named in each row.

| task / model | comparator | gain vs comparator (90 % CI) | curve reading |
|---|---|---|---|
| 10 Track B D-MPNN (10 targets), `track_b/retrain/READING.md` | v1 D-MPNN; descriptor trees | beats v1 on tm_K only (+3.20 K, Holm p 0.015); beats trees on dcn only (+1.70, Holm p 0.044); trees beat the GNN on flash_C, dhfus, logvisc293 | multi-task final fits: best epoch median 198 (range 98 to 623) of a median 298 run (patience stop); val/train loss 2.5 at the best epoch, an overfitting gap |
| 11 MolPool mixture (5 heads), `mixture/READING.md` | best of physics/t_bias/desc_gbm/lin_molar; v6 MolPool v1 | none beats its best comparator; beats v1 on boehm ln_visc (+0.019) and kaust ron (+0.58); worse than v1 on leenhouts mon (-1.40) | per-head curves in `mixture/retrain/*/curves_*.csv` (not summarised here) |
| 11 flash oxygenates, `flash_oxy/READING.md` | desc_gbm (4.08), v1 (3.44) | +0.47 [0.03, 0.91] vs desc_gbm, Holm p 0.162 (not significant); -0.17 vs v1; hydrocarbon non-inferiority fails | curves in `flash_oxy/retrain/arm*/` |
| 12 set models, 6 cells, `set_models/READING.md` | LightGBM pair baselines | no cell passes the gate; best viscosity st +0.0022 [-0.0002, 0.0047]; density gnn significantly worse | best epoch 4 to 27 (of up to 500, median stop 54 to 77); val/train 0.5 (pds, gnn viscosity) to 2.0 (st viscosity), 1.3 to 1.5 (density) |
| 13 viscosity_pair | refit original | -0.0005 [-0.0021, 0.0012] | best iteration 28 (stop 127) |
| 13 density_pair | refit original | -0.0295 [-0.0839, 0.0238] | best iteration 67 (stop 167) |
| 13 ysi_lgbm | refit original (comparison on 169 dev compounds) | +35.9 [11.9, 64.2] | best iteration 32 (stop 131) |
| 13 fusion_estimate | refit original | -4.9 [-16.4, 5.2] (family-macro) | best iteration 97 (stop 197) |
| 13 freezing_step3 | refit original | +0.247 [-0.111, 0.607] | best iteration 6 (stop 105) |
| 13 mixture_gbr (5 units) | refit original | 0 of 5 units above 0; boehm ln_visc significantly worse | best iteration 1 to 464 by unit |

Summary: after tuning with curves and nested search, the only intervals above 0 against a comparator are tm_K (Task 10 vs v1), dcn (GNN vs trees), boehm ln_visc and kaust ron (MolPool vs v1), and the ysi_lgbm tree (whose comparator is not the shipped ridge). Trees remain better than the GNN on three targets. Best epochs are early everywhere and the validation curves rise afterwards on most targets (val/train above 1), which supports stopping on a held-out split rather than long fixed budgets. No model is eligible for adoption on this evidence.

Full test suite at the Phase 3 exit: 835 passed, 2 skipped, 1 deselected (`tests/test_s4_splits.py::test_committed_splits_match_rebuild`). That test compares the text of `splits/locked.json` byte for byte with a rebuild; the working copy has CRLF line endings (the file was already modified in the working tree before this work and its git diff is line-ending only), so the comparison fails for a reason unrelated to v7. The file was not touched. Earlier baseline: 741 passed, 11 skipped.
