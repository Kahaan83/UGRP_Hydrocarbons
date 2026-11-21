# Phase 2 readings and Phase 3 trial budgets (dev rows only; all post-selection results exploratory)

Budget rule (predeclared in the plan): a model family whose Phase 2 evidence shows no headroom (Task 5 flat, Task 8 slope CI including 0, Task 6 comparator equal or better) gets the minimum budget, 20 trials per outer fold; all others get 60.

| Family | Task 5 | Task 6 | Task 7 and debug | Task 8 (0.8 to 1.0 slope CI) | Budget |
|---|---|---|---|---|---|
| Track B D-MPNN | not flat: median best epoch 136 (long cap) vs 36.5 (v1 cap); outer MAE lower with long cap on 8 of 10 targets | trees equal or better on 9 of 10 targets; GNN earns its place only on dcn under the v1 cap and on none after the long cap (Holm) | flash_C gap to Chemprop closes with epochs 1000, patience 100 (H2); stop rule lifted for Task 10 with that cap | CI includes 0 on 9 of 10 targets; ron excludes 0 (-3.99, -0.06). Only 2 seeds and 300 epochs, so wide | 60 (Task 5 not flat, so the rule's no-headroom condition is not met) |
| Descriptor/LightGBM trees (Task 13 trees) | n/a | the comparators | n/a | CI excludes 0 (negative) on all 10 targets: more data still helps | 60 |
| Freezing LOGO ridge and LGBM residual | n/a | n/a | n/a | CI includes 0 for both (n = 52 entities, wide) | 20 |
| MolPool mixture heads, flash oxygenate arm (Task 11) | n/a | n/a | n/a | no Phase 2 curve | 60 |
| v5 set models (Task 12) | n/a | n/a | n/a | no Phase 2 curve | 60 |
| v5 viscosity/density pair residuals (Task 13) | n/a | n/a | n/a | not run in Task 8 (see its predeclaration) | 60 |

Time projections (projected from Task 5 measurements, not measured for Phase 3): D-MPNN 0.17 to 0.20 s per epoch per model on the RTX 3050; an inner fit stops near best epoch plus patience (about 240 epochs, about 45 s), three inner folds per trial, so about 2 to 3 min per trial before pruning. 60 trials x 5 outer folds is about 10 to 15 h unpruned; `timeout_s` per outer fold is set to 5400 s (60 trials) and 1800 s (20 trials), so the study stops on time even if fewer trials finish. Trials completed are recorded in `best_params.csv`.

Caveats: Task 8 slopes rest on one fraction step and 15 units per target; flat is "no evidence of benefit", not saturation. Task 5 and long-cap comparisons are post-hoc.
