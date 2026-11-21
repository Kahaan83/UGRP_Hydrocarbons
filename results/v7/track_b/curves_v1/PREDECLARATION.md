# Task 5 predeclaration: epoch curves of the declared v1 D-MPNN config

Rows: v6 Track B pool (`results/v6/gnn/pool.csv`), folds from the full pool (generic-Murcko GroupKFold(5)), then ysi_holdout labels masked (dev only).
Configs: both v1 `CFGS` from `scripts/v6/gnn/train_eval.py` (h200 d4 do0.1; h300 d3 do0.2), 5 seeds, lr 1e-3, wd 1e-5, bs 32, no scheduler.
Regimes: v1-cap (epochs 150, patience 25, as v6) and long-cap (epochs 1000, patience 100). Config chosen per outer fold on the inner group-held-out validation loss (15 % of train groups), as v6.
Metric: outer-fold MAE per task (OOF over all folds) and the per-model `best_epoch_`.
Reading rule: if the median `best_epoch_` of the long-cap run is below 0.5 x 150 = 75 and the val loss is flat or rising afterwards, longer training is not the limit. Outer MAE of v1-cap and long-cap are reported side by side; no claim of gain without both.
