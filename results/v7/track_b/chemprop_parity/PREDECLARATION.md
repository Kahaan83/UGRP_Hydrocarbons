# Task 7 predeclaration: Chemprop parity

Rows/folds: Track B pool, Task 4 folds (full-pool generic-Murcko GroupKFold(5)), ysi_holdout masked. Targets: `flash_C`, `ysi`, single-task.
In-house: `fuelblend.ml.dmpnn` hidden 200, depth 4, dropout 0.1, v1 cap (150 epochs, patience 25), 5 seeds. Chemprop 2.3.1: BondMessagePassing d_h 200, depth 4, dropout 0.1, mean aggregation, FFN hidden 200 x 1, 150 epochs, early stopping patience 25 on val_loss with best checkpoint restored, standardised targets, 5 seeds, batch 32, lr defaults. Inner val: 15 % of train groups (same split for both).
Metric: OOF MAE; gain = MAE(chemprop) - MAE(in-house), paired group bootstrap 90 % CI (4000 draws).
Reading: parity holds if the CI contains 0 for both targets. If in-house is worse with the CI excluding 0 on either target, stop Phase 3 Track B tuning and open a debugging task.
