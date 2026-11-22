# Task 12 reading: v5 set models retrained with validation split, early stopping and curves
Protocol: `PREDECLARATION.md`. Dev rows and v5 folds only; locked systems not scored. Outputs: `readings_dev.json`, `oof_corrections.csv.gz`, `curves.csv`, `run_log`-equivalent `scratch/v7_t12.log`. Gain = pair-baseline MAE minus candidate MAE (system-macro, paired system bootstrap, 90 % CI); gate passes only if the interval is wholly above 0.

| property / model | pair MAE | v7 MAE | gain 90 % CI | gate | v5 (CPU, 60 ep) MAE, gate |
|---|---|---|---|---|---|
| viscosity pds | 0.05989 | 0.06083 | [-0.0034, 0.0017] | fail | 0.0570, fail |
| viscosity gnn | 0.05989 | 0.06234 | [-0.0051, 0.0003] | fail | 0.0611, fail |
| viscosity st | 0.05989 | 0.05773 | [-0.0002, 0.0047] | fail | 0.0562, pass |
| density pds | 0.9578 | 0.9427 | [-0.066, 0.104] | fail | 0.9866, fail |
| density gnn | 0.9578 | 1.118 | [-0.282, -0.044] | fail (significant loss) | 1.0568, fail |
| density st | 0.9578 | 0.9936 | [-0.136, 0.062] | fail | 0.9418, fail |

Curves (median over 5 folds x 5 seeds; val loss is the weighted 20 % held-out-system loss in scaled units):
- Best-epoch location is early: viscosity pds 4, gnn 7, st 27; density pds 5, gnn 14, st 8. Median stop epoch 54 to 77 (patience 50). No fit reached 500 epochs, so the longer budget was not used; the cosine schedule over 500 epochs keeps the learning rate high at these early epochs, so this is a statement about this schedule, not about the optimum.
- Val/train loss at the best epoch: viscosity pds 0.48, gnn 0.56 (val below train, group-held-out systems are easier than the training mix); viscosity st 2.04 and density 1.3 to 1.5 (overfitting gap).

Status:
- The viscosity set transformer passed the v5 CPU gate (interval [0.0012, 0.0064]) but does not pass it in v7 (interval [-0.0002, 0.0047], lower bound just below 0), because 20 % of the training systems are held out for stopping and the interval moved; point gain 0.0022 vs 0.0037. This is a loss of the one v5 pass, not an improvement, and it comes with fewer training rows (stated in the predeclaration).
- No candidate passes the gate in v7. Density: no model beats the pair baseline; the GNN is significantly worse. Nothing is eligible for designation, nothing is adopted into `predict_fuel()`.
