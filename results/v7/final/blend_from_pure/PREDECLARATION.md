# Predeclaration: blend properties from pure-compound predictions (exploratory, dev rows only)

Question: can the deployed pure-compound models (`fuelblend.api.predict_pure`) give a usable blend value through a linear mole-fraction rule, for properties where no blend model exists today?

Estimator: y_hat(mixture) = sum_i x_i p_i, x = mole fractions, p_i = deployed pure prediction for component i (derived cetane, MON, RON: D-MPNN long cap; per `best_methods.PURE`). No fitting on mixtures. p_i is the cross-fitted out-of-fold Track B prediction when the component is in the Track B pool (so a component never predicts itself from its own label); otherwise the final fitted model. Rows with any component that cannot be predicted (unparseable SMILES) are dropped and counted.

Tasks (dev rows; rows in `mixture_lock`, `mixture_holdout51` or `splits_locked` removed; mixture rows from `results/v6/mixture/rows.csv` with role == dev):
1. dcn_leenhouts: derived cetane, Leenhouts mixtures in `measurements_dedup.csv` (binary and multi).
2. mon_leenhouts: MON, Leenhouts mixtures (rows.csv).
3. mon_kaust: MON, KAUST mixtures (rows.csv).
4. ron_kaust: RON, KAUST mixtures (rows.csv).
Cetane has no mixture data and cannot be tested; it is not added to blends.

Comparators: C0 = leave-one-system-out median of the task's mixture values (a constant). Reference only, not a gate: the same linear rule with measured pure labels where every component has one.
Metric: MAE. Gain = MAE(C0) - MAE(model). Paired group bootstrap over systems (4000 draws, seed 0), 90 percent interval. Holm correction across the 4 tasks at 0.10.
Reading rule: a task is "usable, exploratory" only if the Holm-adjusted interval excludes 0 with gain > 0. Otherwise "no skill shown" and the property is not added to `predict_fuel`.
Consequence rule: usable tasks are added to `predict_fuel` as extra keys, labelled exploratory, uncertified, no lo/hi, with the task's MAE as evidence. No existing key, value or margin changes. Dev-only results cannot reopen a v5 or v6 verdict. Locked sets are not scored.
