# Task 7 reading (dev rows only; 5 folds x 5 seeds, v1 cap, 1171 s on the GPU)

flash_C (n=837): MAE in-house 17.19, Chemprop 16.83; gain (Chemprop minus in-house) 90 % CI [-0.75, -0.06]. The CI excludes 0 with in-house worse: predeclared stop rule triggered (difference 0.36 degC, about 2 %).
ysi (n=168): MAE in-house 125.4, Chemprop 121.6; CI [-16.0, 1.8] includes 0: parity holds.
Predeclared consequence: stop Phase 3 Track B tuning (Task 10) and open a debugging task on the in-house D-MPNN versus Chemprop. Task 10 is therefore not run until that task is resolved. Phase 3 tasks that do not tune Track B (11 to 13) are not covered by the rule.
