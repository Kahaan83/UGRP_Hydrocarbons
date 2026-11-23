# v11 predeclaration: tuning the LightGBM pair models

Written before any tuned result exists. Question: do tuned LightGBM settings beat the fixed v5 settings (400 trees, lr 0.05, 15 leaves, min child 40, colsample 0.8, lambda 1) for the viscosity pair residual?

* Data: the v5 dev systems only (`role == "dev"`). Locked systems are not read, scored or inspected.
* Search: 30 random configurations (seed 11) over n_estimators {200,400,800}, learning_rate {0.02,0.05,0.1}, num_leaves {7,15,31}, min_child_samples {20,40,80}, colsample_bytree {0.6,0.8,1.0}, reg_lambda {1,5,20}. Score = out-of-fold per-system macro MAE of the corrected ln viscosity on dev rows at T >= 283.15 K, GroupKFold(5) on system_id (the v5 split).
* Confirmation (to remove selection optimism): the best configuration and the v5 default are re-scored on 3 other GroupKFold(5) partitions of the dev systems (systems shuffled with seeds 1, 2, 3).
* Adopt only if, pooled over the 3 confirmation partitions, the per-system MAE gain (default minus tuned) has a 95 % system-bootstrap interval above 0 and the mean gain is at least 0.002 ln units. Otherwise v5 settings stay and the result is a negative result in `docs/v5_negative_results.md`.
* If adopted: retrain, re-run the R-MLV gates, re-score the scoreboard, margins and designs, rewrite BASE and the manifest (user decision 2026-11-23), update docs and tests. If not adopted: BASE is unchanged.
* Density pair model: same protocol, run only if the viscosity result is adopted or if the user asks.
