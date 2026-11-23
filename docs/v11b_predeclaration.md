# v11b predeclaration: nested tuning of the viscosity and density pair models

Written before any v11b result exists. v11 (30 random configurations, reshuffled-partition check) stays recorded as the first attempt: gain 0.0010 ln, interval includes 0, not adopted. v11b is a second, stricter analysis.

**Correction to an earlier statement.** v9 added 4 density binary systems (209 to 213, `docs/final_numbers_v9.md` section A) and no viscosity systems. The v5 density delta regenerates its rows from the v4 committed table and `data/processed/v2/thermoml_clean_hc.csv`, not from `data/unified_v9/`; the 4 new systems are not in the density dev pool here (4 of about 2 hundred systems; adding them needs a new row builder and is not done). This is a stated limitation.

**Noise floor.** Replicate scatter in `measurements_dedup.csv` (binary and multi rows with n_obs >= 2): viscosity mean half-range 0.0007 ln (5,331 of 16,662 rows have replicates), density 0.11 kg/m3 (386 of 13,644). Both are far below the model errors (0.060 ln; about 1.35 kg/m3), so measurement noise does not rule tuning out.

## Method
* Data: v5 dev systems only; locked systems not read.
* Outer: GroupKFold(5) on system_id (the v5 folds). Inner: GroupKFold(3) inside each outer training split. Only outer test folds are scored.
* Tuner: Optuna TPE, 60 trials per outer fold, seed 11. Search: learning rate fixed 0.03; n_estimators 100 to 1200; num_leaves 4 to 48; min_child_samples 10 to 100; colsample 0.4 to 1; subsample 0.6 to 1; reg_lambda 0.1 to 50; loss in {l2, huber (delta 0.3 to 3), l1}; pair weight exponent p in {0, 1, 2} (w = (x_i x_j)^p, normalised per system; p = 2 is v5).
* Inner objective: per-system macro MAE of the corrected quantity on inner validation rows (viscosity: dev rows with T >= 283.15 K, ln units; density: dev rows, kg/m3).
* Baseline: the v5 settings on the same outer folds.

## Adoption rule (per model)
1. Per-system macro MAE gain (v5 minus tuned) on outer folds: mean at least 0.002 ln (viscosity) or 0.04 kg/m3 (density), and the one-sided system-bootstrap p-value (2000 resamples) below 0.05 after Holm correction across the two models.
2. Guard (viscosity): the conformal margins of the three temperature bins (R-MLV-6 definition) are not larger than with v5 settings.
3. Density: dev ternary rows are reported separately. Even if adopted, density stays a labelled variant unless the ternary failure clears.

## Consequences
If a model passes: retrain, re-run gates (a locked reading is then a second look, not a fresh hold-out, because R-MLV-2 and R-MLD-2 already read those systems), rebuild scoreboard, margins and designs, rewrite BASE and the manifest, update tests and docs. If neither passes: BASE unchanged, one negative-results entry.
