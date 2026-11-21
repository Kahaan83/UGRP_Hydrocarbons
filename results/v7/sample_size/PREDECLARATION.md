# Task 8 predeclaration: sample-size learning curves

Design: for each outer fold (or LOGO held-out entity), only the training groups are subsampled at fractions 0.1 0.2 0.4 0.6 0.8 1.0, 3 repeats (nested across fractions for a fixed repeat). The test fold is untouched. Group = scaffold group (Track B) or entity (freezing). MAE is scored on the full outer test fold, per target.
Models in this run (dev rows; ysi_holdout masked; no holdouts touched):
- `dmpnn_trackb` (GPU): v1 config 0 (hidden 200, depth 4, dropout 0.1), 2 seeds, epochs 300, patience 50, early stopping on 15 % of the subsampled training groups.
- `gbm_trackb` (CPU): RDKit-descriptor LightGBM with fixed parameters (num_leaves 15, lr 0.05, 300 trees, min_child_samples 5, subsample and colsample 0.8), no inner tuning.
- `freezing_ridge` and `freezing_lgbm` (CPU): `fp_ridge` and `make_lgbm("pred_e")` from `scripts/v6/freezing/step3_models.py`, LOGO over the 52 dev entities.
Not in this run: the v5 viscosity and density pair residuals. Their training rows come from long regeneration pipelines (`scripts/v5/base/delta_ml_*.py`) that Task 13 has to refactor for its reproduction gate; their sample-size curves are run there. This is a stated deviation from the plan text.
Reading rule: slope of MAE vs log(n_train rows) over the top two fractions (0.8, 1.0), 90 % bootstrap CI over (repeat x fold) units. A CI that excludes 0 with a negative slope: more data is expected to help. A CI that includes 0: the curve has flattened. No forecast beyond the data. Reported per target; the freezing CIs rest on few units and are read as wide.
