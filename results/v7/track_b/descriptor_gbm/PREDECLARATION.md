# Task 6 predeclaration: descriptor GBM and RF comparators for Track B

Rows: v6 Track B pool, same outer folds as Task 4 (full-pool generic-Murcko GroupKFold(5)), ysi_holdout labels masked. GNN comparator = v6 `results/v6/gnn/oof.csv` restricted to the same reduced rows.
Features: RDKit 2D descriptors (all `Descriptors.descList`), constant columns dropped and NaN median-imputed inside each fold (RF); LightGBM takes NaN natively.
Models per target (single task): LightGBM (inner Optuna, 30 trials, inner val = 15 % of train groups; num_leaves, min_child_samples, learning_rate, n_estimators with early stopping) and RandomForest (500 trees, defaults).
Metric: OOF MAE per target. Gain = MAE(better of GBM, RF) - MAE(GNN), paired group bootstrap (4000 draws, group = scaffold).
Family: 10 targets, Holm correction on bootstrap two-sided p-values at level 0.10.
Reading rule: "GNN earns its place on target t" only if the Holm-adjusted gain interval excludes 0 with gain > 0. Otherwise report "no evidence the GNN beats a descriptor tree model" for t.
