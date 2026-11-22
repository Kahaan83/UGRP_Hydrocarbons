# Post-hoc addendum to Task 5: long-cap D-MPNN vs Task 6 trees (labelled post-hoc; chosen after seeing the Task 5 and Task 6 readings)

The first Task 5 run did not save OOF predictions, so the long-cap regime (epochs 1000, patience 100, both v1 configs, 5 seeds, config chosen per fold on inner val) is rerun with `--long-only`, which saves `oof.csv`. Same rows, folds and inner splits as Task 5. GPU nondeterminism means the MAE may differ slightly from the first run.
Comparison: `scripts/v7/track_b/long_vs_trees.py`. Gain = MAE(better of GBM, RF from Task 6) - MAE(long-cap GNN), paired group bootstrap (4000 draws), Holm over 10 targets at 0.10, same rule as Task 6. Because the long-cap setting was chosen after seeing Task 5, the result is exploratory and reported as such.
