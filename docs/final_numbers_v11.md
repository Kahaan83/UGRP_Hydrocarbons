# v11 final numbers: tuning the viscosity pair model (negative result)

Source: `results/v11/tune_viscosity/summary.json` and `search.csv`; protocol `docs/v11_predeclaration.md`; script `scripts/v11/tune_pair_viscosity.py`. Dev systems only; locked systems not read.

| item | value |
|---|---|
| dev systems | 264 |
| configurations searched | 30 random plus the v5 default |
| v5 default, per-system macro MAE (ln), v5 folds | 0.05989 |
| best tuned (800 trees, lr 0.02, 15 leaves, min child 20, colsample 0.8, lambda 20), same folds | 0.05847 |
| confirmation on 3 reshuffled partitions: mean gain (default minus tuned) | 0.0010 |
| 95 % system-bootstrap interval of the gain | -0.0002 to 0.0023 |
| adoption rule (interval above 0 and gain at least 0.002) | not met |

Reading: tuning does not give a gain that survives reshuffled partitions. The v5 settings stay, `results/v5/BASE/` and its manifest are unchanged, and the density pair model was not tuned (predeclared: only if viscosity is adopted). Status: no improvement shown.
