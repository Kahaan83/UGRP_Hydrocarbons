# v11c: minor optimisations (record of rules)

Rules for items 3 and 5 were written in the script headers (`scripts/v11/seed_avg.py`, `pure_ensemble.py`, `tune_rf_pure.py`) before each was run; this file collects them. Item 1 is a data correction, not a tuned quantity.

| item | question | rule |
|---|---|---|
| 1 | o-/m-xylene eutectic composition transposed | correction applied if ideal Schroeder-van Laar at the published eutectic temperature (paper's own Tm and dHfus) reproduces the corrected order and not the entered one; checked by consistency, not against the page image |
| 2 | add v9 density systems to the density dev pool | applicable only to v9 systems not in v5 and not in `results/v9/lock.json` |
| 3 | average 5 LightGBM seeds in the pair models | paired system-bootstrap 95 % interval of the gain above 0 and mean gain at least 1 % of the single-seed MAE |
| 4 | viscosity: tuned model in warm bins only (hybrid) | post hoc (suggested after seeing the v11b margin-guard failure) and adopted at the user's request. Cold bins (T below 283.15 K) keep the v5 model unchanged; BASE `delta_ml_viscosity` regenerated; design effect checked by re-running pool P0 with and without the warm model |
| 5 | ensemble the deployed pure-compound models; tune the random forest | per property: gain interval above 0 and Holm-adjusted one-sided p below 0.05 over 10 properties |
| 6 | make the fresh v7 flash lock scorable | needs verified Antoine (40 to 100 C window) and LFL for the 34 missing compounds |
