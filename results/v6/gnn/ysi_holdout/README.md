# Track B YSI locked holdout (declared 2026-11-21, before any fit on the split)

Question: does the multi-task D-MPNN keep its YSI advantage over the v4 ridge group-contribution model and the constant on compounds it has
never seen? In Track B CV the GNN had MAE 69.4 vs ridge 106.5 (n 224); that CV had no locked holdout.

## Lock (`lock.json`, `scripts/v6/gnn/ysi_holdout.py --stage lock`)
Unit: generic Murcko scaffold group (acyclics bucketed by carbon count), same grouping as Track B. Outcome-blind rule: groups in
md5(`v6ysi:`+group) order are locked while locked YSI compounds stay <= 25 % of the 224. Result: 56 compounds in 17 of 36 groups. The locked
compounds are removed from training for ALL ten targets (not just YSI). `pool.csv` hash is stored in the lock.

## Candidate and comparators
- GNN: v1 Track B recipe unchanged (5-model ensemble, up to 150 epochs, early stopping on a group-held-out 15 % validation split of the
  training groups, two configurations chosen on that split). GPU/CPU float differences are possible; v2 settings are not used.
- Constant: median YSI of the training compounds. Ridge: v4 `ysi_model` GC ridge fitted on the same training compounds.

## Gate
90 % paired scaffold-group bootstrap (2000 draws) of MAE gain over the constant AND over the ridge both have lower bound > 0. Scored once
(`scored.flag`). One-sided over-prediction mean is reported; no multiplicity correction (single head).
Margin: 90th percentile of absolute out-of-fold YSI error from the existing Track B CV (`oof.csv`) on the non-holdout compounds, applied to
the holdout and its coverage reported. Disclosure: those CV folds trained on holdout compounds, so the margin is a rough one, not a clean
conformal calibration.

## Expected outcome
A positive gain over the constant is near certain. Versus the ridge a tie or win is plausible but 17 locked groups give a wide CI, and the
GNN's CV gain could shrink on unseen scaffolds. A failure is logged in `docs/v6_negative_results.md`.

## Result (scored once, 2026-11-21 on GPU; `readings.json`, `holdout_pred.csv`)
n 56 compounds, 17 scaffold groups. MAE: GNN 59.6, ridge 133.6, constant 366.0 (YSI units). Gain over ridge 90 % CI [32.6, 154.0]; over
constant [223.9, 435.2]. **Gate passes.** Chosen config hidden 300 / depth 3 / dropout 0.2. Mean one-sided over-prediction 49.2. Rough 90 %
absolute margin 190.1 (from non-holdout CV errors) covers 91 % of the holdout. Caveats: one head, 17 groups, no lab data beyond the pool; the
constant is a weak bar (the YSI values span orders of magnitude); the margin is not a clean conformal calibration (see above). Model-level
result only; YSI is not exposed in `predict_fuel`.
