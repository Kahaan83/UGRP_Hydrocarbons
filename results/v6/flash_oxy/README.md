# v6 flash arm: oxygenate mixtures, result (scored 2026-11-21; declaration and Amendment 1 in `PREDECLARATION.md`, committed before any fit)
Files: `rows.csv`, `lock.json`, `dev_oof.csv`, `locked_pred.csv` (scored once), `readings.json`, `chosen_dev.json`. Model-level only.

## Primary claim (oxy, arm B MolPool vs gate comparator = descriptor GBM, chosen on dev OOF)
| | MAE comparator (desc GBM) | MAE MolPool | gain 90 % CI (systems bootstrap) | n_comp >= 3 MAE (comp / MolPool) |
|---|---|---|---|---|
| dev OOF (57 systems) | 4.08 | 3.44 | [0.16, 1.04] | 3.40 / 2.25 |
| locked (21 systems, once) | 4.48 | 3.88 | [-0.47, 1.38] | 4.82 / 3.52 |
Other comparators (dev / locked MAE): linear molar 7.20 / 13.19, blending index 4.15 / 4.54. Le Chatelier is not evaluable on oxygenates
(0 % of dev oxy rows have Antoine and LFL constants). **Primary gate fails**: the dev interval clears zero, the locked interval does not.
Point estimates favour MolPool on both sets, so the result is "not established", not "worse".
Safety: p90 over-prediction MolPool 5.16 / 8.72 vs GBM 6.20 / 9.26 (dev / locked); mean signed error on locked is +1.7 for both, so
both over-predict on average on the locked oxygenates. No veto.

## Secondary claim (hydrocarbon non-inferiority, arm B minus arm A MAE)
Locked (14 systems, 31 rows): A 8.99, B 7.88, 90 % upper bound 0.65 <= 1.0, so the declared criterion is met. Dev OOF (17 systems): A 3.94,
B 5.50, upper bound 3.29, which fails. The locked pass rests on 14 systems of one-row-sized groups and on arm A being unusually bad there
(mean signed +8.6, i.e. strong over-prediction, comparable to linear molar 8.99); the dev result contradicts it. Treat non-inferiority as
**not established**; the oxygenate model is not to be used for hydrocarbon fuels.

## Reading
Adding oxygenate data plausibly helps oxygenate blends over a descriptor GBM (consistent sign on dev and locked, larger gain at
n_comp >= 3), but the locked set (21 systems) cannot confirm it. Reporting-only columns (benchmark splits, random 10-fold,
leave-component-out) were not run. Status: no skill shown against the strongest comparator.
