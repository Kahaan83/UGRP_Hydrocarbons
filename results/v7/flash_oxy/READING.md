# Task 11 reading: flash-point oxygenate retrain (arm B all dev rows, arm A hydrocarbon rows)
Rows are the v6 flash_oxy dev rows (already seen), exploratory; the flash_oxy_lock was not scored. Protocol: `PREDECLARATION.md`; outputs `retrain/cv_summary.csv`, `retrain/oof.csv`. Gain = comparator MAE minus v7 MAE (degC), 90 % group bootstrap, Holm over two comparisons.

- Oxygenate rows: 641 rows, 57 systems. v7 arm B MAE 3.608 degC; v6 MolPool v1 3.439; gate comparator (desc_gbm) 4.077.
- Gain vs comparator: +0.469 [0.027, 0.908], Holm p 0.162, so not significant after correction (`beats_comparator` False). Gain vs v1: -0.169 [-0.645, 0.314], Holm p 0.638 (no evidence of change; v1 not beaten).
- Multi-component rows (209): v7 3.221 vs comparator 3.403.
- Safety (signed over-prediction, the dangerous direction): mean over-prediction v7 1.730 degC, comparator 2.329, v1 1.770; p90 signed error v7 5.121, comparator 6.203, v1 5.156. `unsafe_tail` False (no worse tail than the comparator).
- Hydrocarbon non-inferiority (95 rows): arm A MAE 4.210, arm B 6.398, B minus A 2.188 [-0.692, 4.859], upper bound above the 1.0 degC margin, so `hc_noninferior` False. Mixing oxygenates into the training pool still costs hydrocarbon accuracy at this sample size.
- Status: the gate (beat the comparator with Holm significance and hydrocarbon non-inferiority) is not met. No skill shown over the comparator; v6 flash_oxy verdict unchanged. Nothing adopted into `predict_fuel()`.
