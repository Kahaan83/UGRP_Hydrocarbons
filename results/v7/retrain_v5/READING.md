# Task 13 reading: CPU tree models, reproduction gate, curves and nested Optuna
Protocol: `PREDECLARATION.md`. All six reproduction gates passed (`<model>/repro_gate.json`; exact for five, freezing_step3 within the 4-decimal precision of the v6 files, fusion_estimate with rows pinned to the BASE compounds, ysi_lgbm gate on all 222 rows and comparison on the 169 dev compounds). Dev rows only; nothing locked was scored; `src/` and `results/v5/BASE/` hashes unchanged in every run. Gain = refit-original MAE minus v7 MAE, paired group bootstrap 90 % CI; rule: "v7 improves" only if the interval excludes 0 with a positive gain. All of this is exploratory tuning on already-seen dev rows; nothing is adopted into `predict_fuel()` (a separate user decision). No multiplicity correction across models.

| model | metric | v7 | refit original | gain [90 % CI] | reading |
|---|---|---|---|---|---|
| viscosity_pair | macro MAE ln visc, 264 systems | 0.0604 | 0.0599 | -0.0005 [-0.0021, 0.0012] | no evidence of improvement |
| density_pair | macro MAE kg/m3, 167 systems | 0.9872 | 0.9578 | -0.0295 [-0.0839, 0.0238] | no evidence of improvement |
| ysi_lgbm | OOF MAE YSI, 169 compounds, 21 groups | 220.1 | 236.8 | +35.9 [11.9, 64.2] | improves on the v5 tree; constant baseline 234.1; the shipped model is the ridge, which this does not change |
| fusion_estimate | LOFO Tm MAE (K), 109 compounds, 9 families | 50.48 | 54.69 | -4.9 [-16.4, 5.2] (family-macro) | no evidence of improvement; pooled MAE is 4.2 K lower but the family-macro gain is negative, so the pooled difference is driven by a few families. C8 to C20: 43.41 vs 47.75, gain -7.3 [-17.0, 1.8]. Joback baseline 60.1 K |
| freezing_step3 | dev LOGO MAE degC, 52 groups | 3.068 | 3.315 | +0.247 [-0.111, 0.607] | no evidence of improvement; without fuel 64: 2.766 vs 3.015; constant 6.635; over-prediction mean 1.27 vs 1.58 |
| mixture_gbr | dev OOF MAE, 5 units | see below | | | 0 of 5 units with interval above 0 |

mixture_gbr per unit (v7 vs refit original MAE, gain [CI]): boehm density 3.686 vs 3.885, +0.178 [-0.399, 0.746]; boehm ln_visc 0.0625 vs 0.0551, -0.0053 [-0.0103, -0.0003] (v7 significantly worse); kaust mon 2.849 vs 3.154, -0.159 [-0.682, 0.391]; kaust ron 2.023 vs 1.486, -0.728 [-1.405, 0.043]; leenhouts mon 2.779 vs 3.084, +0.219 [-0.160, 0.657].

Curves (final fits, median over outer folds; LightGBM curves hold the validation l2 only, so no train/val gap is read):
- Best iteration is small: viscosity_pair 28 (stop 127), density_pair 67 (167), ysi_lgbm 32 (131; one fold hit the 5000 cap), fusion_estimate 97 (197), freezing_step3 6 (105). mixture_gbr: boehm density 1, boehm ln_visc 20, kaust ron 45, kaust mon 159, leenhouts mon 464. So the models stop within a few hundred rounds at low learning rates and the 5000-round cap is not the binding limit; with best iterations of 1 to 30 on several models the validation curve is flat or noisy, consistent with little signal to learn beyond the physics baseline.
- Stated limits: freezing_step3 (52 rows) and the Boehm/KAUST units have few groups; the gain intervals are wide and one-sided claims are not supported.

Status: tuning gives no reliable improvement over the v5/v6 tree models except ysi_lgbm, whose comparator is not the shipped model. No change to any v5 or v6 verdict.
