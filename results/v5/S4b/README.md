# S4b heat of combustion on measured rows only (v5)

Status: SCORED 2026-11-12 (readings declared before scoring; see Results). Declared 2026-11-12. Do not edit the readings after scoring; add a dated note. Self-declared in the same session as the code.

## Question

S4 found that heat of combustion from bulk descriptors beats a constant, the only T3 positive result. The Fable audit (docs/audit_fable_v5_2026-11-11.md, B3) notes that many AJFTD heat values are method D3338, an estimate computed from density, aromatics and distillation, i.e. from the model's own features; and that the constant was a mean over mixed methods (D4809 mean about 0.34 MJ/kg above D3338). Does the result survive on measured (calorimeter) values with a within-method constant?

## Data fixed in advance

* Same AJFTD table, features and splits as S4: `data/real_fuels/v5/ajftd_fuels_wide.csv` / `ajftd_fuels_long.csv`, `results/v5/S4/splits_ajftd_v5.json` (source groups, dev/locked roles). No new split is created.
* Method classes from the `method` field of the heat rows: MEASURED = D4809, D240 (bomb calorimeter); CALCULATED = D3338 (and D4529 if present); UNCLASSIFIED = codes whose meaning cannot be established from the repo (400C/D/E, 231A/B, others). The agent determines each code's class only from documentation inside the repo or the AJFTD page metadata already on disk, records the evidence per code in `method_classes.csv`, and leaves a code UNCLASSIFIED if none exists. No class is decided from the data's agreement with the model.
* Targets: net heat of combustion (MJ/kg). Features exactly as S4 (D86 points, aromatics, density) so the model comparison is S4's, only the rows change.

## Models (all leave-source-group-out on dev, exactly S4's procedure)

* C_all: mean over all-method dev training rows (S4's constant).
* C_within: mean over training rows of the SAME method class as the scored row.
* Ridge and HGB as in S4, trained on the subset named in each reading.

## Readings (fixed before scoring)

R-S4b-1 (measured only): train and test on MEASURED rows only. A model has "skill" if the paired bootstrap (over source groups, 5000, seed 0) 95% interval of MAE(C_within) minus MAE(model) excludes 0 on the positive side, with at least 8 source groups. Report n fuels and groups.
R-S4b-2 (calculated only): same on CALCULATED rows only. Expectation (not a reading): high skill because the label is a function of the features. If skill appears here and not in R-S4b-1, S4's heat result is attributed to the label construction.
R-S4b-3 (mixed, S4 as published): re-report S4's dev result against C_within instead of C_all, to show how much of the original gain was a method offset.
R-S4b-4 (locked): the locked rows were already scored once in S4. Report MEASURED-only locked rows (n and groups) as a second read, labelled "not a clean test". No new claim rests on it.
R-S4b-5: verdict sentence: "heat skill survives on measured values" only if R-S4b-1 shows skill for ridge or HGB; "not supported" if no model has skill with at least 8 groups; "not evaluable" if fewer than 8 measured groups exist.
R-S4b-6: nothing is written to `configs/`. The S4 flash correlations (coefficients quoted from memory, unverified) are NOT re-run or defended here; this note records that they remain unverified.

A reading that cannot be evaluated is stated "not evaluable" with the reason.

## Results

Scored 2026-11-12 by `scripts/v5/s4b_heat_measured.py` (tests: `tests/test_s4b_heat.py`). Files: `s4b_reproduction_check.csv`, `method_classes.csv`, `s4b_predictions.csv`, `s4b_metrics.csv`, `s4b_bootstrap.csv`, `s4b_readings.csv`. Statistic: pooled MAE(C_within) - MAE(model) with a paired bootstrap over source groups (cells), 5000 resamples, seed 0, 95% percentile; positive = model better. MAE in MJ/kg.

Step 0, reproduction: S4's published heat numbers (dev LOGO n 1683, 33 groups, const 0.214, hgb 0.155, ridge 0.160, with identical CIs; locked n 215, 8 groups) reproduce with C_all to < 5e-4 on every quantity (`s4b_reproduction_check.csv`, all ok).

Method classes (evidence per code in `method_classes.csv`): MEASURED = D4809 (also written "D 4809", "D4809 (SwRI[32])"; 364 rows, 228 fuels) and D240 (1 row); CALCULATED = D3338 ("D 3338" too; 773 rows, 719 fuels); D4529 has no heat rows in `ajftd_fuels_long.csv`. UNCLASSIFIED = 400C (1082 rows), 231B (267), 400D (46), 231A (36), 400E (8), 400IPA/B (3), empty method (27). The saved AJFTD pages carry only the code in the "Test Method" column, no legend, and the repo has no definition, so none was classified. The MEASURED and CALCULATED evidence is the repo's own audit note (docs/audit_fable_v5_2026-11-11.md, B3), not an independent standard text. Fuels with rows of two classes use only the rows of the class being scored (67 fuels have both D3338 and D4809 rows); the fuel target is the median of per-raw-fuel class medians.

| Reading | Model | n fuels | groups | MAE model | MAE C_within | diff (C_within - model) | 95% CI | outcome |
|---|---|---|---|---|---|---|---|---|
| R-S4b-1 measured | ridge | 172 | 16 | 0.260 | 0.347 | 0.087 | 0.027 to 0.158 | skill |
| R-S4b-1 measured | hgb | 172 | 16 | 0.203 | 0.347 | 0.144 | 0.091 to 0.213 | skill |
| R-S4b-2 calculated | ridge | 605 | 19 | 0.098 | 0.158 | 0.060 | -0.117 to 0.199 | no skill shown |
| R-S4b-2 calculated | hgb | 605 | 19 | 0.043 | 0.158 | 0.115 | 0.087 to 0.226 | skill |
| R-S4b-3 mixed (S4) | C_all | 1683 | 33 | 0.214 | 0.191 | -0.023 | -0.046 to -0.002 | C_within is the better constant |
| R-S4b-3 mixed (S4) | ridge | 1683 | 33 | 0.160 | 0.191 | 0.031 | 0.008 to 0.058 | beats C_within |
| R-S4b-3 mixed (S4) | hgb | 1683 | 33 | 0.155 | 0.191 | 0.036 | 0.004 to 0.073 | beats C_within |
| R-S4b-4 locked measured | C_all | 47 | 5 | 0.582 | 0.414 | -0.169 | -0.240 to -0.058 | not a clean test |
| R-S4b-4 locked measured | ridge | 47 | 5 | 0.276 | 0.414 | 0.138 | 0.086 to 0.217 | not a clean test |
| R-S4b-4 locked measured | hgb | 47 | 5 | 0.252 | 0.414 | 0.161 | 0.101 to 0.246 | not a clean test |

* R-S4b-1: evaluable (16 groups >= 8). Both ridge and hgb show skill over the within-method constant on measured values. The S4 constant (mean over all methods) was a worse constant on these rows (MAE 0.407 vs 0.347).
* R-S4b-2: hgb shows skill, ridge does not (CI includes 0). Skill appears on measured values too, so the "attributed to label construction" condition in R-S4b-2 is not triggered.
* R-S4b-3: in S4's mixed set the within-method constant is better than C_all (0.191 vs 0.214), so about 0.023 of S4's 0.059 (hgb) / 0.054 (ridge) MAE gain over C_all was a method offset; the gain over C_within is 0.036 (hgb) and 0.031 (ridge), both CIs above 0. Here the C_within class is the set of row classes of a fuel (e.g. "CALCULATED+MEASURED", "UNCLASSIFIED"); UNCLASSIFIED codes are pooled as one class, which they may not be.
* R-S4b-4: 47 locked fuels with measured heat rows, 5 locked cells; reported as a second read, "not a clean test"; no claim rests on it, and no group threshold or skill label is applied.
* R-S4b-5 verdict: heat skill survives on measured values.
* R-S4b-6: nothing was written to `configs/`. The S4 flash correlations (coefficients quoted from memory) remain unverified and were not re-run or defended.

Caveats (added after scoring, not part of the declared readings): the measured subset is small (172 fuels, 16 cells) and its features are the same bulk descriptors; a held-out cell is a source-year cell, and leakage across cells of one source_group is possible (see S4 README section 2). If a cell has no training rows of its class the training mean is used (R-S4b-3 only; count in `s4b_metrics.csv` note); R-S4b-1/2 train on the same subset so no fallback arises.
