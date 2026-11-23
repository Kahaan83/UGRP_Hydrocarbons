# Revision computations (post v11c)

No BASE file was changed. Scripts are in `scratch/` (rm03_locked.py, rm03_04.py, rm05_ternary.py, rm07.py); outputs in `results/v11/claim_risk/` (rm03_locked_bootstrap.json, rm04_rowcount.json, rm05_ternary.json, rm07_step.json). Models are refit on dev systems only; locked sets are scored with the final dev fit. Bootstrap: 5000 system-level resamples, seed 42, percentile 95 % intervals.

## RM-03 locked paired bootstrap (deployed tuned models)
Viscosity, 70 locked systems, T >= 283.15 K (macro MAE, ln):
- Physics base 0.0544; untuned pair 0.0311; deployed hybrid 0.0290 (reproduces BASE).
- Tuned vs physics base: gain 0.0254, CI [0.0182, 0.0326].
- Tuned vs untuned pair: gain 0.0021, CI [-0.0002, 0.0049] (includes zero).
- Untuned vs physics base: gain 0.0233, CI [0.0163, 0.0302].
- Locked coverage under dev margin 0.1195: 63 of 70; under nested margin 0.1258: also 63 of 70.

Density, 41 locked systems (kg/m3):
- Ideal 1.2008; untuned pair 0.8972; deployed tuned 0.8389 (reproduces BASE).
- Tuned vs ideal: gain 0.3619, CI [0.1034, 0.6098].
- Tuned vs untuned: gain 0.0583, CI [0.0008, 0.1222] (barely excludes zero).
- Untuned vs ideal: gain 0.3036, CI [0.0241, 0.5706].

Q1, the 3 "untraceable" locked viscosity systems: ids 43649dc66e54 (ternary, 45 rows), 846723718924 (ternary, 45 rows), f1d5f2f85452 (binary, 15 rows), all source `boehm2022`, 233.15-288.15 K. They come from `data/real_fuels/v5/boehm2022_viscosity.csv`, which `delta_ml_viscosity.py` concatenates with `measurements_dedup.csv`; they are not in the dedup table by design. All 70 systems, including these 3, are in the locked MAE and in the 63/70 count.

## RM-04 row-count confounding (dev OOF, deployed hybrid, rows_scored.csv)
Spearman between rows per system and the system worst residual:
- T >= 283.15 K: n=264, mean 45.3 rows, rho 0.288, p 2e-6. Margin 0.1195; systems with <=10 rows (n=54): 0.1132; >10 rows (n=210): 0.1258; per-system cap of 5 rows: median 0.1004; cap 10 rows: 0.1074; mean-residual score: 0.0626.
- 263-283 K: n=16, mean 5.3 rows, rho 0.294, p 0.27. Margin 0.6084.
- T < 263 K: n=20, mean 5.65 rows, rho 0.594, p 0.006. Margin 0.2875; cap 3 rows: median 0.1532; cap 5: 0.2635.
Reading: the worst-residual score grows with row count, so the margin is set mainly by data-rich systems. It is conservative for sparse systems; it is not optimistic for them.

## RM-05 ternary density (19 ternary dev systems, 615 rows, one system with 216 rows), out-of-fold with deployed params
| model | macro MAE (19 sys) | row MAE | macro excl. 216-row system (18) | row MAE excl. |
|---|---|---|---|---|
| ideal | 2.2074 | 3.6439 | 1.9589 | 2.0002 |
| v5 pair | 2.2044 | 3.9020 | 1.9341 | 2.1862 |
| tuned pair | 2.0544 | 3.6193 | 1.7919 | 1.9089 |
The tuned pair model beats the ideal rule on ternaries under row-level and macro scoring, with and without the 216-row system. The 2.207 to 2.054 macro figures are ideal vs tuned, not v5 vs tuned. All dev: macro 1.3335 / 0.9578 / 0.9076 (ideal / v5 / tuned).

## RM-07 283.15 K step (dev-fit models, rows 283.15-293.15 K, 1601 rows, 169 systems)
|warm correction - cold correction|: mean 0.0119 ln (1.27 %), median 0.0047, p95 0.0456 (4.7 %), p99 0.0779, max 0.627 (87 %). Mean absolute cold correction 0.0446, warm 0.0369. The mean step is about 5.6 times the 0.0021 ln nested gain.

## RM-08 conflict sensitivity (Q4)
- The lab viscosity selection (`lab_selection`, `~conflict`) and lab density selection (`eval_lab_known_composition.py` line 149) already exclude flagged conflicts. Viscosity and density margins are therefore unaffected by the 1,221 flags.
- Flag breakdown by kind: 806 real-fuel (310 freezing, 46 flash), 385 pure (241 flash, 92 vapor pressure, 51 viscosity), 30 binary. 739 of the 1,221 have two observations.
- Real-fuel freezing: 15 of the 50 headline fuels map to AJFTD ids; 5 of those are flagged. Excluding them (45 fuels): MAE 7.30, bias +7.30, constant baseline MAE 5.81, 95.6 % under the 13.843 C margin. Still no skill. The mapping covers only the 15 AJFTD-id fuels.
- Not done: bicyclohexyl dHfus and the o-/m-xylene source check need primary literature. The xylene composition was corrected in commit 341a056.

## RM-13 flash class diagnostics (35 mesitylene-excluded systems)
- The 10.60 C case is system 39a914ff0579 (1,4-dimethylcyclohexane + n-dodecane), a single measured point (Leenhouts). It is not bicyclohexyl or cis-decalin, so the registry overlay is not implicated.
- The 2.207 C margin is the third largest of 35 per-system worst over-predictions (10.60, 3.79, 2.21). Margin guarantee for cyclic blends rests on one point.
- n-Dodecane Antoine: refit_measured, covers the window, so the AGENTS trap 3 is not the cause. Cause not established.
- The "61 unscored systems (658 rows)" in the review: dedup holds only 32 binary/multi flash systems, all scored, so this figure was not reproduced.
