# v4 Task 11 - viscosity: nested CV and physics-structured residuals

Files: `viscosity_leaderboard_v4.csv` (all rows, CIs, paired deltas), `viscosity_predictions_v4.csv`, `selected_hyperparameters_v4.csv`.
Code: `scripts/v4/build_viscosity_leaderboard_v4.py`, `src/fuelblend/mixture/structured_residual.py`, tests `tests/test_v4_viscosity.py`.
Tables below were recomputed from the CSVs (generator not committed).

## Protocol
- Data: 57 CRC 647 jet fuels, kinematic viscosity at +20, -20, -40 C (ASTM D455, `raw_file` provenance from `jet_properties_v4.csv`) = 171 obs.
- **Dev = 45 fuels, locked = 12 fuels** (`splits/v4/real_fuels_v4_splits.json`). Locked fuels are never used for fitting, selection or tuning. Dev rows are leave-one-fuel-out (LOO); hyperparameters (bilinear rank/lambda, ridge alpha, GP length-scale/noise) are chosen by inner grouped 5-fold CV inside each training fold. Locked rows are predicted by models trained on all 45 dev fuels (hyperparameters from inner CV on dev only). Zero-parameter rows (raw Refutas, raw Walther) and the lab zero-shot row involve no CRC fitting.
- Physics baseline: mass-fraction Refutas from pure-component eta(T), rho(T); residual models learn delta ln(nu) vs raw Refutas. Walther (D341) is anchored on Refutas at +20/-20 C and extrapolated to -40 C (v3 definition; the -40 C row is the only one that differs from its anchor model).
- Structured residual: delta ln(nu) = a0(z) + a1(z)*1000/T, z = 5 class mass fractions + mean/std carbon number; exactly linear in 1000/T (second difference test). Ridge: unpenalised intercept + 1000/T slope, penalised z-dependence. GP: OLS recalibration then zero-mean GP with kernel RBF(z)*(1+t t'); 90% interval = +-1.645 predictive sd (stage-1 coefficient uncertainty ignored).
- Lab-trained row: ridge structured residual fit on lab binary/multi hydrocarbon mixtures (dedup, non-conflict, 253-340 K; residual vs the same Refutas), applied zero-shot to CRC fuels.
- Bootstrap: cluster bootstrap by fuel (existing `cluster_bootstrap_ci`, 2000 resamples) and paired fuel-cluster bootstrap of each model minus recalibrated Refutas. A model is only called better/worse where the paired CI excludes 0.

## Headline (dev, LOO by fuel)
**dev (LOO by fuel), overall (all T)** (N = 45 fuels, 135 obs)

| Model | log-RMSE [95% CI] | log-MAE | MAPE % | bias | paired dRMSE vs recal [95% CI] | CI excl. 0 |
|---|---|---|---|---|---|---|
| Raw Refutas | 0.4267 [0.3959, 0.4549] | 0.3281 | 25.53 | -0.3011 | +0.2916 [0.2562, 0.3195] | True |
| Recalibrated Refutas | 0.1351 [0.1157, 0.1550] | 0.1042 | 10.45 | -0.0000 | - | - |
| Walther (D341) | 0.3872 [0.3567, 0.4152] | 0.3020 | 24.09 | -0.2750 | +0.2521 [0.2173, 0.2794] | True |
| Walther + recalibration | 0.1633 [0.1384, 0.1898] | 0.1221 | 13.42 | +0.0701 | +0.0281 [0.0066, 0.0486] | True |
| Bilinear pairwise | 0.0901 [0.0799, 0.1014] | 0.0706 | 7.17 | +0.0160 | -0.0451 [-0.0630, -0.0267] | True |
| Structured ridge | 0.0826 [0.0731, 0.0949] | 0.0637 | 6.38 | -0.0019 | -0.0525 [-0.0734, -0.0323] | True |
| Structured GP | 0.0888 [0.0772, 0.1012] | 0.0682 | 6.86 | +0.0007 | -0.0464 [-0.0621, -0.0302] | True |
| Lab-trained structured residual (zero-shot) | 0.3952 [0.3652, 0.4226] | 0.3021 | 23.87 | -0.2786 | +0.2601 [0.2256, 0.2872] | True |

**dev (LOO by fuel), T = -40 C** (N = 45 fuels, 45 obs)

| Model | log-RMSE [95% CI] | log-MAE | MAPE % | bias | paired dRMSE vs recal [95% CI] | CI excl. 0 |
|---|---|---|---|---|---|---|
| Raw Refutas | 0.6735 [0.6247, 0.7178] | 0.6519 | 47.14 | -0.6519 | +0.4947 [0.4499, 0.5310] | True |
| Recalibrated Refutas | 0.1788 [0.1454, 0.2113] | 0.1458 | 14.04 | -0.0461 | - | - |
| Walther (D341) | 0.5977 [0.5494, 0.6415] | 0.5736 | 42.85 | -0.5736 | +0.4188 [0.3751, 0.4540] | True |
| Walther + recalibration | 0.2391 [0.1995, 0.2795] | 0.1995 | 22.95 | +0.1644 | +0.0602 [0.0139, 0.1062] | True |
| Bilinear pairwise | 0.0775 [0.0617, 0.0915] | 0.0589 | 5.78 | -0.0165 | -0.1014 [-0.1253, -0.0760] | True |
| Structured ridge | 0.0801 [0.0630, 0.0954] | 0.0634 | 6.08 | -0.0498 | -0.0988 [-0.1308, -0.0671] | True |
| Structured GP | 0.0918 [0.0712, 0.1119] | 0.0703 | 6.77 | -0.0447 | -0.0871 [-0.1141, -0.0598] | True |
| Lab-trained structured residual (zero-shot) | 0.6254 [0.5775, 0.6688] | 0.6031 | 44.54 | -0.6031 | +0.4466 [0.4029, 0.4824] | True |

## Locked fuels (untouched; 12 fuels, wide CIs)
**locked (untouched), overall (all T)** (N = 12 fuels, 36 obs)

| Model | log-RMSE [95% CI] | log-MAE | MAPE % | bias | paired dRMSE vs recal [95% CI] | CI excl. 0 |
|---|---|---|---|---|---|---|
| Raw Refutas | 0.3806 [0.3135, 0.4564] | 0.2915 | 23.49 | -0.2428 | +0.1866 [0.1070, 0.2693] | True |
| Recalibrated Refutas | 0.1939 [0.1341, 0.2594] | 0.1491 | 15.19 | +0.0584 | - | - |
| Walther (D341) | 0.3427 [0.2740, 0.4189] | 0.2656 | 21.94 | -0.2169 | +0.1488 [0.0706, 0.2303] | True |
| Walther + recalibration | 0.2408 [0.1845, 0.3029] | 0.1954 | 21.43 | +0.1284 | +0.0469 [0.0129, 0.0768] | True |
| Bilinear pairwise | 0.1670 [0.0987, 0.2445] | 0.1060 | 10.51 | +0.0397 | -0.0269 [-0.0612, -0.0061] | True |
| Structured ridge | 0.1497 [0.0810, 0.2289] | 0.0905 | 8.70 | +0.0113 | -0.0442 [-0.0911, -0.0174] | True |
| Structured GP | 0.1517 [0.0860, 0.2296] | 0.0920 | 8.90 | +0.0181 | -0.0423 [-0.0871, -0.0163] | True |
| Lab-trained structured residual (zero-shot) | 0.3523 [0.2857, 0.4268] | 0.2653 | 21.63 | -0.2214 | +0.1584 [0.0828, 0.2381] | True |

**locked (untouched), T = -40 C** (N = 12 fuels, 12 obs)

| Model | log-RMSE [95% CI] | log-MAE | MAPE % | bias | paired dRMSE vs recal [95% CI] | CI excl. 0 |
|---|---|---|---|---|---|---|
| Raw Refutas | 0.5903 [0.4881, 0.7099] | 0.5654 | 42.44 | -0.5654 | +0.4159 [0.3424, 0.4873] | True |
| Recalibrated Refutas | 0.1744 [0.1105, 0.2440] | 0.1459 | 14.56 | +0.0404 | - | - |
| Walther (D341) | 0.5161 [0.4120, 0.6351] | 0.4878 | 37.80 | -0.4878 | +0.3417 [0.2665, 0.4120] | True |
| Walther + recalibration | 0.3027 [0.2507, 0.3444] | 0.2849 | 33.28 | +0.2505 | +0.1283 [0.0320, 0.2016] | True |
| Bilinear pairwise | 0.1133 [0.0712, 0.1542] | 0.0911 | 9.17 | +0.0316 | -0.0611 [-0.0925, -0.0295] | True |
| Structured ridge | 0.0711 [0.0411, 0.0960] | 0.0565 | 5.50 | -0.0268 | -0.1033 [-0.1520, -0.0531] | True |
| Structured GP | 0.0700 [0.0359, 0.0969] | 0.0507 | 4.94 | -0.0174 | -0.1044 [-0.1502, -0.0577] | True |
| Lab-trained structured residual (zero-shot) | 0.5442 [0.4426, 0.6619] | 0.5181 | 39.68 | -0.5181 | +0.3698 [0.2967, 0.4393] | True |

Other temperature slices (+20 C, -20 C; N = 45 dev / 12 locked obs each) are in the CSV.

## Structured GP: 90% interval coverage on held-out fuels
| split | condition | N obs | coverage | mean width (log) |
|---|---|---|---|---|
| dev (LOO by fuel) | overall (all T) | 135 | 0.770 | 0.212 |
| dev (LOO by fuel) | T = 20 C | 45 | 0.956 | 0.217 |
| dev (LOO by fuel) | T = -20 C | 45 | 0.600 | 0.204 |
| dev (LOO by fuel) | T = -40 C | 45 | 0.756 | 0.214 |
| locked (untouched) | overall (all T) | 36 | 0.694 | 0.247 |
| locked (untouched) | T = 20 C | 12 | 0.917 | 0.250 |
| locked (untouched) | T = -20 C | 12 | 0.333 | 0.243 |
| locked (untouched) | T = -40 C | 12 | 0.833 | 0.248 |

Nominal 0.90: the GP is under-covering (overall 0.77 dev, 0.69 locked), especially at -20 C; its intervals should not be used as calibrated.

## Reading the results (only where paired CI supports it)
- Raw Refutas and raw Walther are strongly biased low (-40 C bias about -0.65 and -0.57 in ln); recalibration removes most of it. Walther + recalibration is worse than recalibrated Refutas (CI excludes 0, dev and locked overall).
- Bilinear, structured ridge and structured GP all beat recalibrated Refutas on dev overall and at -40 C (paired CI excludes 0). On locked fuels the paired CI also excludes 0 for these three at overall and -40 C, but with only 12 fuels. Differences **among** bilinear / ridge / GP were not tested pairwise and their CIs overlap; no ranking among them is claimed.
- Lab-trained structured residual applied zero-shot does **not** fix the cryogenic collapse: it is essentially raw Refutas (-40 C log-RMSE about 0.63 dev, 0.54 locked vs recalibrated about 0.18, 0.17). In this ridge the intercept and pure 1000/T column are unpenalised, so the lab row applies a lab-fitted global offset and slope plus a near-zero z-correction (alpha = 1000, the grid maximum, shrinks only the z columns); dev bias moves from -0.301 (raw) to -0.279. So the "cannot blow up" property holds (no worse than raw physics) but it does not supply the missing offset. This is a negative result for lab-to-fuel transfer, not for the structure.
- Nested-CV selection varies by outer fold for ridge (alpha 0.01 x20, 0.1 x18, 1 x6, 10 x2), GP and bilinear (the latter two pick from a narrow part of the grid: bilinear always rank 2; GP always length-scale 8, noise 0.1 or 0.3). Counts are in `selected_hyperparameters_v4.csv`. Several choices sit at the grid edge: GP length-scale 8 = maximum (46 of 46 fits); GP noise ratio 0.1 = minimum (20 of 46); ridge alpha 0.01 = minimum (20 of 45); bilinear rank 2 = minimum (46 of 46) and bilinear weight-decay 1e-4 = minimum (30 of 46). The grids may be too narrow.

## Caveats
- **Compositions are reconstructed** (`provenance = reconstructed`, formula-derived from measured aromatics + D86 via carbon-number bins), not measured GC. Every error here mixes reconstruction error and physics error, and the structured/bilinear residuals can absorb reconstruction bias. The descriptors z are themselves functions of the same measured aromatics/D86, so the structured models are partly learning a fuel-level mapping from those inputs.
- Components are mapped to registry compounds with valid pure eta/rho by (class, carbon number) (as `fuel_hindcast.Reps`); perdeuterated entries use the parent class; the unresolved lump and unmappable mass (mean 0.24%, max 0.56%) are dropped and renormalised.
- CRC647 has only 3 temperatures per fuel, so per-T slices are 45 dev / 12 locked obs each; locked CIs are wide. Only 3 temperatures means a1(z) is barely identified beyond a line; the -40 C row is the cryogenic stress test.
- Bilinear predictions with `seed=0` are deterministic; its 46 fits chose rank 2 only.
- Numbers differ from v3 R1 (e.g. v3 recalibrated overall log-RMSE 0.1411 over all 57 with LOO training on every other fuel) because v4 trains dev folds on dev fuels only and uses v4 reconstructed compositions; not directly comparable.

## Claims-register rows this would update (not edited here)
Viscosity rows that cite the v3 R1 leaderboard (cryogenic collapse of lab-trained neural residuals; recalibrated Refutas / bilinear "winner" claims): the v4 result confirms recalibrated Refutas >> raw physics, weakens any "bilinear wins" claim to "bilinear, structured ridge and GP all beat recalibration, not tested pairwise (no ranking claimed among them)", and shows lab-trained structured residual zero-shot does not repair -40 C. Task 19 should map to exact CLM ids.
