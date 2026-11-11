# S7 viscosity model loose ends (v5)

Status: SCORED 2026-11-11 (readings declared before scoring; see Results). Declared 2026-11-11. Do not edit the readings after scoring; add a dated note.

## Questions

A. v4 compared bilinear, structured ridge and structured GP each against recalibrated Refutas only. Do they differ from each other?
B. The GP 90% intervals under-cover (0.77 dev, 0.69 locked). Can they be recalibrated, or should they be dropped?

## Data fixed in advance

Existing v4 outputs only, no refit: `results/v4/V/viscosity_predictions_v4.csv`, `viscosity_leaderboard_v4.csv`, and the README there. Dev = 45 fuels (leave-one-fuel-out), locked = 12 fuels, as in `splits/v4/real_fuels_v4_splits.json`. Conditions: overall, -20 C, -40 C (and others if present). Error metric: log-RMSE of ln(pred/meas) over fuel-resampled sets.

## Readings (fixed before scoring)

R-S7-1: pairwise paired bootstrap over fuels (5000 resamples, seed 0) of the log-RMSE difference for each pair (bilinear-ridge, bilinear-GP, ridge-GP), per split and condition. A pair "differs" if the 95% interval excludes 0. Report every pair's interval. No ranking is claimed among models whose pair intervals include 0. Dev and locked are never pooled.
R-S7-2: recalibrate GP intervals by split-conformal at fuel level: scale s = the ceil((n+1)*0.9)-th smallest of per-fuel max |residual|/predictive sd over dev fuels (leave-one-fuel-out predictions), per condition if at least 9 dev fuels, else one global factor. Intervals become pred +/- s*sd. Locked fuels are scored once with the dev-derived s.
R-S7-3: the recalibrated GP interval is "adopted" if dev coverage is at least 0.88 AND locked coverage is at least 0.80; mean interval width is reported next to the ridge interval built by the same recipe (no threshold). Otherwise the GP interval is dropped and ridge/bilinear intervals from the same recipe are reported instead.
R-S7-4: nothing is written to `configs/`.

A reading that cannot be evaluated is stated "not evaluable" with the reason.

## Results

Scored 2026-11-11 with `scripts/v5/s7_viscosity_pairwise.py` (existing v4 outputs, no refit). Files: `s7_pairwise_bootstrap.csv`, `s7_recalibration_factors.csv`, `s7_coverage_widths.csv`, `s7_verdict_rs73.csv`, `s7_reproduction_check.csv`. Tests: `tests/test_s7_viscosity.py` (3 pass).

Reproduction check: all 40 logged v4 numbers (GP 90% coverage and width per split/condition, and the paired-vs-recalibrated-Refutas dRMSE and CIs for bilinear/ridge/GP) were recomputed exactly (40/40 match; coverage 0.77 dev / 0.694 locked). The predictions file has no GP sd column, but gp90_lo/hi are symmetric in ln space around the GP prediction (max deviation 2e-16), so sd = (ln hi - ln lo)/(2*1.645) is recovered exactly; no retraining.

### R-S7-1 (pairwise, d = RMSE(first) - RMSE(second), 95% paired fuel bootstrap, 5000, seed 0)

| split | condition | bilinear - ridge | bilinear - GP | ridge - GP |
|---|---|---|---|---|
| dev | overall | +0.0075 [0.0020, 0.0132] differs | +0.0013 [-0.0080, 0.0091] | -0.0061 [-0.0167, 0.0012] |
| dev | 20 C | +0.0174 [0.0055, 0.0276] differs | +0.0132 [-0.0007, 0.0254] | -0.0043 [-0.0113, 0.0010] |
| dev | -20 C | +0.0108 [0.0010, 0.0232] differs | +0.0073 [-0.0049, 0.0200] | -0.0034 [-0.0173, 0.0059] |
| dev | -40 C | -0.0026 [-0.0179, 0.0133] | -0.0143 [-0.0325, 0.0022] | -0.0117 [-0.0243, -0.0007] differs |
| locked | overall | +0.0173 [0.0086, 0.0349] differs | +0.0153 [0.0065, 0.0310] differs | -0.0020 [-0.0057, -0.0004] differs |
| locked | 20 C | -0.0040 [-0.0294, -0.0016] differs | -0.0047 [-0.0332, -0.0019] differs | -0.0007 [-0.0053, -0.0001] differs |
| locked | -20 C | +0.0357 [0.0209, 0.0504] differs | +0.0294 [0.0125, 0.0464] differs | -0.0063 [-0.0124, -0.0007] differs |
| locked | -40 C | +0.0422 [0.0125, 0.0668] differs | +0.0433 [0.0163, 0.0659] differs | +0.0011 [-0.0031, 0.0079] |

Verdict: dev overall, ridge beats bilinear (intervals exclude 0); bilinear-vs-GP and ridge-vs-GP include 0, so no ranking between GP and either other model. At dev -40 C ridge beats GP (small, +-0.0117). Locked (12 fuels): bilinear is worse than ridge and GP overall, at -20 C and at -40 C; ridge and GP are separated by only 0.002 overall (ridge lower, interval barely excludes 0) and are indistinguishable at -40 C. The locked 20 C rows are unreliable: the interval for a negative-looking difference lies mostly below the point estimate (-0.0040 vs [-0.0294, -0.0016]), i.e. the bootstrap distribution is skewed by one or two fuels; do not read it as bilinear beating ridge. Net: ridge <= GP <= bilinear in point estimate in most cells, but only "bilinear is worse than ridge" holds on both dev (overall, 20, -20 C) and locked (overall, -20, -40 C) with intervals excluding 0; ridge and GP are not separable in any consistent way.

### R-S7-2 (split-conformal factors, dev fuels only, n = 45 per condition, rank k = ceil(46*0.9) = 42)

| interval | 20 C | -20 C | -40 C |
|---|---|---|---|
| GP (s multiplies sd) | 1.569 | 3.508 | 2.434 |
| ridge (s in ln units, sd=1) | 0.0745 | 0.1508 | 0.1498 |
| bilinear (s in ln units, sd=1) | 0.1035 | 0.1681 | 0.1334 |

The raw GP needed 1.57-3.51 times its stated sd (stated 1.645 sd). Ridge and bilinear have no predictive sd, so "the same recipe" is read as a constant-width interval pred +- s (score = |residual|) per condition; that is my interpretation, not specified in the readings.

### R-S7-3 (coverage / mean width in ln units, overall all T)

| interval | dev cov | locked cov | dev width | locked width |
|---|---|---|---|---|
| GP before (v4) | 0.770 | 0.694 | 0.212 | 0.247 |
| GP after (conformal) | 0.933 | 0.972 | 0.320 | 0.375 |
| ridge after | 0.933 | 0.833 | 0.250 | 0.250 |
| bilinear after | 0.933 | 0.694 | 0.270 | 0.270 |

Per condition, locked GP-after coverage: 0.917 (20 C), 1.000 (-20 C), 1.000 (-40 C). GP before at locked -20 C was 0.333.

Verdict: threshold met (dev 0.933 >= 0.88, locked 0.972 >= 0.80), so the recalibrated GP interval is "adopted" under the declared rule. Caveats that the rule does not capture: (i) dev coverage of 0.933 = 42/45 holds by construction per condition, so only the locked number is evidence; (ii) locked = 12 fuels, 36 correlated rows; 0.972 = 35/36, and the binomial uncertainty on a 90% target at n=12 fuels is large (0.80 threshold is ~one fuel); (iii) the GP interval is wider than the ridge interval (0.375 vs 0.250 locked) and ridge reached 0.833 (also passes the 0.80 locked bar, no threshold was set for it) while bilinear failed (0.694); the extra GP width is partly because GP sd varies by fuel and condition (-20 C s=3.5 sd). The ridge constant-width interval is a cheaper alternative with comparable locked coverage and 33% narrower; nothing in the readings picks between them, no choice made here.

### R-S7-4
Nothing was written to `configs/`.

### Caveats
Dev predictions are leave-one-fuel-out (each dev fuel predicted by a model fit on the other 44), so conformal exchangeability is approximate: the dev calibration residuals and the locked residuals come from different fits (locked from the all-45-dev model). Compositions are reconstructed, so interval widths include reconstruction error. Only 3 temperatures (20, -20, -40 C) exist; the conformal factor is not validated at other temperatures. Percentile bootstrap CIs with 12 fuels are crude.
