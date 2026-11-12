# C1 cold-viscosity correction from AV-33-22 (v5)

Status: SCORED 2026-11-12 (readings declared before scoring; see Results). Declared 2026-11-12. Do not edit the readings after scoring; add a dated note. These readings are self-declared in the same session as the code; git cannot prove they predate scoring.

## Question

With measured composition, raw Refutas under-predicts real-fuel kinematic viscosity at low temperature (S3 corrected: ln bias -0.27 at -20 C, -0.65 at -40 C; lab mixtures show about -0.11 below 263 K). AV-33-22 has measured viscosity at 9 temperatures (-40 to 100 C) for 38 fuels with GCxGC composition. Can a correction fitted on AV-33-22 remove most of this bias on other real fuels (AJFTD GCxGC fuels with cold viscosity)?

## Data fixed in advance

* Fit/within-source set: AV-33-22 (38 fuels) from `data/real_fuels/v5/crc_av3322_properties.csv` (all viscosity temperatures; where a temperature has both lab and certificate values use the median, as in S3, and report which). Composition mapping as in the corrected S3 (`results/v5/S3/corrected_2026-11-12_*`; zero-mass representatives dropped; non-finite predictions excluded and counted).
* Transfer set: AJFTD GCxGC fuels with a measured kinematic viscosity at -20 C and/or -40 C (S3 headline groups; the 11 S4-locked fuels reported apart, not in the headline). Double-count groups from the corrected S3 groups file.
* Base prediction M0: raw Refutas from measured composition at each temperature (the S3 route). No refit of any pure-component property.

## Corrections (all fitted on training fuels only)

* K0: none (M0).
* K1: per-temperature constant ln offset (mean of ln(meas/pred) at that temperature over training fuels).
* K2: ln offset linear in 1000/T, fitted across all of the training fuel's temperatures (a + b*1000/T), usable at any T.
* K3 (reference only, not a predictive method): Walther/ASTM D341 line through the fuel's own measured viscosities at 15-100 C extrapolated to the cold temperature. Uses measured viscosity of the scored fuel, so it is an upper bound and is never used in a verdict.

## Readings (fixed before scoring)

R-C1-1 (transfer): fit K1 and K2 on all AV-33-22 fuels, apply to AJFTD headline groups at -20 C and at -40 C separately. A correction "works" at a temperature if (a) the paired bootstrap (groups, 5000, seed 0) 95% interval of MAE(K0) minus MAE(K) excludes 0 on the positive side AND (b) |bias| after correction is at most 0.10 ln. Require at least 8 groups at that temperature; otherwise "n too small".
R-C1-2 (within source): leave-one-fuel-out inside AV-33-22 at -20 and -40 C, same rule. Report with the transfer result; within-source success alone does not count as working.
R-C1-3 (reverse check): fit K1/K2 on AJFTD headline groups, apply to AV-33-22; report the same metrics. A correction counts as source-robust only if both directions work.
R-C1-4 (margin context): report the fraction of fuels inside the S1 cold margins (0.4154 for T < 263 K, 0.6084 for 263-283 K) for K0 and for each K, per source. Nothing is written to `configs/`; any margin change is a proposal only.
R-C1-5 (decomposition): report K0 bias by temperature for AV-33-22 (9 points) so the warm-to-cold trend is visible, and K3's error at -20/-40 C as a reference for how far Refutas-at-40C-plus-slope could go.
R-C1-6: report the pure-component cause check if feasible without a download: for the 5 most abundant representatives in AV-33-22, compare the provider's pure ln nu at -20/-40 C with the Boehm 2022 pure data where it exists, stating whether the provider is extrapolating there; "not evaluable" if no overlap.

A reading that cannot be evaluated is stated "not evaluable" with the reason. Fuels already read (AV-33-22 in S2/S3/S5; AJFTD GCxGC in G1/L1/S3/S5) are labelled exploratory: this is a correction study, not a clean holdout.

## Results

Scored 2026-11-12 with `scripts/v5/c1_cold_viscosity_correction.py` (tests: `tests/test_c1_cold_viscosity.py`). Exploratory: every fuel here was already read in S2/S3/S5/G1/L1; this is a correction study, not a clean holdout. Sign: err = ln(pred) - ln(meas); the gain is MAE(K0) minus MAE(K), positive = correction better; paired bootstrap over groups, 5000 resamples, seed 0, percentile 95 %. The bootstrap resamples scored groups only; the uncertainty of the fitted offsets themselves is not propagated. Files in `results/v5/C1/`: `c1_metrics.csv`, `c1_bootstrap.csv`, `c1_readings.csv`, `c1_predictions.csv` (scored rows with K0/K1/K2/K3 errors), `c1_base_predictions.csv` (all fuel x T M0 predictions), `c1_fits.csv`, `c1_by_temperature.csv`, `c1_av3322_temperature_methods.csv`, `c1_pure_check.csv`.

### What was run (and choices the README left open)

* M0 is the S3 Refutas-from-measured-composition route (`predict_comp`, zero-mass representatives dropped), run at every measured temperature. AV-33-22: 38 fuels in 35 S3 groups, 9 temperatures. AJFTD: the 112 S3-eligible GCxGC fuels (unmapped mass <= 10 %) that have a measured viscosity; kinematic viscosity as in `targets(exclude_censored=True)` (reported temperatures only, snapped to the AV-33-22 temperature within 0.3 K). 32 AJFTD fuels have a usable viscosity row at a snapped temperature.
* Non-finite predictions: 1 of 423 fuel x T rows (AV-33-22 fuel CRC_AV3322_28 at 373.15 K); it is excluded from fitting and scoring and counted in `n_excluded_nonfinite`. It does not touch the cold cells.
* Median rule (S3): where a temperature has several values the median is used. Only -20 C has this (37 of 38 fuels have more than one value; methods D445, D7042, D7945, one "Not specified"); all other AV-33-22 temperatures are D7042 only (`c1_av3322_temperature_methods.csv`). I did not split lab from certificate values.
* Transfer set (R-C1-1): the S3 headline rule is one fuel per double-count group with the S4-locked groups removed, but I did NOT require the S3 reconstruction route (R), which cut AJFTD to 21 fuels in S3. So the AJFTD sets are the fuels with a measured cold viscosity: 8 groups at -40 C, 19 at -20 C. Per temperature the fuel with the most viscosity rows represents its group. Locked fuels (S4) are scored apart: 2 at -40 C, 9 at -20 C.
* Fits: transfer uses all 38 AV-33-22 fuels (341 rows). Within-source: leave-one-group-out, the whole S3 group of the scored fuel is held out (35 fits); AV-33-22 scored on one fuel per group. Reverse: fitted on AJFTD non-locked groups, one fuel per group per temperature (19 fuels, 50 rows; K1 has an offset only at the AJFTD temperatures: 8 groups at -40 C, 19 at -20 C), scored on the 35 AV-33-22 groups. K3 needs at least 2 measured viscosities at 15-100 C for the scored fuel.
* Success rule applied as declared: n >= 8 groups, bootstrap interval of the gain excludes 0 on the positive side, and |mean err after correction| <= 0.10 ln.

### R-C1-1 transfer: fitted on AV-33-22, applied to AJFTD headline groups

| T | corr | groups | bias K0 | bias K | MAE K0 | MAE K | gain (95 % CI) | in S1 margin K0 -> K | works |
|---|---|---|---|---|---|---|---|---|---|
| -40 C | K1 | 8 | -0.672 | -0.053 | 0.672 | 0.094 | +0.577 (0.520, 0.619) | 0.00 -> 0.875 | yes (n = 8, the minimum) |
| -40 C | K2 | 8 | -0.672 | -0.281 | 0.672 | 0.281 | +0.390 (0.390, 0.390) | 0.00 -> 0.875 | no: abs bias 0.28 > 0.10 |
| -20 C | K1 | 19 | -0.279 | -0.025 | 0.357 | 0.134 | +0.223 (0.166, 0.254) | 0.79 -> 0.947 | yes |
| -20 C | K2 | 19 | -0.279 | -0.006 | 0.357 | 0.125 | +0.232 (0.171, 0.269) | 0.79 -> 0.947 | yes |

The K2 -40 C interval has zero width because every error stays negative after correction, so the gain is the same shift for every group. S4-locked fuels (apart, same fits): -20 C n = 9, K1 bias -0.065 and K2 -0.047 (both meet the rule; K0 bias -0.319); -40 C n = 2 "n too small" (K1 bias -0.155, K2 -0.384).

### R-C1-2 within AV-33-22 (leave-one-group-out, 35 groups)

| T | corr | bias K0 | bias K | MAE K0 | MAE K | gain (95 % CI) | in margin K0 -> K | works |
|---|---|---|---|---|---|---|---|---|
| -40 C | K1 | -0.637 | -0.018 | 0.637 | 0.111 | +0.526 (0.477, 0.570) | 0.09 -> 1.00 | yes |
| -40 C | K2 | -0.637 | -0.247 | 0.637 | 0.248 | +0.389 (0.386, 0.391) | 0.09 -> 0.89 | no (abs bias 0.25) |
| -20 C | K1 | -0.264 | -0.010 | 0.264 | 0.062 | +0.203 (0.175, 0.227) | 0.97 -> 1.00 | yes |
| -20 C | K2 | -0.264 | +0.008 | 0.264 | 0.062 | +0.202 (0.171, 0.230) | 0.97 -> 1.00 | yes |

### R-C1-3 reverse: fitted on AJFTD headline groups, applied to AV-33-22 (35 groups)

| T | corr | bias K0 | bias K | MAE K0 | MAE K | gain (95 % CI) | in margin K0 -> K | works |
|---|---|---|---|---|---|---|---|---|
| -40 C | K1 | -0.637 | +0.035 | 0.637 | 0.111 | +0.527 (0.466, 0.581) | 0.09 -> 1.00 | yes |
| -40 C | K2 | -0.637 | -0.142 | 0.637 | 0.163 | +0.474 (0.453, 0.492) | 0.09 -> 1.00 | no (abs bias 0.14) |
| -20 C | K1 | -0.264 | +0.015 | 0.264 | 0.060 | +0.204 (0.171, 0.233) | 0.97 -> 1.00 | yes |
| -20 C | K2 | -0.264 | +0.074 | 0.264 | 0.088 | +0.177 (0.134, 0.216) | 0.97 -> 1.00 | yes |

### Verdict per correction x temperature x direction (success rule)

| correction, T | transfer AV->AJFTD | within AV (LOGO) | reverse AJFTD->AV | source-robust (transfer and reverse) |
|---|---|---|---|---|
| K1, -20 C | works (n = 19) | works | works | yes |
| K1, -40 C | works (n = 8, at the minimum) | works | works | yes (the AJFTD side rests on 8 groups) |
| K2, -20 C | works (n = 19) | works | works | yes |
| K2, -40 C | not met (bias -0.28) | not met (-0.25) | not met (-0.14) | no |

### R-C1-4 margin context

Fraction of fuels inside the S1 two-sided margin (abs err <= 0.4154 for T < 263 K; 0.6084 for 263-283 K), K0 -> K1 / K2:
* AJFTD headline: -20 C 0.79 -> 0.947 / 0.947 (n = 19); -40 C 0.00 -> 0.875 / 0.875 (n = 8).
* AV-33-22: -20 C 0.97 -> 1.00 / 1.00; -40 C 0.09 -> 1.00 (K1) / 0.89 (K2) within-source; reverse -40 C 1.00 for both. K0 is nearly always inside the margin at -20 C for AV-33-22, so the raw margin was not the problem there; at -40 C raw Refutas is outside the margin for 91 % of AV-33-22 fuels and every AJFTD headline fuel.
* Nothing was written to `configs/`. A margin change is only a proposal and is not made here: the corrected-model residual spread would need its own conformal calibration on data not used for the offset.

### R-C1-5 decomposition: K0 bias by temperature (AV-33-22, 35 groups; `c1_by_temperature.csv`)

| T (C) | -40 | -20 | 0 | 15 | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| bias ln | -0.637 | -0.264 | -0.055 | +0.019 | +0.035 | +0.066 | +0.058 | +0.035 | -0.007 (34 groups) |

So Refutas is slightly high (up to +7 %) from 15 to 80 C and turns sharply low below 0 C. The warm-to-cold change is not linear in 1000/T, which is why K2 (a single straight line in 1000/T over all nine temperatures) leaves -0.25 at -40 C in the leave-one-group-out run, while per-temperature K1 does not. AJFTD K0 shows the same shape where it has data: -0.672 (-40 C, 8 groups), -0.279 (-20 C, 19), -0.075 (0 C, 3), +0.011 (40 C, 10).

K3 (reference only, never a verdict): Walther line through the fuel's own 15-100 C viscosities. AV-33-22 (35 groups): error +0.136 at -40 C and +0.075 at -20 C (it over-predicts). AJFTD: only 4 groups at -40 C (error +0.194) and 6 at -20 C (+0.070) have two or more warm viscosities, so it is not comparable there. Its MAE is about as large as K1's (0.111 and 0.062 within source) while using the scored fuel's own measurements, so Refutas plus a fitted offset is not worse than this extrapolation upper bound.

### R-C1-6 pure-component check (`c1_pure_check.csv`)

The 5 most abundant AV-33-22 representatives (mean mass share 11.8, 8.9, 5.9, 5.5, 5.2 %): decalin, 2-methyldodecane-type iso-alkane (SMILES CCCCCCCCC(C)CC), butylcyclohexane, tetralin, hexylcyclohexane. Boehm 2022 pure kinematic viscosity (`boehm2022_viscosity.csv`, on disk; nothing downloaded) exists for only two of them: decalin at -20 C and butylcyclohexane at -20 and -40 C; the other cells are "not evaluable" (no Boehm pure data at that temperature). Provider ln(pred/Boehm): decalin -20 C -0.075; butylcyclohexane -20 C -0.177 and -40 C -0.460 (provider 4.39 cSt vs Boehm 6.96 cSt). The provider is a fit to measured pure data whose range starts at 253.2 K for butylcyclohexane (so -40 C is 20 K outside the range, inside its 25 K tolerance) and at 243.1 K for decalin; for the iso-alkane (range from 283.1 K), tetralin and hexylcyclohexane (293.1 K) both -20 and -40 C lie more than 25 K below the data, so the provider switches to its ML predictor there, i.e. it is extrapolating. The pure-component error at -40 C is large (-0.46 ln) for the one abundant compound that can be checked there, in the same direction as the fuel bias; three data points cannot establish the cause of the fuel-level bias (-0.64), and three of the five dominant compounds cannot be checked at all.

### Caveats

* Exploratory (flag `exploratory = True` in the metrics file). The fuels were already read, and the cold bias of AV-33-22 and AJFTD is similar, which is what K1 exploits.
* K1 removes the mean bias at -40 C and -20 C by construction (the offset is the mean raw error at that temperature). The informative parts are that the spread also falls (MAE 0.64 -> 0.11 within source at -40 C) and that the offsets transfer across sources in both directions (AV-33-22 offsets leave bias -0.053 at -40 C and -0.025 at -20 C on AJFTD; AJFTD offsets leave +0.035 and +0.015 on AV-33-22).
* K1 cannot be used at a temperature absent from the fit; K2 can, but it fails the rule at -40 C.
* The AJFTD -40 C cell has exactly 8 groups; in the reverse direction the AJFTD K1 offset is fitted on the same 8.
* Offsets are not stored in `configs/`; the margin and bias numbers above are not applied anywhere.

### Notes (dated)

* 2026-11-12: AJFTD sets here are not restricted to the S3 R-paired fuels (S3 needed a reconstructed composition; this study does not). The README says "S3 headline groups"; I read it as S3's grouping and locked rules, not its R requirement.
