# delta_ws6_gcxgc (WS6): real fuels from GCxGC class x carbon-number distributions

Declared 2026-11-14 before `scripts/v5/base/delta_ws6_gcxgc.py` exists; this commit contains no result file.
Worktree branch `worktree-agent-ab5aceb174f266a9b`, HEAD 6204538. Label for every number below: `exploratory / dev` (all T2 fuels were read before this declaration;
the 12 CRC 647 "locked" fuels have been read five times).
Nothing in `configs/`, `predict_fuel`, the engines or `verification:` flags is changed. Any margin or model designation is "needs user decision".

## Question
Does a model fed the measured class x carbon-number matrix of a real fuel beat the current physics (mapped representatives), beat physics plus a constant offset, and beat a constant,
for density, flash, kinematic viscosity at -20 and -40 C, freezing and smoke? And does it give a usable one-sided cold-viscosity or freezing margin on real fuels?

## Current status being tested (checked against the repo, not taken from v4 text)
* Cold viscosity: raw Refutas is biased by about +0.6 to +0.7 ln at 233-253 K (scoreboard MAE 0.65 at 233.15 K). The `inf` is the yaml `raw_refutas` margin below 283 K
  (`cells_t2_bulk.py`: `m90 = inf if T-bin != T_ge_283K`), shown in every T2 scoreboard cell at 233.15 K. The K1 cold offset (a constant) already brings MAE to 0.11 ln
  (`refutas_raw+k1`, 48 groups). `delta_t2_residual` found that a ridge residual wins on cold viscosity only through the intercept. So the K1-style offset, not raw physics,
  is the comparator that matters for viscosity.
* Freezing: the v4 statement "no skill over a constant" is out of date. `delta_freezing` / S5: F1x (2-parameter map from first-species-to-freeze temperature) has pooled MAE 3.07 C (n = 49)
  and beats the leave-one-group-out constant pooled and on AV-33-22 (interval above 0); AJFTD alone has 12 groups (not evaluable). The ideal SvL designated route does not (S5 `I`). F1x is the declared physics for freezing.
  Its `a`, `b` are fitted leave-one-group-out in S5, so its predictions are out of sample for the scored fuel.
* Smoke: AJFTD bulk regression (`smoke_point_from_bulk`, aromatics vol %) was fitted on 100 AJFTD fuels; the YSI route returns no value on AV-33-22 (coverage < 0.90).

## Data (fixed now)
* Fuels: T2 entities with `chosen_for_group = True` in `results/v5/BASE/t2_fuels.csv` (146 groups: AJFTD GCxGC 113 entities, AV-33-22 38, FuelLib 4; the S3 double-count groups, so the shared fuels
  10325, 10289, 4658, 11769 count once). Targets: `t2_targets.csv` medians (censored values excluded already). `low_coverage` fuels excluded (one fuel).
* Properties and temperatures: density at 288.15 K; flash (all methods, in one pooled cell, method family as stratum for the constant); ln kinematic viscosity at 253.15 K and 233.15 K;
  freezing (D5972 etc., uncensored); smoke (D1322; AV-33-22 and FuelLib from `t2_targets.csv`, AJFTD from `ajftd_smoke_naph.csv` rows with no qualifier).
* Physics predictions are reused, not recomputed: density `ideal_volume` / designated, viscosity `refutas_raw` / designated (`rows_t2_bulk.csv`), flash `pred_designated` (`rows_t2_flash.csv`),
  freezing `F1x` (`s5_predictions.csv`), smoke `smoke_point_from_bulk(measured aromatics vol %)`. A spot check that the reused values reproduce the scoreboard MAEs is part of step A.
* Smoke in-sample rule: AJFTD fuels whose `fuel_id` appears in the file the bulk regression was fitted on (`data/processed/v4/sooting/smoke_point_ajftd.csv`) are in-sample for the physics baseline
  and are excluded from the smoke gate (counted and listed). Smoke fuels with aromatics outside 8-23 vol % are scored but flagged out of domain.
* CRC 647 (57 jet fuels) has class types and a reconstructed carbon-number split, not GCxGC. It never enters training. It is an application tier only (step T5), labelled `fold_b / dev`, and skipped
  if the per-fuel physics predictions in `results/v4/R3` cannot be aligned.

## Fuel family grouping
Metadata check: `ajftd_fuels_wide.csv` / `ajftd_fuels_long.csv` carry no fuel-type column (all 39 AJFTD fuels with targets have `source_group = NUMERIC_ID`; `category` is a metric type; AV-33-22 has Jet A / Jet A-1 only).
Fallback, declared: family = source, with FuelLib merged into AJFTD (same AFRL GCxGC program; three of its four fuels are AJFTD fuels).
* **Gate split (G):** repeated grouped K-fold, K = 5, 20 repeats, unit = S3 double-count group, both sources in the training folds. A fuel's CV error is the mean over repeats of its out-of-fold error.
* **Transfer split (T, secondary, not a gate):** leave-one-source-out, AV-33-22 <-> AJFTD(+FuelLib). Largely run already by `delta_t2_residual`; reported for comparison.
Standardisation, offsets, PLS components, GP hyperparameters and NN weights are fitted inside the training fold only.

## Feature schema (fixed now)
Nine harmonised classes: n_paraffin, iso_paraffin, monocycloparaffin, dicycloparaffin, tricycloparaffin, alkylbenzene, cycloaromatic, diaromatic, alkene (AV-33-22 names n-alkane, iso-alkane, monocycloalkane,
dicycloalkane, tricycloalkane map one-to-one; AJFTD `diene` is added to alkene). Ten carbon bins per class: <= 8, 9, 10, ..., 16, >= 17 (open-ended "plus" bins and C > 16 go to >= 17).
Mass % renormalised over the mapped total to 1. 90 features (p >> n), so compression is mandatory: PLS (below). Fuels with mass outside the nine classes or a closure outside 98-102 % are excluded and counted.

## Models (fixed now; one primary model per property)
For target y (density, flash, freezing, smoke raw; ln nu for viscosity) with physics prediction p:
* C0 constant: training-fold mean of y (flash and smoke strata: none; density and viscosity are single-temperature cells).
* C1 physics: p.
* C2 physics + offset: p + mean(y - p) over training fuels (the K1-style offset; no composition information).
* **M1 (primary): C2 + PLS residual.** PLS (scale = True) on the 90 features against (y - p), components chosen from 1..5 by inner 5-fold group CV on the training fold, applied to the centred residual.
* M2 (secondary): C2 + GP residual on the PLS scores (ARD RBF x constant + white noise, 2 restarts, hyperparameters by marginal likelihood on the training fold).
* M3 (secondary): C2 + small NN residual on the PLS scores (MLP, one hidden layer of 8 tanh units, L2 alpha 1.0, lbfgs, max 500 iterations, mean of 10 seeds). No tuning on any test fold.
* M0 (secondary, no physics): direct PLS regression of y on the features.

## Safety directions (fixed now)
* Flash: over-prediction unsafe (rate of error > +5 C, worst over-prediction).
* Cold viscosity (-20, -40 C): under-prediction unsafe. Margin = one-sided group conformal margin (`fuelblend.uq.conformal.one_sided_margin`, direction `under`, alpha 0.10) on out-of-fold errors;
  `inf` if too few groups. Unsafe rate = share of fuels under-predicted by more than the margin of the same comparator (and by 0.4154 ln, the S1 margin below 263 K, as a fixed yardstick).
* Freezing: spec is a maximum, so under-prediction of the freezing point is unsafe (same margin routine, direction `under`).
* Smoke: spec is a minimum, so over-prediction is unsafe (direction `over`).
* Density: two-sided (|error|).

## Readings and gate (fixed now)
For each cell and each model Mi: d_k = MAE(Ck) - MAE(Mi) for k in {C0, C1, C2}, paired bootstrap over fuels, 5000 resamples, seed 0, 95 % percentile interval (positive = Mi better).
* **R-WS6-1 (gate, per primary cell):** M1 passes iff all three intervals (vs C0, C1, C2) lie above 0 on split G, **and** for flash the over-prediction rate (> +5 C) is not higher than C2's,
  **and** for the one-sided cells the unsafe rate at the cell's own out-of-fold margin is not higher than C2's.
* **R-WS6-2 (multiplicity):** six primary cells (density, flash, nu(-20), nu(-40), freezing, smoke). A pass counts as established only if the one-sided bootstrap p-value of its weakest comparator also survives
  Holm over the six cells at 0.05. Unadjusted passes are reported as such. M2, M3 and M0 are secondary: they need the same three intervals at 98.3 % (Bonferroni over three models) and are never called established.
* **R-WS6-3 (cold-viscosity margin):** "a real-fuel out-of-fold one-sided margin is finite and the unsafe-side coverage on split T (margin set on one source, applied to the other) is at least 90 % in both directions",
  evaluated for C2 (the offset, what could ship) and M1, at 253.15 K and 233.15 K. If it holds it is a proposal only. `proposed_margins.yaml` is written in this folder.
* **R-WS6-4 (freezing):** "a composition-aware model (M1) beats F1x + offset (C2) and the constant on split G with intervals above 0". Failure is a valid, reported outcome.
* **R-WS6-5 (physics vs offset):** descriptive: for each cell, whether C2 beats C1, i.e. how much of the physics error is a constant.

## Power (stated before any model is fitted)
Step A computes, from baseline errors only (C2 vs C0 per fuel, which involves no ML model), the standard deviation of the per-fuel paired difference of |error|, and from it the minimum detectable
MAE difference for 80 % power at two-sided 5 % (2.8 * sd / sqrt(n)), per cell. Written to `power.csv` before `oof_predictions.csv` exists. The README results state which gate differences were below that floor.
Samples are small (n about 40-75 per cell, 40-50 for freezing and smoke); a null here means "not shown at this n", not "no effect".

## Published comparison (only what is already in the repo or openly fetchable)
In repo: Bell, Boehm and Heyne 2025 freezing MAE 4.4 C (roadmap text, from the paper), against F1x LOGO MAE 3.07 C on 49 fuels (different data and split; not like-for-like). Illinois deep learning on GCxGC (106 samples),
Dayton GCxGC-VUV (ten properties) and DLR iso-alkane subgrouping numbers are NOT in `data/raw`; they are compared only if an open PDF is fetched and its numbers read (recorded in the results), otherwise "not compared".

## Needs user decision (nothing here is adopted)
* Whether any M1/M2/M3 variant becomes a shipped or designated variant of `predict_fuel`.
* Whether any proposed real-fuel cold-viscosity or freezing margin replaces the `inf` yaml margin.

---

## Result 2026-11-14 (appended after the run; the declarations above are unchanged)

Run: `scripts/v5/base/delta_ws6_gcxgc.py` (`table`, `base`, `models`, then `summary`: a bug in the secondary-interval line stopped `models` after the 15-minute fits, so `oof_G.csv` / `oof_T.csv` were reused, no refit),
`delta_ws6_crc647.py` (application tier), `delta_ws6_posthoc.py` (post-hoc control). Files in this folder: `fuel_table.csv`, `features.csv`, `power.csv`, `summary_G.csv`, `summary_T.csv`, `margins.csv`, `readings.json`,
`crc647_application.csv`, `posthoc_source_offset.csv`, `proposed_margins.yaml`. All numbers are dev / exploratory.

### Data actually used
145 chosen fuels (the declaration said 146; one low-coverage fuel is dropped), 142 with a usable 9-class matrix. FuelLib 10264, 10289 and 10325 are excluded and counted: FuelLib stores only coarse classes (cycloalkane, aromatic), so under 97 % of their mass
maps to the nine-class schema (their groups AJFTD:39 / AJFTD:122 are therefore absent from the model cells). Cell sizes: density 15 C 61 (35 AV-33-22, 26 AJFTD), flash 67, nu(-20 C) 62, nu(-40 C) 45 (only 10 AJFTD), freezing 54, smoke 59.
Smoke in-sample rule: 0 chosen AJFTD fuels overlap the 100 fuels the bulk regression was fitted on, so nothing was excluded. 14 of 59 smoke fuels have aromatics outside 8-23 vol % (flagged; every model errs by about 8 mm there).
Clarification, not a gate change: a fuel's CV error is the mean over the 20 repeats of its absolute out-of-fold error (signed mean for bias and margins). Spot check of the reused physics: C1 MAE for nu(-40 C) is 0.646 ln (scoreboard pooled 0.650) and flash 2.56 C, consistent with BASE.

### Split G (gate). MAE, then interval of MAE(comparator) - MAE(M1) (95 %; positive = M1 better)

| cell | n | C0 const | C1 physics | C2 phys+offset | M1 PLS | M1 vs C0 | M1 vs C1 | M1 vs C2 | unsafe side, M1 vs C2 | R-WS6-1 |
|---|---|---|---|---|---|---|---|---|---|---|
| density 15 C (kg/m3) | 61 | 16.39 | 4.75 | 3.50 | 2.28 | +14.1 (11.5, 16.9) | +2.47 (1.77, 3.17) | +1.22 (0.66, 1.78) | two-sided | PASS (Holm survives) |
| flash (C) | 67 | 5.13 | 2.56 | 2.39 | 2.06 | +3.07 (2.07, 4.18) | +0.50 (-0.18, 1.17) | +0.33 (-0.35, 1.03) | over >5 C: 7.5 % vs 6.0 % | FAIL |
| ln nu, -20 C | 62 | 0.208 | 0.296 | 0.088 | 0.064 | +0.144 (0.108, 0.183) | +0.232 (0.189, 0.265) | +0.024 (0.004, 0.041) | under >0.4154: 0 vs 0 | PASS (weakest p 0.011, Holm survives) |
| ln nu, -40 C | 45 | 0.223 | 0.646 | 0.110 | 0.047 | +0.176 (0.131, 0.222) | +0.599 (0.561, 0.638) | +0.063 (0.042, 0.085) | under >0.4154: 0 vs 2.2 % | PASS (Holm survives) |
| freezing, F1x (C) | 54 | 6.11 | 3.55 | 3.63 | 4.30 | +1.81 (0.08, 3.74) | -0.75 (-1.56, -0.05) | -0.68 (-1.50, 0.02) | not better | FAIL |
| smoke (mm) | 59 | 4.44 | 3.22 | 3.40 | 2.93 | +1.51 (0.34, 2.66) | +0.29 (-0.62, 1.22) | +0.47 (-0.34, 1.24) | not better | FAIL |

* The declared clause "unsafe rate at the cell's own out-of-fold margin not higher than C2's" is near-vacuous (about 10 % by construction, 6-9 % at these n). The effective guards are the fixed yardsticks in the table and the margin size (below).
  Flash: M1's over-prediction rate above 5 C is higher (7.5 % vs 6.0 %) but its worst case is shorter (+7.7 C vs +18.0 C for C2); flash fails the interval part anyway.
* R-WS6-2: three cells pass with Holm (density, nu at -20 C, nu at -40 C); flash, freezing, smoke do not. M2 (GP) and M3 (NN) pass the 98.3 % intervals on the same three cells only. M0 (PLS without physics) is worse than C2 on flash, freezing and smoke.
* R-WS6-5: C2 beats C1 on density, flash and both viscosities (viscosity almost entirely: the raw Refutas bias is -0.27 ln at -20 C and -0.65 ln at -40 C). C2 does not beat C1 on freezing or smoke.
* Freezing status correction: F1x beats the constant, MAE 3.55 vs 6.11 C, interval of the difference +2.56 (0.89, 4.64), 54 fuels. The "no skill over a constant" statement is true of the v4 designated route, not of F1x.
  M1 is significantly worse than F1x itself (-0.75, interval below 0).
* Smoke (descriptive): on the 45 in-domain fuels M1 MAE is 1.24 mm against physics 1.73 and offset 2.09; the 14 out-of-domain fuels carry about 8 mm error for every model.

### Does the gain survive a change of source? (split T, secondary, not a gate)
No. Training on one source and testing on the other, M1 is worse than C2 (offset only) in 10 of 12 source x cell comparisons (density, flash, both viscosities, freezing, on both sources). Examples: nu(-40 C) 0.202 vs 0.110 (AJFTD) and 0.215 vs 0.112 (AV-33-22);
density 9.7 vs 6.2 and 6.4 vs 3.4 kg/m3; flash 3.3 vs 1.4 C on AV-33-22. M1 is better only for smoke (5.20 vs 5.57 AJFTD, 2.57 vs 4.86 AV-33-22), and there the in-domain physics alone gives 1.33 mm on AV-33-22.
GP / NN transfer a little better in places (nu(-40 C) on AV-33-22: GP 0.085) but not consistently. Same pattern as `delta_t2_residual`.
Post-hoc control (not declared; labelled post hoc): against C2s = physics + per-source offset fitted in-fold, M1 still wins on split G for density (+0.79, 0.29 to 1.34), nu(-20 C) (+0.025, 0.006 to 0.041) and nu(-40 C) (+0.068, 0.047 to 0.089).
So the split-G gain is not source identity alone. What remains is within-population composition signal or near-duplicate batches inside a source (S3 groups remove only identical shared fuels); split T gives no support for a transferable signal.

### CRC 647 application tier (fold_b / dev; reconstructed r3 compositions; models trained on all 142 GCxGC fuels; 57 fuels)
Flash MAE physics 3.97, offset 4.10, M1 4.00 C (over >5 C: 15.8 %, 21.1 %, 14.0 %). nu(-20 C): physics 0.363, offset 0.103, M1 0.116. nu(-40 C): physics 0.739 (all 57 fuels under-predicted by more than 0.4154), offset 0.128, M1 0.286 (12 % under by more than 0.4154), constant 0.223.
The real-fuel cold offset transfers to reconstructed compositions; the learned composition term does not. Not a locked test (the 12 locked fuels have been read five times).

### Cold-viscosity margin and freezing (R-WS6-3, R-WS6-4)
* The `inf` is the yaml raw-Refutas margin below 283 K. On real GCxGC fuels the raw physics margin (one-sided, 90 %) is finite but large: 0.405 ln at -20 C and 0.859 ln at -40 C, because of the -0.27 / -0.65 ln bias.
* With the K1-style offset (C2): pooled one-sided margin 0.135 ln (-20 C) and 0.215 ln (-40 C). At -40 C the margin set on AV-33-22 (0.215) covers 90 % of AJFTD groups and the margin set on AJFTD (0.432) covers 100 % of AV-33-22, so R-WS6-3 HOLDS for C2 at -40 C (AJFTD has only 10 groups).
  At -20 C it does not hold: the AV-33-22-set margin (0.123) covers 77.8 % of the 27 AJFTD groups; the AJFTD-set value 0.186 covers all AV-33-22 groups. M1 gives smaller margins (0.053, 0.074) but fails the transfer check (88.9 %, 80 %). Proposals are in `proposed_margins.yaml`, not in `configs/`.
* Freezing (F1x, under-prediction unsafe): offset margin 5.3 C pooled, 7.1 C set on AJFTD, transfer coverage 80 % (not shown to transfer); M1 margin 8.2 C. R-WS6-4 fails.

### Power
`power.csv` was written before any model was fitted (SD of the per-fuel paired difference of C0 and C2 errors): minimum detectable MAE difference at 80 % power, two-sided 5 %: density 4.2, flash 2.0, nu(-20 C) 0.046, nu(-40 C) 0.045, freezing 2.8, smoke 0.69.
That proxy is too coarse for M1 vs C2 (differences between similar predictors have a smaller SD; density and viscosity differences below it were still resolved). The practical resolution is the half-width of the M1-vs-C2 interval: about 0.7 C flash, 0.8 C freezing, 0.8 mm smoke, 0.6 kg/m3 density, 0.02 ln viscosity.
Flash, freezing and smoke nulls mean "no gain larger than about that shown at this n", not "no effect". Cell sizes 45-67; AJFTD at -40 C has 10 fuels.

### Published numbers
Qualitative only. A web search found the open Dayton paper (Heyne et al., Fuel 2022; DOT ROSAP 64804/64805); its numbers were not read from the PDF. Per its abstract-level summary, density and flash point (small isomer variance) are predicted better from GCxGC than freeze point and viscosity.
Our split-G result matches that ordering for flash/density versus freezing, but not for viscosity, where a constant offset of the physics already captures most of the error. Illinois and DLR numbers were not compared. Bell, Boehm and Heyne freezing MAE 4.4 C (roadmap text) vs our F1x 3.55 C (54 fuels): different data and splits, not like-for-like.

### Readings summary
R-WS6-1 passes for density, nu(-20 C), nu(-40 C); fails for flash, freezing, smoke (Holm leaves the three passes standing). R-WS6-3 holds for the offset model at -40 C, not at -20 C, not for M1. R-WS6-4 fails. R-WS6-5: the constant offset explains most of the physics error for viscosity.
The robust finding is the constant cold offset and its margin, not the composition model: M1's split-G density and viscosity gains do not transfer across sources or to reconstructed CRC 647 compositions.

### Needs user decision (nothing adopted; `configs/`, `predict_fuel` and the engines are untouched)
1. Whether to replace the `inf` cold-viscosity yaml margin (and the K1 variant's two-sided 0.172 / 0.293) with the one-sided real-fuel offset margins in `proposed_margins.yaml` (0.135 / 0.215 ln pooled; 0.186 / 0.432 ln if the conservative AJFTD-set value is preferred; -20 C does not pass the transfer check).
2. Whether any of M1-M3 becomes a variant of `predict_fuel`. Evidence here argues no (does not transfer across sources).
3. Whether to correct the v4 "freezing: no skill over a constant" lines for F1x (noted in the `docs/model_decision.md` WS6 section only).
