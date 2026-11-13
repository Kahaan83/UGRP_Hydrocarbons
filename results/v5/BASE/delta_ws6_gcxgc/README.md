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
