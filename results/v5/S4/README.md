# S4: held-out evaluation of bulk-descriptor models on AJFTD (v5)

**DECLARED BEFORE SCORING.** This file and `splits_ajftd_v5.json` were written before any model was fitted or scored.
The `dev` stage records the SHA-256 of both; the `locked` stage refuses to run if either changed, and refuses to run
twice (`locked_scores.csv` / `locked_run.json` guard). Hyper-parameters below are fixed; there is no tuning.

Inputs: `data/real_fuels/v5/ajftd_fuels_wide.csv` (+ `_long.csv` for viscosity temperature flags and flash method).
Provenance class `real_fuel_bulk_descriptors`: no molecular composition, bulk descriptors only.

## 1. Features, targets, models
* Features (only these): D86 T0, T10, T50, T90, T100 (C), aromatics vol%, density at 15 C. Sparse D86 points (T5, T20, T95 ...) are not used.
  No other target and no other-temperature viscosity is ever a feature.
* Targets: flash point (C), freezing point (C), kinematic viscosity at -20 C (cSt), net heat of combustion (MJ/kg).
  Eligible row = target present and at least 3 of the 7 features present. All models are scored on the same rows per target.
* **Viscosity target is rebuilt from the long table**: `T_K == 253.15` and `T_source in {reported, reported_parsed_text}`, median per fuel.
  A fuel with no such row but with `assumed_spec_-20C` rows is "assumed-T"; those fuels (578 in dev, 0 in locked) are excluded from fitting
  and from every verdict and are scored separately (`*_assumed_visc.csv`, models fitted on measured dev rows only).
  Verdict scale for viscosity: MAE in cSt (models fit ln(cSt) and are back-transformed; ln-scale not used for the verdict).
* Models: `const` (training mean, raw scale); `ridge` (fold-internal median imputation + missing indicators + standardise, alpha = 1);
  `hgb` (sklearn HistGradientBoostingRegressor, 200 iter, lr 0.05, depth 3, min_samples_leaf 20, l2 1, random_state 0).
  Published correlations, unfitted: flash from D86 T10 by the API Technical Data Book 2B7.1 form
  `1/T_F = -0.024209 + 2.84947/T10 + 0.0034254 ln T10` (K, T10 = D86 10 % point; coefficients quoted from memory, UNVERIFIED against
  the book; sanity check T10 = 175 C gives 53 C) and Butler (1956) `FP[F] = 0.69 T10[F] - 118.2` (also from memory). Rows with no T10 fall back to the training mean (count reported).
  **No correlation for viscosity at -20 C** (standard density/distillation viscosity correlations are defined at 37.8/98.9 C, extrapolation to -20 C is not defensible)
  **and none for heat of combustion** (no form I can cite verbatim: ASTM D3338/D4529 need inputs not in this table). Freezing: no published bulk-descriptor correlation used.
* **Smoke point is not in the wide table: skipped.** v4's smoke null (CLM-30) is therefore not revisited here.

## 2. Splits (`splits_ajftd_v5.json`)
* "Source group" for splitting and bootstrapping = **source_group x test_year ("cell")**, 45 cells. The 10 raw source_groups cannot support the
  analysis: a locked 20 % would be 2 groups, a bootstrap over them is meaningless, and conformal needs > 9 groups. This is a declared deviation;
  consecutive IDs within a source across adjacent years may still be correlated, so leakage across cells of the same source_group is possible.
  A dev-only sensitivity (`logo_coarse_source_group`) holds out whole raw source_groups.
* Locked: cells sorted by `sha256(seed|cell)`, first ceil(0.2 x 45) = 9 cells, seed `S4-v5-split-seed-1`, no re-rolling.
  Resulting locked set is **264 of 2468 deduplicated fuels (10.7 %)** and **69 % of it is one cell (JA|2011, 182 fuels)**; JA|2011 is a single
  source/year, so locked conclusions are about that regime. JAA, JFK, LAB_REPORT 2015 and most of the large cells are dev.
  Per target the locked set has 7-8 cells; a locked verdict with fewer than 5 groups would read "not decidable".
* Dedup (108 raw fuel ids merged): identical numeric vectors (>= 6 values) within a source_group collapse to the lowest fuel_id; `X_<digit>` labels whose base label
  exists in the same source_group collapse onto it; remaining members merged by median. Identical vectors spanning dev/locked would be forced to dev (0 found).
* Leave-group-out CV (LOGO) on dev cells = the verdict scheme. Within-group error = random 5-fold over dev fuels (repeat batches leak; reported to show the inflation).
* Time split (dev only): test = latest years until >= 25 % of dev fuels with a known year (2015, 2016, 2017, 2021); train earlier years.
  Confounded with source: JFK, LAB_REPORT and PAREN_BATCH dominate the late years.
* JFK freezing is 292 of 353 identical -43.5 C (D2386) and JFK heat 43.2 MJ/kg (D3338) is near constant: probable reporting defaults.
  A sensitivity subset `excl_JFK` is reported for dev and locked. The main verdict uses all rows.

## 3. Declared reading (the only criterion)
> A model beats the constant-mean baseline for a target iff the paired-bootstrap 95 % CI (5000 resamples of source groups = cells, seed 0) of
> MAE(model) - MAE(const) **excludes 0 on the negative side on leave-group-out dev AND on the locked set**.

* Dev-only success or locked-only success is not "beats". A locked verdict with fewer than 5 groups is "not decidable".
* v4 stated "freezing models show no skill over a constant" (CLM-16, CLM-38, CLM-41; on 57 reconstructed-composition CRC fuels). That null is **overturned
  for bulk-descriptor models only if the criterion above holds for freezing** with ridge or hgb. Otherwise it stands.
* Correlation models are judged by the same criterion (flash only).

## 4. Uncertainty
Split-conformal two-sided margin (`src/fuelblend/uq/conformal.py::two_sided_margin`, alpha = 0.10) from **dev LOGO out-of-fold** residuals; the conformal group is the
deduplicated fuel (assumes fuel exchangeability, violated across sources). Margins are in the target's raw units and are fixed before the locked run. Coverage on the locked set is
reported with `conformal.coverage` (counts) and as a fraction; target 0.90.

## 5. Order
1. `splits` (done) -> 2. this README -> 3. `dev` (LOGO, within-group, time split, coarse-group sensitivity; assumed-T viscosity) -> 4. `locked`, exactly once -> 5. `report`.
