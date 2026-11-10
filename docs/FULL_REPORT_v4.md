# UGRP Hydrocarbons - Full Report v4 (integrity and upgrade)

Branch `v4-integrity`. Every number in this report is transcribed from a committed file under `results/v4/`, `configs/` or
`splits/v4/` and cited in brackets (path relative to the repo root). Nothing here is recomputed from memory.
Reproduce everything with `python scripts/v4/reproduce_all.py` (section 10). Claim-level status is in `docs/claims_register.md`.

## 0. Read this first

1. **The v3 corpus and design optimizer were formula-derived and are quarantined; v3 claims that depend on them are retracted.**
   The provenance audit [`results/v4/P1/AUDIT.md`, `results/v4/P1/provenance_ledger.csv`, 13,546 ledger rows] shows that v3 real-fuel
   properties and compositions were largely target-derived, formula-derived or placeholder constants (for example
   freezing = cloud - 3, pour = freezing - 4, diesel compositions back-solved from DCN and density, 804.0 / 795.0 kg/m3 jet densities).
   Of 153 v3 fuels, 92 are not in the v4 tables; the 61 kept fuels are all jet [`docs/claims_register.md` CLM-29,
   `data/quarantine/dropped_fuels.csv`, 92 rows]. Retracted: CLM-20 (octane lump model), CLM-21 (DVPE), CLM-22 (diesel cetane/cold flow),
   CLM-23 (v3 designs), CLM-24 (robustness check). Quarantined data: `data/quarantine/`. Pinned v3 tests moved to `tests/quarantine/`.
2. **What v4 can actually claim is narrow.** The real-fuel evidence is 57 CRC 647 jet fuels whose *compositions are reconstructed*
   (formula-derived from measured aromatics and D86; ruling R10), plus 210 gasoline blend rows (octane only) and a handful of GC jet fuels
   (only POSF-10325 is a headline fuel, n = 1, no statistics allowed). **Diesel has zero qualifying fuels.** Every error on CRC fuels mixes
   reconstruction error with physics error and cannot be separated [`results/v4/V/README.md`, `results/v4/B/README.md`].
3. **Result of the design re-run: 0 of 13 designs certified; the 6 v3 recipes re-scored: 0 of 6** [`results/v4/D/README.md`].
   This does not mean all fail physically (5 fail at the point estimate, 8 cannot be assessed).

## 1. Data and provenance mix

| corpus | N | provenance | source |
|---|---|---|---|
| v4 fuel inventory | 61 fuels (all jet: 57 CRC 647 reconstructed + 4) | none verified yet (`verified = False` for all 61) | `data/real_fuels/v4/inventory_v4.csv` (61 rows); `docs/claims_register.md` CLM-29 |
| v4 compositions | 2,835 kept rows; 518 composition rows quarantined (463 from compositions_v3.csv and 55 from compositions_measured.csv; exact-duplicate collapse only) | R15 gate: kept fraction >= 0.95 | `data/real_fuels/v4/compositions_v4.csv`; `data/quarantine/compositions_quarantined.csv` (518 rows) |
| CRC 647 targets | 57 fuels with `raw_file` flash, density, freezing, pour, kinematic viscosity (3 T), net HoC | targets raw_file; compositions `reconstructed` | `results/v4/P2/jet_corpus.md` (2,735 rows, 797 target rows) |
| GC jet fuels | 7 listed; 5 complete; 1 headline (POSF-10325, n = 1, D2887-split, partly formula-derived) | headline set allows no statistics | `results/v4/P2/jet_corpus.md` |
| Gasoline blends | 210 rows: 66 Co-Optima bRON (65 with RON, 66 with MON as targets), 144 FACE (never targets) | RON/MON `raw_file`; FACE compositions unresolved | `results/v4/P2/gasoline_corpus.md` |
| Lab mixtures (flash) | 126 records, 31 systems (123 binary / 30 systems; 3 ternary / 1 system) | measured, literature benchmark | `results/v4/F/README.md` |
| Lab mixtures (viscosity) | 14,237 dynamic-viscosity rows, 322 systems | measured | `results/v4/E/README.md` |
| ThermoML SLE | 778 scored points, 59 systems (57 binary, 2 ternary) | measured | `results/v4/SLE/README.md` |
| AJFTD smoke point | 100 fuels (849 in a broad sensitivity set) | measured, no molecule list | `results/v4/SP/README.md` |
| Pure constants | 411 rows audited: verified 216, disputed 78, unverified 65, not_comparable 26, reference_only 26 (counts of the status column in data/processed/v4/pure_props_v4.csv, 411 rows) | external-library cross-check | `data/processed/v4/pure_props_v4.csv`, `results/v4/P1/pure_constant_diffs.csv` |

Splits: 276 fuels in 210 groups; locked 56 fuels (42 groups), dev folds 41/40/45/48/46; the locked set has 12 CRC jet target fuels for each
CRC property and 12 RON / 12 MON fuels from only 4 blendstock CAS groups; **diesel has 0 qualifying fuels**
[`results/v4/P2/coverage_matrix_summary.md`, `splits/v4/real_fuels_v4_splits.json`]. The 12-fuel and 4-group locked sets give wide or
unreliable CIs; locked numbers below are point estimates unless a CI is shown.

## 2. Viscosity (kinematic, jet fuels, -40 to 20 C)

* Data: 57 CRC jet fuels x 3 temperatures = 171 observations, dev 45 fuels (135 obs) / locked 12 fuels (36 obs), compositions reconstructed
  [`results/v4/V/README.md`].
* Outer CV: leave-one-fuel-out over dev with inner grouped 5-fold CV for hyperparameters; locked fuels never used for fitting or selection.

Dev, all temperatures, log-RMSE [95 % fuel-cluster CI] [`results/v4/V/viscosity_leaderboard_v4.csv`, `results/v4/V/README.md`]:

| model | dev (45 fuels) | locked (12 fuels) |
|---|---|---|
| Raw Refutas | 0.4267 [0.3959, 0.4549] | 0.3806 [0.3135, 0.4564] |
| Recalibrated Refutas | 0.1351 [0.1157, 0.1550] | 0.1939 [0.1341, 0.2594] |
| Bilinear pairwise | 0.0901 [0.0799, 0.1014] | 0.1670 [0.0987, 0.2445] |
| Structured ridge | 0.0826 [0.0731, 0.0949] | 0.1497 [0.0810, 0.2289] |
| Structured GP | 0.0888 [0.0772, 0.1012] | 0.1517 [0.0860, 0.2296] |
| Lab-trained structured residual (zero-shot) | 0.3952 [0.3652, 0.4226] | 0.3523 [0.2857, 0.4268] |

* **What is supported:** bilinear, structured ridge and structured GP each beat recalibrated Refutas (paired fuel-cluster CI excludes 0 on
  dev overall and at -40 C; on locked overall and -40 C also, with 12 fuels) [`results/v4/V/README.md`].
* **Not claimed:** any ranking among bilinear / ridge / GP. Pairwise differences among them were not tested and their CIs overlap
  [`results/v4/V/README.md`]. The ridge has the lowest point estimate only.
* **GP under-covers:** 90 % interval coverage 0.770 dev / 0.694 locked overall (0.600 dev, 0.333 locked at -20 C); its intervals must not be
  used as calibrated [`results/v4/V/README.md`].
* **Lab zero-shot does not fix -40 C:** -40 C log-RMSE 0.6254 dev / 0.5442 locked, essentially raw Refutas (0.6735 / 0.5903) versus
  recalibrated 0.1788 / 0.1744 [`results/v4/V/README.md`].
* Hyperparameters sit at grid edges (GP length-scale 8 in 46 of 46 fits; bilinear rank 2 in 46 of 46) so grids may be too narrow
  [`results/v4/V/selected_hyperparameters_v4.csv`, `results/v4/V/README.md`].
* **Conformal margin:** structured ridge 0.162 ln-ratio (x1.176) over 45 groups, locked 10/12; recalibrated Refutas 0.347 (x1.414), locked
  10/12; at -40 C 0.150 and 0.347 [`configs/uncertainty_v4.yaml`; `results/v4/U/README.md`].
* **Set encoder (Task 16) lost to the ridge:** dev log-RMSE 0.0862 [0.0758, 0.0983] vs ridge 0.0826 (LOO) / 0.0833 (same folds); paired dRMSE
  +0.0036 [-0.0004, +0.0088], CI includes 0; not promoted. CPU-only, under-trained, 5 seeds x 40 epochs; no GPU result exists. Fuel-level
  conformal coverage 41/45 = 0.911 dev, **locked 8/12 = 0.67** (row coverage 0.970 is inflated) [`results/v4/E/README.md`,
  `results/v4/E/promotion_decision.json`].
* **Changes vs v3:** v3 R1 reported Bilinear LOO 0.1242 and recalibrated Refutas 0.1411 over all 57 fuels [`docs/claims_register.md` CLM-14,
  `docs/model_decision.md` section 2]; v4 trains dev folds on dev fuels only with v4 reconstructed compositions and is **not directly
  comparable** [`results/v4/V/README.md`]. v4 weakens "bilinear wins" to "three learned models beat recalibration, unranked".
  The bilinear init is now seed-deterministic (CLM-19a), which moved the bilinear v3 numbers (0.1157 vs 0.1223 committed v3) [`docs/claims_register.md`].
* **Still unverified:** compositions are reconstructed, not GC; learned residuals can absorb reconstruction bias; only 3 temperatures per fuel;
  nothing is validated for diesel/gasoline viscosity or for 40 C; whether the encoder could win with tuning or GPU training.

## 3. Flash point

* Data: 126 lab mixture records / 31 systems; 57 CRC 647 D93 fuels (dev 45 / locked 12; flash `raw_file`, composition reconstructed); D56:
  one fuel (POSF-10325, n = 1, no statistics); diesel and gasoline: 0 fuels with raw flash [`results/v4/F/README.md`].
* Models (all zero-parameter; nothing fitted): ideal Le Chatelier vs original-UNIFAC Le Chatelier [`results/v4/F/flash_summary.csv`].

| slice | model | N | MAE C | bias C | worst over C |
|---|---|---|---|---|---|
| lab, 31 systems | ideal | 126 | 3.74 | -3.30 | +10.60 |
| lab, 31 systems | UNIFAC | 126 | 4.20 | -3.91 | +11.10 |
| CRC D93 dev | ideal / UNIFAC | 45 | 9.44 / 10.05 | -9.44 / -10.05 | -2.53 / -2.89 (all under-predicted) |
| CRC D93 locked | ideal / UNIFAC | 12 | 8.44 / 9.03 | -8.44 / -9.03 | -5.68 / -6.17 |

* **Best model:** ideal. UNIFAC is worse: system-cluster paired bootstrap dMAE +0.45 C, 95 % CI [+0.11, +0.78] (UNIFAC better in 16 systems, worse in
  13, equal in 2) [`results/v4/F/README.md`, `results/v4/F/flash_lab_bootstrap.json`].
* **Margin:** lab ideal 3.699 C over 31 groups (k = 29); D93 reconstructed CRC computed -6.547 C over 45 groups, applied 0.00 C (negative because every
  CRC fuel is under-predicted; floor); D56 `inf` (1 group) [`configs/uncertainty_v4.yaml`, `results/v4/U/README.md`]. The D93 entry is an offset
  artefact, not a safety property; novel blends must use the larger of lab and D93 applied margins (never 0).
* **Suspect-input sensitivity (binding disclosure).** Two lab compounds have flagged anchored vapour pressures
  [`results/v4/F/README.md`]: mesitylene (anchored Psat 0.52x the external reference, 9 records, 1 system) and n-hexylbenzene (anchored Psat 1.20x the reference, 9 records, 1 system; its 'disputed' status in pure_props_v4.csv refers to the hand-listed registry value, not necessarily the Antoine set used here, per results/v4/F/README.md). Recomputed by `scripts/v4/flash_input_sensitivity.py` -> `results/v4/F/flash_input_sensitivity.csv`:

| subset (ideal model) | records | systems | MAE C | bias C | worst over C | 90 % margin C |
|---|---|---|---|---|---|---|
| all | 126 | 31 | 3.74 | -3.30 | 10.60 | 3.70 |
| excluding mesitylene records | 117 | 30 | 3.83 | -3.50 | 10.60 | **0.31** |
| excluding n-hexylbenzene records | 117 | 30 | 3.46 | -2.98 | 10.60 | 3.70 |
| excluding both | 108 | 29 | 3.53 | -3.17 | 10.60 | 0.31 |
| mesitylene records only | 9 | 1 | 2.66 | -0.69 | 3.70 | inf (1 group) |
| n-hexylbenzene records only | 9 | 1 | 7.44 | -7.44 | -5.97 | inf (1 group) |

  **Reading:** the headline lab margin of 3.70 C is *set by the single mesitylene system* (its maximum over-prediction is 3.699 C); without it the
  margin is 0.31 C. Because mesitylene's Psat input is itself flagged (0.52x reference), the lab margin is only as good as that input: if the
  reference Psat is right the model over-predicts that system because of a bad input; if the anchored value is right the margin stands. This
  analysis cannot tell which. n-hexylbenzene does not change the margin (its records are all under-predicted) but removing it lowers MAE from 3.74 to 3.46.
  UNIFAC gives the same margin pattern (1.558 -> 0.314) [same CSV]. The conservative 3.70 C is retained downstream.
* **Other unverified inputs:** LFL is Jones-estimated (0.55 x Cst) for 82 of 98 species; Antoine from `compounds_anchored_vp.csv` (63 species from
  unreviewed external correlations) [`results/v4/F/README.md`]. CRC compositions are reconstructed; closed-cup D93 vs model offset (about -9.4 C) is not separable from reconstruction error.
* **Changes vs v3:** CLM-12's 2.28 C figure belongs to the learned R2 row; the zero-shot ideal anchored baseline reproduces at 3.74 C [`results/v4/F/README.md`].
  CLM-26 D93 bias reproduces at -9.23 C on 57 fuels (v3 mixed set -8.12 C). The earlier "42.2 C / +4.2 C margin for D56/general" recommendation has no D56 support in v4 (n = 1).

## 4. Density (kg/m3, jet fuels)

* Data: 57 CRC fuels at 15.56 C (398 rows over all temperatures), ideal volume mixing with ThermoML pure densities (no constant fallback); compositions reconstructed
  [`results/v4/B/README.md`, `results/v4/B/density_summary_v4.csv`].
* Result: dev MAE 11.81 [9.65, 14.30], bias -11.11 [-13.98, -8.54] (N = 45); locked MAE 8.39 [5.37, 11.92], bias -7.93 (N = 12). The bias is one-sided and the CI
  excludes 0 in every slice [`results/v4/B/README.md`].
* Margin: 22.25 kg/m3 (two-sided, 15.56 C) and 23.61 (all T), 45 groups, locked 11/12 [`configs/uncertainty_v4.yaml`]. Replaces the v3 2.50 kg/m3 (CLM-25/28) that the v3 residuals could not support.
* v3 vs v4: v3 G10 MAE 16.86 kg/m3 on 82 jets with 77 placeholder densities; v4 11.09 on 57 raw fuels - **not a model improvement**, a target-quality change [`results/v4/B/README.md`].
* **Still unverified:** the **cause of the -11 kg/m3 bias is unknown** (candidates: reconstructed class split, tier-2 densities for 20 % of mass, dropped tail, deuterated stand-ins; the
  deuterium hypothesis is contradicted by r = +0.50) [`results/v4/B/README.md`]. No model beyond ideal volume mixing was evaluated in v4.

## 5. Net heat of combustion (MJ/kg)

* Data: 57 CRC fuels; group additivity [`results/v4/B/hoc_summary_v4.csv`].
* v2-comparable (constant 13.8 wt % H): dev MAE 0.191 [0.151, 0.230], locked 0.154 [0.099, 0.214], all 57 0.183 [0.148, 0.218] (bias -0.112). Composition-H variant 0.153 dev / 0.111 locked.
  Using *measured* H makes it worse (0.238 dev) [`results/v4/B/README.md`].
* Margin 0.194 MJ/kg (v2 variant) / 0.224 (composition-H), 45 groups, locked 12/12 [`configs/uncertainty_v4.yaml`]; replaces the hand-typed v3 0.28.
* v3 vs v4: v2 MAE 0.1792 reproduced at 0.183 on v4 inputs (CLM-18). The measured-H result suggests part of the apparent accuracy is error cancellation (not isolated).
* Unverified: Tier A measured formation enthalpies cover only 53 % of mass; the ridge residual (CLM-07) was not re-run in v4.

## 6. Freezing point and solid-liquid equilibrium

* **CRC fuels (57; dev 45 / locked 12):** the ideal WAT is biased warm by +18.04 C (dev) / +16.94 C (locked) so the legacy 14.58 offset is stale; the dev-fit offset gives locked MAE 3.61
  [1.90, 5.73], bias -1.10; Coutinho 4.25 locked / 5.57 dev; legacy offset 4.12 locked [`results/v4/B/README.md`, `results/v4/B/freezing_summary_v4.csv`].
  **No model shows skill over predicting the dev mean:** constant-mean MAE 3.32 C locked vs 3.61 for the offset model; 4.90 in-sample dev vs 4.40 [`results/v4/B/README.md`].
  Note the library Coutinho puts *mass* fractions into the Schroeder-van Laar log, and its parameters were calibrated on CRC/Dauphin data, so locked is not a clean holdout.
* **ThermoML SLE (new):** ideal SvL with verified constants has MAE 2.02 K (bias -1.10) on **only 8 systems / 103 points** (n-alkane blends, C6-C18); library Coutinho 9.97 K on 44 systems / 525 points;
  on the common 8 systems ideal 2.06 K vs Coutinho 3.85 K [`results/v4/SLE/README.md`, `results/v4/SLE/sle_summary.csv`]. Outside n-alkane blends Coutinho errs by tens of K.
* **Margins:** CRC offset model 10.115 C (45 groups), constant-mean comparator 6.914 C, Coutinho library 4.050 C; SLE Coutinho library 38.27 K over 44 systems (freezing_point_C.sle_alkane_blends.coutinho_lib_all); ideal-SvL margin is `inf` (8 groups < 9); indicative non-conformal point quantiles q90 4.05 K,
  q95 7.36 K with bootstrap 90th-percentile interval [+0.99, +7.40] K [`configs/uncertainty_v4.yaml`, `results/v4/SLE/README.md`].
* **Still unverified:** freezing models carry no demonstrated composition information on real fuels; for jet design the freezing alternates disagree (Coutinho vs ideal SvL) so the result is not robust to model form
  [`results/v4/D/README.md`]; no SLE validation set for aromatic or heavy systems.

## 7. Smoke point (new in v4)

* Data: AJFTD 100 fuels (smoke point 19.3-25.0 mm); 28 hand-typed pure-compound smoke points (23 uncensored) [`results/v4/SP/README.md`].
* Bulk route, grouped by inferred T-code (most conservative): selected-in-fold MAE 1.29 mm vs constant mean 1.47 mm, difference CI [-0.24, -0.11] (N = 100); grouped by fuel_id 0.90 vs 1.14
  (inflated by near-duplicate fuels); **broad set N = 849: 1.68 vs 1.65, CI [-0.09, +0.15], no skill** [`results/v4/SP/smoke_point_eval.csv`].
* YSI route (pure compounds, N = 23): MAE 5.07 vs 6.98; **fails on saturates** (five censored alkanes predicted 13-22 mm against > 50 mm). **Not evaluated on real fuels** (no fuel has both speciation and smoke point).
* Margin `inf` (6 T-code groups) [`configs/uncertainty_v4.yaml`]. D1655 limits are from secondary sources only; BIS spec values are UNVERIFIED; the `smoke_point_mm` yaml keys are config only and read by no code
  [`results/v4/SP/README.md`]. No v3 equivalent (previously unmodelled).

## 8. Octane (RON / MON), blends of base gasoline + one blendstock

* Data: Co-Optima bRON, 65 RON / 66 MON target rows, 21 blendstock CAS groups (dev 18 groups / locked 4 groups, 12 rows) [`results/v4/O/README.md`]. Only 2 distinct base gasolines (RON 90.3, 87.7), so base inputs are nearly
  uninformative. FACE rows (144) are never targets.
* Dev, stored 5 folds, RON MAE [95 % CI] (N = 53, 18 groups): base_unchanged 5.10, linear-by-volume prior 6.49, linear-by-mole 6.79, constant mean 3.22, pooled BOV 2.73 [2.02, 3.49], Ghosh lumps 3.19,
  GBM residual 2.55 [1.79, 3.38], GP residual 2.29 [1.61, 3.05]. MON (N = 54): base 3.26, constant 1.62, pooled 1.25, GBM 1.11, GP 1.15 [`results/v4/O/octane_eval_v4.csv`, `results/v4/O/README.md`].
* **Supported:** the learned models beat "base unchanged" (CI excludes 0) and the linear volume/mole priors are *worse* than leaving the base unchanged (dMAE +1.39 [0.80, 1.75] RON). pooled-BOV beats the constant mean on RON in both views (stored folds -0.49 [-0.81, -0.18]; leave-one-blendstock-out -0.56 [-0.93, -0.17]); GBM and GP beat it on stored folds only (-0.67 [-1.12, -0.19] / -0.93 [-1.40, -0.44]); in leave-one-blendstock-out their CIs include 0 (GBM -0.10 [-0.78, 0.69], GP -0.37 [-1.21, 0.58], Ghosh -0.27 [-0.80, 0.28]);
  **no CI-supported winner among pooled-BOV / Ghosh / GBM / GP** (their paired CIs include 0). The locked set (12 rows, 4 CAS groups) gives point estimates only (RON MAE Ghosh 1.48, pooled 1.84, GBM 1.61, GP 1.80).
* Margins are large: RON 6.4-9.4 ON for the learned models (0.83 for base_unchanged, 6.94 for the constant mean) and MON 4.0-8.4 ON for the learned models, all at rank 18 of 18 groups (the largest group over-prediction) [`configs/uncertainty_v4.yaml`, `results/v4/U/README.md`].
* v3 vs v4: v3 reported RON MAE 3.10 on 55 FACE gasolines (CLM-20); that was circular and is retracted. v4 numbers are for base + one blendstock only and are not comparable. Multi-component pools have no validated octane model [`results/v4/D/README.md`].
* Unverified: typed base assumptions (M = 100 g/mol, rho = 745 kg/m3), alcohol/ester/ketone lumps default to the tail prior, several pure RON/MON inputs unverified, no diesel/other fuels.
* **DVPE:** a diagnostic only; units are **UNCONFIRMED** (inferred psi for the NREL rows). On 51 BOB + oxygenate rows ideal Raoult beats the constant mean on dev (dMAE -0.26 [-0.48, -0.11]) but this does not confirm the unit; excluded from the margins yaml and never certified
  [`results/v4/O/README.md`].

## 9. Pour point, cloud point, cetane/DCN

* **No v4 model.** Pour point exists as a target on 57 CRC fuels [`results/v4/P2/jet_corpus.md`] but was not modelled; cloud/pour/DCN are in the `excluded` section of the margins yaml (no valid margin or no model).
  The v3 diesel cetane/cold-flow evaluation is retracted (CLM-22). Diesel evaluation is impossible until diesel fuels with composition and measured properties are added (`docs/v4_manual_downloads.md`).

## 10. Design engine and prediction API

* **Designs:** engine on library physics, specs only from `configs/specs_astm.yaml` / `specs_bis.yaml`, margins only via `applied_margin`. 13 designs: 5 fail at the point estimate, 0 meet-but-uncertifiable, **0 certified**, 8 undetermined;
  v3 recipes re-scored **0 of 6** [`results/v4/D/README.md`, `results/v4/D/designs_v4.csv`, `results/v4/D/comparison_v3_v4.csv`].
* Certification is jet-only (all margins validated on jet-range data); a diesel or gasoline pass using a jet margin is `pass_jet_derived_margin_only` and not counted.
* Freezing alternates disagree: Coutinho passes on 2 jet designs, ideal SvL does not. Independent verification (thermo/chemicals): 17 agree, 8 skipped, 0 disagree - arithmetic only [`results/v4/D/verification_v4.csv`].
* v3 designs and the "5-model robustness" check are retracted (CLM-23/24).
* **API:** `src/fuelblend/api/predict_blend.py` no longer retrains lazily; intervals use `applied_margin` only; only `fuel_type='jet'` can be in domain; viscosity/density/flash margins were validated on CRC jet fuels only; freezing, smoke, RON/MON, DCN, DVPE and pour are
  uncertified (value None, `in_domain=False`) [`docs/model_decision.md` section 5]. `fuel_type='jet'` is trusted from the caller.

## 11. Open items and unverified list

* ThermoML default-mode `phase_of()` reads only `ePhase`; real records use `eVarPhase` / `eConstraintPhase`, so vapour compositions of VLE data may leak into the default-mode extraction and `thermoml_hc.csv` may be contaminated. v4 density tier-1 pure inputs (data/processed/v2/thermoml_clean_hc.csv) come from the same default-mode extraction; pure rows are expected to be unaffected but this is unverified.
  **Recorded, not fixed** (out of v4 scope; SLE mode is unaffected) [`results/v4/SLE/README.md`].
* CRC 647 compositions are reconstructed (formula-derived); locked set is thin (12 CRC fuels, 4 octane CAS groups).
* Diesel: 0 qualifying fuels, nothing validated. Gasoline: octane only; no density/HoC/freezing/DVPE (units unconfirmed).
* D1655 smoke-point limits secondary-source-only; BIS spec values UNVERIFIED; IS 2796 DVPE minimum not citable (R4 not applied).
* Density bias cause unknown; viscosity ranking among bilinear/ridge/GP untested; GP under-covers (0.77 / 0.69); lab zero-shot does not fix -40 C; freezing models show no skill vs constant mean; Task 16 encoder lost (0.0862 vs 0.0833, CI includes 0), no GPU run.
* 26 `reference_only` pure constants and 63 anchored Antoine sets come from unreviewed external correlations [`results/v4/F/README.md`, Task 4 note].
* Not part of this report: v3 claims about a bilinear "winner", 2.28 C flash, 42.2 C flash threshold, density margin 2.50 - all revised in `docs/claims_register.md`.

## 12. Reproduction and determinism

`python scripts/v4/reproduce_all.py [--skip-tests] [--include-encoder] [--verify-only]` runs 19 build/eval steps in dependency order (2026-11-10: see 13.4, the pipeline now has 26 steps), pytest, row-count verification and writes `SHA256SUMS_v4.txt`. Steps guarded and reported as
skipped when not runnable: extras-dependent steps (when `extras/` is absent), encoder retraining (use `--include-encoder`), the one-off E-output docs patch, the Antoine/Cox comparison (separate uncommitted work), network downloads and the ThermoML SLE
re-extraction. Two consecutive full runs produced identical `SHA256SUMS_v4.txt` (see `.superpowers/sdd/2026-11-09-v4-integrity-and-upgrade/task-19-report.md`). The pytest baseline of pre-existing failures (all "deliverable file missing" asserts for v2/v3
files) is in `results/v4/baseline_pytest.txt`.

## 13. Review fixes (2026-11-10, branch `v4-review-fixes`)

Added after the v4 review; sections 0-12 above are unchanged. Claim rows: CLM-38, CLM-39, CLM-40; CLM-25 and CLM-33 carry dated notes.
All numbers are transcribed from `results/v4/R2/README.md`, `results/v4/L/README.md` and `results/v4/D2/README.md`.

**13.1 Realistic CRC 647 reconstruction (r2)** [`results/v4/R2/README.md`, `results/v4/R2/reconstruction_report.md`, CLM-38].
r2 replaces the v4 representatives with realistic ones that also have an Antoine and an LFL entry; the nominal saturates split
(n, iso, cyclo) = (0.1366, 0.3374, 0.5260) comes from 4 GC fuels. Mean fallback mass is 0.0202 (max 0.0538), 65% of it diaromatic mass
(no naphthalene representative); there is no bicyclic cycloalkane representative. Disclosure: the 12 locked fuels were already scored under the v4
reconstruction, so r2 choices were fixed on dev fuels and the GC-fuel prior only. Pre-declared readings on the 45 dev fuels:

| property | v4 | r2 | reading |
|---|---|---|---|
| density bias at 15.56 C (kg/m3) | -11.11 | -4.69 | outside +/-3, so not fully attributed to the old reconstruction (about 6.4 of 11.1 kg/m3 went away); 13 of 32 split draws within +/-3 |
| freezing, ideal SLE (mole) bias (C) | +14.69 | +0.32 | beating the dev-mean constant not established (MAE difference CI [-2.17, +0.25]); 20 of 32 draws within +/-3 |
| viscosity, raw Refutas dev bias (ln) | -0.301 | -0.244 | not below 0.10, so the recalibration was not just correcting the reconstruction |
| flash, raw ideal bias (C) | -9.44 | +4.22 | sign flips to the unsafe (over-prediction) side; 53% of dev fuels over by more than 5 C (locked +6.32); cause not isolated |

The flash flip is the safety-relevant result. The saturates split does not explain it (0 of 32 draws within +/-3 C of zero bias); representative
choice and lump handling were not varied. Dev-fit offset flash MAE rises from 1.71 to 4.12 C (dev) and 1.92 to 3.29 C (locked), against 5.51 / 5.90 for the dev-mean
constant. The Coutinho library and the legacy offset now under-predict by 12-17 C and should not be combined with r2. The compositions remain reconstructed.

Dated note 2026-11-10 (flash recommendation): the `configs/uncertainty_v4.yaml` FLASH_REC rationale ("the D93 reconstructed-CRC entry floors at 0 because every CRC fuel is under-predicted") no longer holds on r2 (raw bias +4.2 C; 53% of dev fuels over-predicted by more than 5 C). Designs are unaffected: they are known compositions judged by max(lab 3.70, D93) and D2 jet flash slack is about 9 C or more; the existing yaml key is unchanged.

**13.2 Known-composition lab margins** [`results/v4/L/README.md`, `configs/uncertainty_v4.yaml`, CLM-39]. Lab mixtures have measured compositions, so their residuals
isolate the physics from reconstruction error. Margins use dev systems only (84 of 419 systems locked). Density, ideal volume mixing: over 4.587 / under 3.250 kg/m3
(167 dev systems). Raw Refutas viscosity: T_ge_283K 0.186 ln (255 groups); T_lt_263K (7 groups) and T_263_283K (2 groups) are `inf` with no `applied_margin`.
Lab data are mostly binaries at 293-323 K; transfer to 100-component jet fuels is not demonstrated.

**13.3 Design re-run with lab margins** [`results/v4/D2/README.md`, CLM-40]. Category counts A/B/C/D are 7/0/0/6 against v4's 5/0/0/8; 0 of 13 are certified. All four jet
-20 C kinematic-viscosity constraints are uncertified, because the cold lab viscosity bin is `inf`. The gasoline E10 designs are now ethanol-free (the hydrocarbon-only density model
is unavailable with ethanol) and fail density. With `certified_kinds=None` the diesel 40 C kinematic-viscosity constraint goes fail -> pass and the gasoline ASTM density constraint goes pass_jet_derived_margin_only -> pass; these are margin-scope changes, not better designs (diesel is unvalidated). Pools, margins and rules were not adjusted.

Cold known-composition viscosity: below about 273 K, where the tier-1 pure densities do not cover the temperature, the known-composition API path is not Refutas; it falls to `gn_mass` / `gn_fallback`. At 253.15 K pure decane gives gn_mass 1.755 mPa.s against 1.966 mPa.s (gn_fallback) for decane/dodecane 0.999/0.001, a 12% jump, both `in_domain=False`, so nothing is certified from it (documented in the `predict` docstring).

**13.4 Open after these fixes.** Flash inputs (four `webbook` Antoine sets, Liaw formulation), the in-sample Coutinho margin, joint-coverage correction, measured-composition fuels and
low-temperature lab viscosity (needed to make the cold bin finite) are out of scope here. Reproduction: `scripts/v4/reproduce_all.py` now also runs the r2 reconstruction and
evaluations, the lab-residual build and the D2 design run, inserted in dependency order. The two default design steps (`run_designs_v4.py` / `verify_designs_v4.py` into `results/v4/D`) are no longer re-run: with the lab margins the default engine would overwrite the committed v4 D evidence (a first full run rewrote its five CSVs), so `results/v4/D` is frozen, verified and hashed, and `results/v4/D2` is regenerated instead. The pipeline now has 26 build/eval steps. `run_designs_v4.py` and `verify_designs_v4.py` now default to `results/v4/D2`; the original v4 D evidence is reproducible by checking out commit 65212cc (no legacy flag, because the direction-reason text format changed, so byte-for-byte reproduction in the current tree is not possible). Two consecutive full runs (`--skip-tests`) produced identical manifests apart from the five `results/v4/D` files that the first run had overwritten; a third run of the final script against itself was not made.
