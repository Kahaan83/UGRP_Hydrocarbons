# S3 composition error versus physics error (v5)

Status: SCORED 2026-11-11 (readings declared before scoring; see Results). Declared 2026-11-11. Do not edit the readings after scoring; add a dated note.

## Question

For real fuels with measured GCxGC composition, how much of the property error is composition guesswork and how much is physics? Each fuel is scored twice for density, kinematic viscosity, flash and freezing: (M) measured composition, mapped to representative compounds; (R) a composition reconstructed only from bulk descriptors available for that same fuel (D86 points, aromatics, and the r2/r3 fixed priors in `src/fuelblend/data/reconstruction_r2.py` and `reconstruction_r3.py`; where Table G.1 hydrocarbon types are not available for a fuel, the r2 route is used and labelled r2).

## Data fixed in advance

* Fuels with measured composition: FuelLib (4 with properties), CRC AV-33-22 (38), AJFTD GCxGC (113 complete, 39 with a measured target). Mapping: `map_fuel()` in `scripts/v5/eval_ajftd_gcxgc_v5.py`; AV-33-22 map from `results/v5/S2/av3322_representative_map.csv`.
* Fuels are grouped by double-count group BEFORE pooling: AJFTD labels 4658, 10289, 10325, 11769 overlap FuelLib/AV-27-18; overlap of AJFTD GCxGC with AV-33-22 is checked by label/sample and reported. Pooled statistics use one fuel per group (the one with the most targets); per-source tables are also given.
* The 11 AJFTD fuels with S4 role "locked" were already read in G1 for composition-based flash, viscosity, density and freezing. They are reported in a separate row and excluded from the pooled headline.
* Targets: density at 15 C; kinematic viscosity at -20 C and 40 C; flash (ideal Le Chatelier, registry as in v4 so numbers are comparable to r2/r3); freezing (ideal SLE). Property models are the existing ones (no refit).
* No parameter is tuned on any scored fuel. The R route has no access to the measured composition of the fuel it reconstructs.

## Readings (fixed before scoring)

R-S3-1 (central table): for each target, report mean bias and MAE for M and R on the SAME fuels, with a paired bootstrap over fuel groups (5000 resamples, seed 0) of the MAE difference (R minus M). "Composition error matters" for a target if the 95% interval of the difference excludes 0 in favour of M; "physics-limited" if the interval includes 0 AND M's |bias| is no smaller than R's. Report n per target; targets with fewer than 8 groups are reported as "n too small", not tested.

R-S3-2 (density): the v4/r3 density bias of -14 kg/m3 is attributed to the reconstruction if M's mean bias at 15 C has |bias| <= 7 kg/m3 AND R's |bias| >= 10 kg/m3 on the same fuels. Otherwise it is not attributed (state which side failed). G1 already found -6.1 kg/m3 for M (n=30, different property provider); this reading is about the paired R comparison.

R-S3-3 (flash safety side): report the fraction over-predicted by more than 5 C for M and R; unsafe bias is removed by measurement if M's fraction is at most half of R's and M's mean bias is at most +2.0 C.

R-S3-4 (unmapped mass): report mean and max unmapped mass fraction per source; fuels with more than 10% unmapped mass are excluded and listed.

R-S3-5 (cold viscosity): report M and R log-bias at -20 C (and -40 C where n>=8) for raw Refutas only; G1 gave -26% at -20 C for M. State whether R is worse than M.

Anything that cannot be evaluated is stated "not evaluable" with the reason.

## Results

Scored 2026-11-11 with `scripts/v5/s3_composition_vs_physics.py` (tests: `tests/test_s3_composition.py`). Files: `s3_central_table.csv`, `s3_readings.csv`, `s3_predictions.csv` (per fuel and target, M and R), `s3_fuels.csv`, `s3_r_compositions.csv`, `s3_groups.csv`, `s3_group_links.csv`, `s3_bulk_overlap_pairs_top60.csv`, `s3_av3322_vs_ajftd_composition_distance_top25.csv`. Errors are predicted minus measured (kg/m3, C; viscosity ln(pred/meas)).

### What was run

* M: AJFTD GCxGC via `map_fuel()` (G1 mapping, 112 of 113 fuels eligible); AV-33-22 via `results/v5/S2/av3322_representative_map.csv` (38); FuelLib via the S2 FuelLib map (4). Both M and R go through one function, `predict_comp()` (ideal Le Chatelier flash with the stored registry, Refutas, ideal volume density, ideal SvL freezing; the Tm/dHfus table is built over the union of M and R compounds).
* R: `reconstruction_r2.reconstruct` with the fixed r2 saturate prior. Inputs per fuel, nothing else: D86 IBP/T10/T20/T50/T90/FBP (median over tests), aromatics vol% (D1319 preferred, else D6379/CoA), naphthalenes vol% (AV-33-22 only; AJFTD has none, so 0), olefins (none in any file, so 0). Density is not an r2 input. No Table G.1 hydrocarbon types exist for these fuels (they are keyed to CRC 647 sample codes), so every R is r2, not r3. For AV-33-22 the inputs came from `crc_av3322_properties.csv` (D86 rows, supplier-CoA aromatics: 33 D6379, 3 D1319, 1 D8305, 1 unstated; D1840 naphthalenes); `avg_mw`, formula and hydrogen content (derived from the composition table) were not used. AJFTD inputs from `ajftd_fuels_long.csv` (D86, D1319/D6379).
* R was available for 33 of 113 AJFTD fuels (78 have no aromatics value, 76 of them because they have no row at all in `ajftd_fuels_long.csv`; 2 have no D86 points), 38 of 38 AV-33-22, and 0 of 4 FuelLib (their file has D86 but no aromatics value; I did not borrow aromatics from the GCxGC bins). All FuelLib fuels are therefore "not evaluable" for the paired comparison, and the paired set is AV-33-22 plus 33 AJFTD fuels, not the 155 fuels in the README.
* Double counts. The `double_count_group` columns are constant per source and were not used. Groups were built from: same POSF label (AJFTD 10289/10325 with FuelLib), the S4 `dup_map`, identical mapped-composition hash, and a bulk-property match fixed in advance (at least 3 of density 15 C, aromatics, D86 T10, T90, flash, freezing, all within 1.5 kg/m3, 0.5 vol%, 2, 2, 2, 2 C). Result: 147 groups from 155 fuels. Merged: AJFTD 2 (10325) + FuelLib 10325 (label) and AJFTD 122, 173 and 24568 (labels 12751, 5171, 13785) by bulk match to it (only 4 comparable features, so these three links could be coincidences); FuelLib 10289 + AJFTD 39 (label); AV-33-22 samples 01/12/35 and 08/23 (all six features within tolerance). AJFTD labels 4658 and 11769 have no partner in S3's sources, so they stay single. AV-33-22 versus AJFTD: no bulk-match link, and the smallest total-variation distance between any AV-33-22 and any AJFTD composition is 0.072 (7 % of mass), so no AV-33-22 fuel appears in AJFTD GCxGC. The headline uses one fuel per group (most paired rows, tie by source order AV-33-22, FuelLib, AJFTD) and drops every group containing one of the 11 S4-locked fuels (11 groups; 9 of those 11 fuels have a paired R). Headline = 56 groups (35 AV-33-22, 21 AJFTD). Sensitivity with strict groups (label, dup_map, identical composition only) gives 62 groups; every verdict below is unchanged.
* Prior: the r2 saturate split is the mean of four GC fuels including POSF-10325. One R fuel (label 10325) has a label in that set; sensitivity runs use a leave-out prior for it and, separately, drop its group. No number or verdict changed (the density R bias moved from -4.02 to -3.98 kg/m3).

### R-S3-1 central table (headline, 5000 paired group bootstraps, seed 0)

MAE difference is R minus M; positive means measured composition is better.

| target | groups | bias M | bias R | MAE M | MAE R | R-M (95 % CI) | verdict |
|---|---|---|---|---|---|---|---|
| density 15 C (kg/m3) | 54 | -3.3 | -4.0 | 4.4 | 11.4 | +7.0 (5.1, 9.0) | composition error matters |
| viscosity -40 C (ln) | 44 | -0.643 (-47 %) | -0.591 (-45 %) | 0.643 | 0.591 | -0.052 (-0.078, -0.020) | CI excludes 0 in favour of R; not in the declared labels (see below) |
| viscosity -20 C (ln) | 55 | -0.270 (-24 %) | -0.209 (-19 %) | 0.297 | 0.245 | -0.052 (-0.075, -0.024) | same as above |
| viscosity 40 C (ln) | 11 | +0.014 | +0.018 | 0.049 | 0.117 | +0.068 (0.021, 0.110) | composition error matters (n = 11, AJFTD only) |
| flash (C) | 55 | -0.36 | +11.6 | 2.6 | 11.7 | +9.1 (7.2, 10.8) | composition error matters |
| freezing (C) | 53 | +5.7 | +2.8 | 10.1 | 5.6 | -4.5 (-8.6, -1.6) | CI excludes 0 in favour of R; not in the declared labels (see below) |

Per source (all paired fuels, not deduplicated across sources): AV-33-22 flash MAE 1.5 (M) vs 13.0 (R, bias +13.0); AJFTD flash 3.7 vs 8.9 (R bias +8.4); AV-33-22 density bias -2.1 (M) vs -7.4 (R), AJFTD -6.1 vs +5.4 (the R density bias has opposite sign in the two sources); AJFTD viscosity -40 C and -20 C are "physics-limited" (CI includes 0, |bias M| >= |bias R|). The 9 S4-locked fuels with a paired R are reported apart in `s3_central_table.csv` (flash M bias -1.9 vs R +6.1, same direction as the headline; most targets n < 8, so "n too small").

The declared labels do not cover a CI that excludes 0 in favour of R. I report it as is: for low-temperature viscosity and freezing, R is not worse than M, and its smaller MAE is a reduction of an error both routes share (see R-S3-5), not evidence that reconstruction is more accurate.

### Verdicts

* R-S3-1: composition error matters for density, 40 C viscosity (n = 11 only) and flash. Not met as "physics-limited" for any target at headline level: the viscosity and freezing rows are not composition-limited but M is the worse or equal route there, so neither label applies. Physics-limited holds only for AJFTD alone (-40 C and -20 C viscosity) and the locked -20 C viscosity.
* R-S3-2: NOT attributed. M side holds (bias -3.3 kg/m3, limit 7). R side fails: bias -4.0 kg/m3, below the 10 kg/m3 needed, even though R's MAE is 2.6 times M's (11.4 vs 4.4). The reconstruction is noisy per fuel but unbiased in aggregate on these fuels, so it does not reproduce the v4/r3 -14 kg/m3 bias (different fuels, different route: r2 here). Same result under both prior sensitivities.
* R-S3-3: MET. Fraction over-predicted by more than 5 C: M 0.073 (4 of 55), R 0.891 (49 of 55); M mean bias -0.36 C (limit +2.0). R's over-prediction is large (+11.6 C) and one-sided, larger than the +4.2 C quoted for r2 on CRC 647 (that figure is on different fuels and was not reproduced here). Cause (not tested): the D86-derived carbon distribution has no C5-C8 light ends (for AV-33-22 sample 2 the lightest R compound is C9; flash 59 vs measured 43 C, M 43.4), and the ideal Le Chatelier sum is dominated by the lightest components. Unsafe bias is removed by measuring the composition, with the same flash model.
* R-S3-4: unmapped mass (M), mean / max: AJFTD 0.015 / 0.102 (113 fuels), AV-33-22 0.033 / 0.056 (38), FuelLib 0.014 / 0.025 (4). Excluded (> 10 %): AJFTD 13197 (0.1015). R fallback mass (r2): AJFTD 0.000 / 0.000 (33 fuels), AV-33-22 0.011 / 0.026; none above 10 %. R missing: 80 AJFTD and 4 FuelLib fuels (reasons in `s3_fuels.csv`).
* R-S3-5: raw Refutas is strongly low at cold temperature for both routes. -20 C: M -0.270 ln (-24 %; G1 gave -26 % on its set), R -0.209 (-19 %), 55 groups. -40 C (n = 44 >= 8): M -0.643 (-47 %), R -0.591 (-45 %). R is not worse than M; it is slightly better (MAE 0.245 vs 0.297 and 0.591 vs 0.643). The cold-viscosity failure is physics-limited: measuring the composition does not fix it.

### Caveats and surprises

* The paired set is small and selective: 56 headline groups, only those with D86 and aromatics in their own bulk data. AJFTD contributes 21 of them. FuelLib contributes none (R impossible without a non-GC aromatics value).
* "R minus M" mixes composition error with representative choice: M uses `SubReps` (r3 sub-kinds) and R uses `RealisticReps` (r2 kinds). This is inherent to the declared route and was not adjusted.
* Viscosity at 40 C exists only for AJFTD (11 headline groups); AV-33-22 has -40, -20 C viscosity only.
* The freezing M result (ideal SvL on the measured composition, bias +5.7 C, MAE 10.1) is worse than R's; a likely cause (not tested) is the ideal SvL liquidus being set by heavy n-paraffin tails. This is a physics-model limit, consistent with G1's finding that freezing has no skill, not a statement about reconstruction.
* The bulk-match links rest on tolerances fixed before the errors were seen; the AJFTD 122/173/24568 merge into the 10325 group rests on four features and may be a coincidence (strict-group sensitivity does not change any verdict).
* Mapped unmapped mass is renormalised away for both routes (diaromatic and tricycloparaffin classes in M, fallback lumps in R), so neither composition carries those classes.
