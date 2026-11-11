# S5 freezing with measured n-alkane distributions (v5)

Status: SCORED 2026-11-11 (readings declared before scoring; see Results). Declared 2026-11-11. Do not edit the readings after scoring; add a dated note.

## Question

Freezing in jet fuel is driven by the n-paraffin share, which only GCxGC measures. Does a first-species-to-freeze model (Bell, Boehm and Heyne 2025, doi 10.1021/acs.energyfuels.4c06091) built on the measured n-alkane distribution beat (a) a constant and (b) the ideal Schroeder-van Laar (SvL) model on the full measured composition? G1 found ideal SvL MAE 15.2 C against 7.3 C for a leave-one-out constant (n = 31, AJFTD).

## Data fixed in advance

* Fuels with measured GCxGC and a measured freezing point: AJFTD GCxGC fuels with freezing (31 scored in G1), CRC AV-33-22 (38, supplier CoA freezing; qualified values such as "<" are censored, excluded from the headline and listed), FuelLib (4 or fewer). Compositions from `data/real_fuels/v5/ajftd_gcxgc_composition.csv`, `crc_av3322_composition.csv`, `fuellib_gcxgc_compositions.csv`.
* Grouping and locking reuse `results/v5/S3/s3_groups.csv` logic: one fuel per double-count group; the 11 S4-locked AJFTD fuels are reported apart and excluded from the headline.
* n-alkane pure data: Tm and dHfus from the audited table used by `scripts/v5/eval_sle_extended_v5.py` / v4 SLE (no new constants; n-alkanes missing from the table are listed and their mass counted as non-solid-forming).
* Freezing-point test methods differ across sources (D2386 vs automated); stated as a caveat, not corrected.

## Models (all leave-one-group-out; no fitting on the scored group)

* C: mean of training freezing points.
* I: ideal SvL on the full mapped composition (as in G1).
* F0 (no fit): for each n-alkane species i present, the temperature T_i at which its ideal solubility ln x_i = -(dHfus_i/R)(1/T - 1/Tm_i) is reached with x_i its mole fraction in the whole liquid; predicted freezing = max_i T_i (first species to freeze). Non-n-alkane mass is an inert solvent.
* F1: linear map freezing = a + b * T_first, a and b fitted on the training groups only.
* P: ridge (alpha by inner CV) on [n-alkane mass fraction grouped C9-C12, C13-C16, C17+, total n-alkane share, aromatic share].

## Readings (fixed before scoring)

R-S5-1: F0, F1 and P each "beat the constant" only if the 95% interval (paired bootstrap over groups, 5000 resamples, seed 0) of MAE(C) minus MAE(model) excludes 0 on the positive side, with at least 20 groups.
R-S5-2: F0 and F1 versus I, same test.
R-S5-3: "the freezing null is turned" only if F1 or P beats the constant (R-S5-1) on AJFTD and on AV-33-22 separately, both with positive differences. One source alone is "source-specific", not a turned null.
R-S5-4: report bias, MAE and Spearman correlation of T_first against measured freezing, per source and headline, with n.
R-S5-5: locked AJFTD fuels reported apart; no claim from them (n < 12).

If the null stands, say so and put the T3 result (S4: freezing from bulk descriptors has no skill) beside it. A reading that cannot be evaluated is stated "not evaluable" with the reason.

## Results

Scored 2026-11-11 with `scripts/v5/s5_freezing_nalkane.py` (tests: `tests/test_s5_freezing.py`, 5 pass). Files: `s5_predictions.csv` (per fuel, all models), `s5_metrics.csv`, `s5_bootstrap.csv`, `s5_readings.csv`, `s5_excluded.csv`, `s5_nalkane_constants.csv`, `s5_fuels.csv`. Error = predicted minus measured, deg C. Bootstrap difference = MAE(reference) minus MAE(model), positive means the model is better; interval is the 95% percentile interval over 5000 group resamples, seed 0.

### Sets actually scored

* Headline: 50 groups (AJFTD 13, AV-33-22 34, FuelLib 3). Locked AJFTD (apart): 8 groups. One fuel per double-count group, S3 source order AV-33-22, FuelLib, AJFTD.
* Excluded (`s5_excluded.csv`, 16 rows): 8 fuels whose only freezing values are censored ("<"; 7 AJFTD, of which 1 locked, and AV-33-22 sample 02), 8 double-count duplicates. G1's n = 31 included censored values and locked/duplicate fuels, so its 15.2 C is not comparable; the G1 ideal-SvL predictions are reproduced exactly (max difference 1e-12 C), so model I here is the G1 model.
* AJFTD:325 (freezing -81 C) has no n-alkane in its GCxGC composition: F0, F1 and I are undefined for it (n = 49 for those, 50 for C and P). AJFTD freezing is D5972 throughout; AV-33-22 is supplier CoA (IP 529, D5972, D2386 mixed, sample 03 reports -80 C by D2386); FuelLib is a table value. Test-method differences are not corrected.

### Headline numbers (LOGO, deg C)

| model | n | bias | MAE | Spearman |
|---|---|---|---|---|
| C constant | 50 | 0.0 | 6.41 (5.94 on the 49 shared with F0/F1/I) | n/a |
| I ideal SvL | 49 | +7.99 | 7.99 | 0.86 |
| F0 | 49 | +7.00 | 7.95 | 0.86 |
| F1 | 49 | -0.20 | 3.46 | 0.85 |
| P ridge | 50 | +1.20 | 4.80 | 0.80 |

### Readings

* R-S5-1 (beats the constant, headline, n >= 20): F0 NOT MET (MAE 7.95 vs 5.94; difference -2.02, interval -4.05 to +0.32, n = 49). F1 MET (+2.47, +0.44 to +4.52, n = 49). P MET (+1.61, +0.30 to +2.87, n = 50). Per source: AV-33-22 (n = 34) F1 +2.93 (+1.64 to +4.39) and P +2.28 (+1.07 to +3.54) MET, F0 not met; AJFTD (n = 12 or 13) NOT EVALUABLE (F1 +1.79, -5.40 to +9.67; P +0.71, -2.82 to +4.18).
* R-S5-2 (versus ideal SvL I): F0 NOT MET on the headline (+0.03, -1.10 to +0.72, n = 49); on AV-33-22 alone F0 is better (+0.50, +0.18 to +0.87, n = 34), a small gain. F1 MET (+4.52, +2.68 to +5.81, n = 49; AV-33-22 +4.62, +3.64 to +5.57); AJFTD alone NOT EVALUABLE (n = 12). F1 is fitted (two parameters, leave-one-group-out) and I is not, so this gain is mostly removal of the +7 to +8 C bias that I and F0 share.
* R-S5-3 (null turned): NOT EVALUABLE. The rule needs F1 or P to beat the constant on AJFTD and on AV-33-22 separately; AJFTD has 12 to 13 groups, below the 20 required. AV-33-22 alone satisfies the condition for both F1 and P. The null is therefore not shown to be turned; the result is source-supported (AV-33-22) only, and the AJFTD point estimates are positive but their intervals are wide (n = 12 or 13). Sensitivity: swapping the group-source priority gives at most 15 AJFTD groups (still below 20); training within source (`C_ws`, `F1_ws`, `P_ws`) leaves AV-33-22 positive (F1 +3.49, +1.92 to +5.32; P +3.41, +1.98 to +4.85) and AJFTD worse than the constant in point estimate (n = 12 or 13, not evaluable).
* R-S5-4 (T_first = F0 against measured freezing): headline n = 49 bias +7.00, MAE 7.95, Spearman 0.86; AJFTD n = 12 bias +6.71, MAE 10.61, Spearman 0.87; AV-33-22 n = 34 bias +7.04, MAE 7.04, Spearman 0.86; FuelLib n = 3 bias +7.74, MAE 7.74, Spearman 0.5. T_first ranks fuels well but sits about 7 C too warm (ideal solubility, no dCp, no solid-solid terms), so unfitted it loses to a constant.
* R-S5-5 (locked AJFTD, n = 8, below 12, no claim): MAE C 5.18, I 9.56, F0 13.18, F1 8.85, P 4.91. F0 is driven by AJFTD:98, where 17.6% of the mass is n-alkane without audited constants and the prediction is -97 C (measured -56).

### Caveats and sensitivities

* Missing n-alkanes in the audited table: C8, C9, C10, C17, C19 (and C21 to C24); their mass is counted as non-solid-forming. F0 loses most where this matters (AJFTD:174, 39% of n-alkane mass without constants, F0 -77 C vs measured -54). The sensitivity F0x/F1x uses the repo's complete Coutinho table (not a declared model): F1x headline MAE 3.07 (n = 49), F1x vs C +2.87 (+1.05 to +5.16); F0x vs C -2.40 (-4.48 to -0.09), i.e. F0x is worse than a constant.
* High-leverage fuel: AJFTD:131 (51% n-alkane, freezing -8 C). Without it (refit, n = 48) F1 vs C is +1.59 (-0.34 to +2.94), no longer excluding 0; P vs C is +1.97 (+0.30 to +3.41) and still holds. So R-S5-1 for F1 depends on this one fuel; P does not. P's LOGO prediction for AJFTD:131 is +50 C (measured -8), a 58 C extrapolation error that inflates P's AJFTD RMSE (17.3).
* The fitted gains come from calibration across a narrow range (most jet fuels freeze between -62 and -40 C); a constant is already within about 6 C. Ridge alpha was 10 in every fold.
* Heavy n-alkanes in trace amount can set T_first (n-C20 at 0.004 mass% in several AV-33-22 fuels) because the ideal model has large dHfus and no solid-solution term.
* The audited dHfus values are fusion-only (e.g. C13 28.5 kJ/mol against C12 36.8 kJ/mol, odd-even); no solid-solid transition enthalpy.

### Verdict

The null is not turned under the declared rule (R-S5-3 not evaluable). With the measured n-alkane distribution, a one-parameter-pair recalibration of the first-to-freeze temperature (F1) and a ridge on n-alkane shares (P) beat a constant on the pooled headline and on AV-33-22 alone, by about 1.6 to 2.5 C MAE (F1 depends on one waxy fuel; P does not), while the unfitted physics (F0, I) do not. The T3 result from S4 (freezing from bulk descriptors has no skill over a constant) stands for bulk descriptors; this result says n-alkane-resolved composition has some skill on AV-33-22 and uncertain skill on AJFTD, and that the part that works is empirical calibration on top of a physics feature rather than the ideal model itself.
