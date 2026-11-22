# Methodology: current model set

This page states the method used for each property. The source of numbers is `docs/final_numbers_v5.md` (blends) and `docs/final_numbers_v7.md` (pure compounds). `Fuel_Blend_Methodology.pdf` has no committed source; this page is the text to carry into its next revision.

## 1. Two scopes

Blend properties (viscosity, density, flash point, freezing point, heat, smoke point, distillation) come from `predict_fuel`. Pure-compound properties (flash point, cetane, derived cetane, density, heat of fusion, log viscosity at 293 K, MON, RON, melting point, YSI) come from `predict_pure`. The method for each property is the one with the lowest committed cross-validated error among the deployable configurations, unless a predeclared gate or a paired bootstrap says the gain is not real, in which case the simpler designated route stays.

## 2. Blend properties

* Viscosity: Refutas blending with an optional LightGBM pair residual, gated on locked systems.
* Density: ideal volume mixing with an optional LightGBM excess-volume residual.
* Flash point: Le Chatelier rule with Raoult's law, anchored Antoine constants and an LFL overlay. Ideal activity coefficients. HANNA and UNIFAC activity coefficients showed no skill against ideal on hydrocarbon binaries; Antoine constants predicted by a D-MPNN head over-predict flash by 16 degC and are not used.
* Freezing point: Schroeder-van Laar with UNIFAC. Real-fuel freezing: status: no skill shown. Learned residual models were tied with simpler ones on development rows and did not reopen that verdict.
* Uncertainty: one-sided 90 percent conformal margins from held-out lab systems. Row-level, family-conditional and shift-weighted margins were studied (coverage by family) and remain exploratory.

No neural mixture model (set transformer, deep sets, MolPool) passed its gate against the v5 route, so none is deployed for blends.

Exploratory blend values from the pure-compound models: derived cetane and RON (and a MON entry marked no skill) are returned by `predict_fuel` as extra keys. The rule is a linear mole-fraction blend of the best pure-compound predictions, with no mixture fitting. Predeclared test on development mixtures: derived cetane MAE 6.90 against 21.57 for a constant (31 systems, 251 rows), RON 6.38 against 13.56 (5 systems), MON usable on KAUST surrogates (5 systems) but no skill on Leenhouts mixtures (32 systems), so MON is not offered. These values are uncertified: no conformal margin, no lab validation, and the RON evidence is thin. Existing keys and margins are unchanged. Cetane has no mixture data and is not blended.

## 3. Pure compounds

Training data: the Track B pool, one row per compound and property, with the locked YSI hold-out labels masked. Evaluation: scaffold-grouped 5-fold cross-validation on those development rows; errors are optimistic because the same rows guided tuning.

| property | deployed method |
|---|---|
| flash point, cetane, heat of fusion, log viscosity at 293 K, melting point, YSI | random forest on RDKit 2D descriptors, 500 trees, median imputation |
| derived cetane, density, MON, RON | multi-task D-MPNN, five-model ensemble, 1000-epoch cap with patience 100 |

Where two methods were tied (cetane, melting point, YSI) the simpler random forest is kept. The D-MPNN wins on derived cetane, MON and RON are nominal (paired-bootstrap interval excludes 0 before correction, Holm p above 0.10). The random forest wins on flash point, heat of fusion and log viscosity are significant after Holm correction. Pure-compound predictions have no conformal margin.

## 4. What was tried and not adopted

Joint pure plus fuel training for freezing (tied with fuel-only), equation-parameter heads for Antoine and Andrade constants (beat a pooled curve, not a descriptor ridge), HANNA activity coefficients in the flash route, set-transformer and deep-sets mixture encoders, oxygenate-augmented flash training (hurts hydrocarbon rows). Each is recorded with its numbers in `docs/v7_errors_by_property.md` and `docs/v6_negative_results.md`.

## 5. Rules that bound every claim

Predeclare before running; dev-only gains cannot reopen a v5 or v6 verdict; locked sets are scored once; paired group bootstrap over systems with Holm correction; no laboratory validation exists, so certification is model-level only.

## 6. v8: ML-filled constants (negative result)

Question: can out-of-fold pure-compound ML values replace measured constants in the v5 physics blend routes, and can a real-fuel margin be earned for represented blendstocks? Gates and tolerances were fixed in `docs/v8_predeclaration.md` before any result. Substituting ML density, viscosity, flash-derived LFL, and melting point with heat of fusion failed the dev gate in all four cases (dMAE upper bounds 11.8 kg/m3, 0.144 ln, 2.07 C, 4.81 C against tolerances 2, 0.02, 2, 2). The designated freezing route on real fuels has no skill over a constant (MAE 8.55 against 6.24 C, bias +5.66 C), so no freezing margin is earned. Designs were not re-run; the certified count stays at 3 of 16 and design certification remains model-level only. Numbers: `docs/final_numbers_v8.md`.

## 7. v9 data

The v9 round added data from sources whose values could be checked against the original files. Boehm 2022 neat real-fuel viscosities (checked cell by cell against the workbook) and ThermoML 2020 files not in v5 (checked against the archive) entered `data/unified_v9/`. Edwards 2017 and the pure-compound drafts could not be verified (the DTIC and ROSAP links returned HTML pages) and were withheld. All seven Antoine refit candidates failed for lack of verified source data.

The lock rule is md5 of the system id modulo 4 equal to 0, applied only to new v9 systems: 4 of 14 were locked, once. No freezing fuel was added and no composition-bearing viscosity system was added, so the 8-fuel floors were not met. Real-fuel freezing point: status: no skill shown. Cold viscosity: status: no skill shown. Numbers: `docs/final_numbers_v9.md`.
