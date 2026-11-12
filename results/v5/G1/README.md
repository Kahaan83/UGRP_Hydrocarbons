# G1: zero-fit physics on AJFTD fuels with measured GCxGC composition (v5)

Status: READINGS DECLARED 2026-11-11, before any score was computed. Results are written by
`scripts/v5/eval_ajftd_gcxgc_v5.py` to `RESULTS.md` in this folder; this file is not edited after scoring (add a dated note).

## Question

With a measured class x carbon-number composition (AFRL GCxGC, `data/real_fuels/v5/ajftd_gcxgc_composition.csv`), how well do
the unfitted physics baselines predict the fuel's measured flash point, kinematic viscosity, density and freezing point?

## Data fixed in advance

* Compositions: the 113 fuels with leaf sum 98 to 102 wt % (`ajftd_gcxgc_fuels.csv`, `composition_complete`).
* Targets: `data/real_fuels/v5/ajftd_fuels_long.csv`. One value per fuel, property and temperature: the median of the tests.
  Viscosity rows with `T_source = assumed_spec_-20C` are not scored.
* Bins are mapped to one representative compound each with the r3 table (`fuelblend.data.reconstruction_r3.SubReps`):
  n-paraffin, iso-paraffin, monocycloparaffin (cyclo_mono), alkylbenzene and alkene use the nearest carbon number within 3;
  dicycloparaffin (cyclo_di), cycloaromatic (indan_tetralin) and diaromatic (aromatic_di) use the nearest available carbon
  number. Tricycloparaffins and dienes have no representative and count as unmapped. The mapping is written to
  `data/real_fuels/v5/ajftd_gcxgc_mapped.csv` before scoring. Mapped mass is renormalised to 1.
* A fuel is eligible if its unmapped mass is at most 10 %. Ineligible fuels are listed, not scored.

## Models (nothing fitted on these fuels)

* Flash: ideal Le Chatelier and UNIFAC (`fuelblend.baselines.flash_unifac.flash_point_unifac_detail`), mole fractions.
* Viscosity: Refutas on the pure-compound provider (`scripts/v4/build_viscosity_leaderboard_v4.refutas_nu`).
* Density: ideal volume mixing from the same provider's pure densities, at the measured temperature.
* Freezing: ideal Schroeder-van Laar liquidus on mole fractions with the library's pure Tm and dHfus
  (`fuelblend.thermo.coutinho_wax.get_pure_thermal_properties`); comparator: leave-one-fuel-out mean of the measured values.

## Readings (fixed before scoring)

R1. Error is predicted minus measured. For each property: n fuels, mean bias, MAE; for flash also the fraction over-predicted
    by more than 5 C and the worst over-prediction. Viscosity is scored as ln(pred/meas), reported per temperature.
R2. Flash, safety reading (same thresholds as S2 R-C2): "measured composition removes the unsafe bias" holds if, for the ideal
    model, mean bias is at most +2.0 C AND the fraction over-predicted by more than 5 C is at most 0.20.
R3. Flash, UNIFAC: reported next to ideal; "UNIFAC helps" holds only if its MAE is lower AND its over-by-5 C fraction is not higher.
R4. Freezing: "the model has skill" holds only if its MAE is below the leave-one-out constant-mean MAE on the same fuels.
R5. Statistics are reported on three sets: all eligible fuels; fuels not shared with another source in this repo
    (labels 10325, 10289, 4658, 11769 are shared with FuelLib, CRC AV-27-18 or v4); and by the S4 split role
    (`results/v5/S4/splits_ajftd_v5.json`). Nothing is tuned, so the locked slice is a read, not a fit; it is still shown apart.
R6. Flash methods D93 and D56 are kept in a column and reported together and apart if both have at least 5 fuels.
R7. Fuels with identical mapped compositions count once in the fuel count note (the number of distinct compositions is reported).

No margin is written to `configs/`. A reading that cannot be evaluated is reported as "not evaluable" with the reason.


## Note, 2026-11-12: freezing without censored values

G1's freezing figure (n = 31, MAE 15.2 C vs leave-one-out constant 7.30 C) used 64 AJFTD freezing rows that carry a "<" qualifier as if they were measurements. Nothing above or in the G1 outputs was changed. Recomputed with uncensored values only (G1 predictions reused, `measured` replaced by the median of uncensored rows; six fuels with only censored values drop out): n = 25, bias +7.4 C, MAE 8.99 C vs leave-one-out constant 5.59 C, so "the model has skill" still does not hold. Subsets: not shared n = 21, MAE 9.25 (constant 6.18); dev n = 17, 8.73 (6.02); locked n = 8, 9.56 (5.29). File: `results/v5/S3/corrected_2026-11-12_g1_freezing_uncensored.csv` (script `scripts/v5/g1_freezing_uncensored.py`).
