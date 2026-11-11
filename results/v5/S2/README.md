# S2 flash: mesitylene input, Antoine range audit, sign-flip test (v5)

Status: SCORED 2026-11-11 (readings declared before scoring; see Results). Date declared: 2026-11-11. Everything below "Results" is filled in only after the scripts run; the readings must not be edited afterwards (add a dated note instead).

## Questions

A. Is the mesitylene vapour pressure used at the predicted flash temperature (about 317-323 K) an extrapolation, and does correcting it change the lab flash margin?
B. Which other pure compounds are used outside the valid range of their Antoine set at the flash temperature actually predicted?
C. Does measured GCxGC composition remove the unsafe (over-predicting) flash bias seen with reconstructed compositions, and does UNIFAC change its sign on aromatic-rich fuels?

## Data fixed in advance

* A: `data/antoine_sources.csv` (stored set, valid 354.64-438.87 K), the anchored set in `results/v4/F/pure_inputs_used.csv`, and any measured mesitylene vapour pressure in the ThermoML-derived tables on disk. `data/processed/v2/thermoml_clean_hc.csv` is known to hold 254 bad vapour-pressure rows: drop rows with phase or pressure flags and list what was dropped. No download.
* B: records in `results/v4/F/flash_lab_loso.csv` and the CRC 647 flash rows scored in v4/r3; predicted flash T per record taken from existing outputs, not re-fitted.
* C: T2 fuels with BOTH measured GCxGC composition and a measured flash point: 4 FuelLib fuels (`data/real_fuels/v5/fuellib_*`) and the 38 CRC AV-33-22 fuels (`data/real_fuels/v5/crc_av3322_*`, supplier CoA flash). Class x carbon-number bins are mapped to representative compounds using the r3 representative table (`results/v4/R3/reconstruction_report.md`); the mapping is written to disk before scoring. Unmapped mass is reported per fuel; fuels with more than 10 percent unmapped mass are excluded and listed.

## Readings (fixed before scoring)

R-A1. Extrapolation is a fact test: the mesitylene set is called extrapolated if the predicted flash T of any lab record using it lies outside 354.64-438.87 K. No threshold.

R-A2. If at least 6 usable measured mesitylene Psat points at or below 345 K exist, refit (same logic as `refit_antoine.py`) and compare Psat at 320 K: stored, anchored, refit. If fewer exist, no refit is made and the result is stated as "no independent low-T measurement on disk" (D7 stays open).

R-A3. "The flash margin no longer hangs on one compound" holds only if, with the best available mesitylene input, the 90 percent system-level margin (ideal and UNIFAC) is finite and differs by at most 1.0 C from the margin with all mesitylene records excluded. Otherwise it still hangs on mesitylene.

R-B1. The audit flags a compound when the predicted flash T of at least one record using it lies outside the stored Antoine validity range. It reports the flagged count, the records and systems affected, and whether each flagged compound has an independent Psat check. The count is the result; there is no pass/fail.

R-C1. Metrics, for ideal Le Chatelier and for UNIFAC, error = predicted minus measured: mean bias, MAE, and the fraction over-predicted by more than 5 C, on all eligible T2 fuels and separately on FuelLib-4 and AV-33-22.

R-C2. "Measured composition removes the unsafe bias" is read as: mean bias at most +2.0 C AND fraction over-predicted by more than 5 C at most 0.20 on the AV-33-22 fuels (ideal model). Reference numbers (different fuels, so a contrast, not a paired test): r3 dev bias +0.44 C and 13 percent over by more than 5 C; r2 +4.22 C and 53 percent. Where a T2 fuel is also in the CRC 647 set (see `v4_overlap`), report the paired difference on those fuels too; with fewer than 3 such fuels, say it cannot be paired.

R-C3. "UNIFAC changes the sign on aromatic-rich fuels": split eligible fuels at the median measured aromatic mass fraction; the claim holds if the mean error of ideal and of UNIFAC have opposite signs in the upper half. Report both numbers either way.

R-C4. Residual check: Pearson correlation of error with D86 T10 (FuelLib has it, AV-33-22 does not). Report n.

Nothing is tuned. No margin is written to `configs/`. A reading that cannot be evaluated is reported as "not evaluable" with the reason, not replaced by another reading.

## Results

Scored 2026-11-11 with `scripts/v5/s2_mesitylene.py`, `scripts/v5/audit_antoine_range_v5.py`, `scripts/v5/s2_t2_flash.py` (tests: `tests/test_s2_flash.py`). All numbers below are copied from the CSVs in this folder. The readings above were not edited.

### Question A: mesitylene

Inputs found. Stored WebBook set (Forziati 1949, `data/antoine_sources.csv`, valid 354.64-438.87 K): A=6.1993, B=1569.62, C=-63.57 (read as log10(P/kPa)). Registry ("anchored") set actually used by the flash model (`compounds_anchored_vp.csv`, source label `webbook`): A=6.14818, B=1632.74, C=-61.27. Measured mesitylene vapour pressure on disk: ONE point (101.33 kPa at 437.83 K, doi 10.1016/j.fluid.2015.09.030), present in `thermoml_hc.csv`, `v2/thermoml_clean_hc.csv`, the v5 re-extraction and `data/unified/measurements.csv` (4 rows, 1 distinct point). Nothing was dropped (phase Liquid, pure, no pressure flag, present in the fixed re-extraction); the known-bad class (rows in the old/clean tables absent from the fixed re-extraction) is 254 rows in `thermoml_hc.csv` and 257 in `thermoml_clean_hc.csv`, none of them mesitylene (`thermoml_bad_vp_rows_summary.csv`, `mesitylene_measured_psat.csv`).

Psat (kPa, `mesitylene_psat_compare.csv`):

| T (K) | stored WebBook | registry (used in v4) | `thermo` library ref (not a measurement) | measured |
|---|---|---|---|---|
| 310.95 | 0.715 | 0.406 | 0.782 | |
| 320.00 | 1.197 | 0.688 | 1.269 | |
| 437.83 | 101.24 | 64.90 | 101.17 | 101.33 |

Registry/stored = 0.568 at 310.95 K, 0.575 at 320 K, 0.641 at 437.83 K. At the one measured point the stored WebBook set is -0.09 percent off, the registry set is -35.9 percent off. Best available input by the fixed rule (smallest log error at the measured point, `mesitylene_best_input.json`): the stored WebBook set (mean |ln error| 0.0008 vs 0.445 for the registry set).

* R-A1: MET (extrapolated). All 9 mesitylene lab records have predicted flash T below 354.64 K: 312.98-323.70 K (ideal) and 312.74-322.79 K (UNIFAC) with the registry set; 312.36-315.34 K with the WebBook set. 9 of 9 outside the 354.64-438.87 K range, by roughly 31-42 K.
* R-A2: NOT MET as a refit trigger: 0 usable measured points at or below 345 K (need 6), so no refit was made. Statement: no independent low-T measurement on disk (D7 stays open). The 320 K comparison is therefore stored vs registry vs the `thermo` library reference only.
* R-A3: HOLDS. 90 percent system-level margin, lab set (`mesitylene_flash_sensitivity.csv`, `mesitylene_readings.csv`):

| mesitylene input | model | all 31 systems | excluding mesitylene (30 systems) | abs diff |
|---|---|---|---|---|
| registry (v4 as run) | ideal | 3.699 C | 0.311 C | 3.39 C |
| registry (v4 as run) | UNIFAC | 1.558 C | 0.314 C | 1.24 C |
| stored WebBook (best) | ideal | 0.311 C | 0.311 C | 0.00 C |
| stored WebBook (best) | UNIFAC | 0.314 C | 0.314 C | 0.00 C |

Both are finite and within 1.0 C, so with the best available input the margin no longer hangs on mesitylene. The v4 dependence was caused by the registry constants: with the WebBook set the 9 mesitylene records change from bias -0.69 C (ideal) to -4.74 C (UNIFAC -1.64 to -5.64 C), i.e. the former worst-system over-prediction (max +3.70 C) becomes an under-prediction (max -3.83 C). Whole-lab-set ideal / UNIFAC MAE and bias move from 3.74 / 4.20 C and -3.30 / -3.91 C to 3.89 / 4.41 C and -3.59 / -4.20 C. The registry rerun reproduced the stored v4 predictions to 1e-14 C before the substitution.

### Question B: Antoine range

* R-B1: 11 compounds flagged out of 79 used (79 compounds in 126 lab records + 58 fuel records, `antoine_range_flags.csv`, `antoine_range_flags_records.csv`). 128 of 184 records use at least one flagged compound; the flagged compounds touch 68 of 89 systems (31 lab systems + 58 fuels). Range basis: 49 compounds use the 273.15-373.15 K anchoring window (3 flagged), 16 use the refit data range (7 flagged), 11 have a WebBook range in `antoine_sources.csv` (1 flagged: mesitylene), 3 have no stored range (hexylbenzene, octylbenzene, dodecylbenzene) and are counted as range-unknown, not flagged.

| compound | records outside (lab / fuel) | systems | max extrapolation (K) | stored / `thermo` ref Psat at median flash T |
|---|---|---|---|---|
| tetralin | 58 (0 / 58) | 58 | 52.9 | 1.04 |
| mesitylene | 52 (9 / 43) | 44 | 41.7 | 0.54 |
| indane | 42 (0 / 42) | 42 | 24.0 | 0.92 |
| n-hexadecane | 20 | 3 | 51.0 | 1.04 |
| n-heptane | 18 | 4 | 5.1 | 1.00 |
| n-octane | 9 | 1 | 5.1 | 1.00 |
| n-propylbenzene | 9 | 1 | 38.4 | 0.95 |
| ethylbenzene | 9 | 1 | 36.3 | 0.98 |
| o-xylene | 9 | 1 | 28.9 | 0.97 |
| 1,2,4-trimethylbenzene | 2 | 1 | 0.7 | 0.97 |
| benz[a]anthracene | 1 | 1 | 3.7 | 2.73 |

Independent Psat check: all 11 have a `thermo`-library reference ratio at their median flash T (column above); only mesitylene (0.54) and benz[a]anthracene (2.73, a solid at those T, so the library reference is questionable) are off by more than 10 percent. ThermoML measured points within 5 K of the flagged flash T exist for 5 of 11 (hexadecane, heptane, o-xylene, 1,2,4-trimethylbenzene, benz[a]anthracene), but for the refit_measured compounds those points are the fit data, so they are not independent. Tetralin and indane (flagged in every CRC fuel record) have no measured point near the flash T.

### Question C: T2 flash

Mapping (`av3322_representative_map.csv`, `fuellib_representative_map.csv`, written before scoring): mapped mass 94.4-98.9 percent for AV-33-22 (mean 96.7) and 97.5-99.6 percent for FuelLib; no fuel exceeded 10 percent unmapped, so all 42 are included (38 AV-33-22 + 4 FuelLib). Unmapped mass is mostly tricycloalkanes (no realistic representative, mean 1.5 percent of AV-33-22 mass), diaromatics (no representative in the r3 table) and chain lengths beyond the 3-carbon gap limit. Primary scoring uses the registry as stored (mesitylene constants as in v4/r3), the same basis as the r2/r3 reference numbers; the corrected-mesitylene variant is in `t2_flash_metrics.csv` for contrast.

* R-C1 (error = predicted minus measured; `t2_flash_metrics.csv`, registry as stored):

| group | model | N | mean bias (C) | MAE (C) | frac over >5 C |
|---|---|---|---|---|---|
| all eligible | ideal | 42 | -0.83 | 1.88 | 0.024 |
| all eligible | UNIFAC | 42 | -1.13 | 2.03 | 0.024 |
| FuelLib-4 | ideal | 4 | +1.36 | 5.82 | 0.25 |
| FuelLib-4 | UNIFAC | 4 | +1.10 | 6.30 | 0.25 |
| AV-33-22 | ideal | 38 | -1.06 | 1.47 | 0.000 |
| AV-33-22 | UNIFAC | 38 | -1.36 | 1.58 | 0.000 |

The single fuel over-predicted by more than 5 C is FuelLib POSF-11498 (an ATJ, 99.6 percent iso-paraffin, no aromatics): +14.3 C ideal, +14.8 C UNIFAC against a measured 49.5 C. The worst AV-33-22 under-prediction is CRC_AV3322_28 (-10.7 C ideal). With the corrected mesitylene input AV-33-22 bias becomes -1.40 (ideal) / -1.73 (UNIFAC) C, still no fuel over by 5 C.
* R-C2: MET on the AV-33-22 fuels (ideal): mean bias -1.06 C (needs at most +2.0) and fraction over-predicted by more than 5 C is 0.000 (needs at most 0.20). Contrast with r3 dev (+0.44 C, 0.13) and r2 (+4.22 C, 0.53), different fuels. Direction caveat: the AV-33-22 error is on the under-predicting (safe-side) side, with a magnitude of about 1 C, not an over-prediction. Pairing with CRC 647: NOT EVALUABLE. 0 T2 fuels carry a CRC 647 `v4_overlap` (the only overlap is FuelLib POSF-10325 = the v4 GC fuel POSF-10325, not a CRC 647 fuel; AV-33-22 files have no overlap field); fewer than 3, so it cannot be paired (`t2_pairing.csv`). For information, POSF-10325 ideal error is -11.88 C in v4 r3 (class x carbon-number derived from D2887, measured 55.0 C) and -2.08 C with the FuelLib GCxGC composition (measured 48.0 C); the two use different measured flash values, so this is not a clean paired difference.
* R-C3: NOT MET. Median aromatic mass fraction of the 42 eligible fuels is 0.2008. Upper half (n = 21): mean error ideal -0.86 C, UNIFAC -1.26 C, same sign (negative); lower half (n = 21): -0.81 / -1.00 C.
* R-C4: FuelLib-4 (declared, n = 4): Pearson r(error, D86 T10) = -0.064 (ideal) and -0.079 (UNIFAC), p = 0.94 / 0.92, not interpretable at n = 4. The premise that AV-33-22 has no D86 was wrong: the file carries D86 T10 for all 38 fuels, so, as an undeclared extra, r = -0.51 (ideal, p = 0.0011) and -0.48 (UNIFAC, p = 0.0026), n = 38. Heavier-distilling fuels are more under-predicted. That extra is not one of the declared readings.

### Caveats

* CoA flash methods differ: 27 AV-33-22 fuels D56 (Tag closed cup), 10 IP 170 (Abel), 1 not specified; the FuelLib values come from different tables (Edwards 2017 Table A-3 for three fuels, NJFCP SwRI for POSF-11498). Bias is about -1 C in each AV-33-22 method group (`t2_flash_metrics.csv`), so method does not drive the result, but no method correction was applied. CoA values are supplier certificates, not repeat lab measurements.
* Compositions are lumped class x carbon-number bins mapped to one representative molecule each (the r3 representatives), so "measured composition" means measured bins, not measured molecules; 0.4-5.6 percent of mass per fuel (tricycloalkanes, diaromatics, long chains) is dropped and the rest renormalised. Mesitylene is the C9 alkylbenzene representative in every fuel (mole fraction 1-12 percent), so the unresolved mesitylene Psat issue touches these predictions; the corrected variant shifts AV-33-22 bias by about -0.3 C.
* n for FuelLib is 4 (one is an ATJ with no aromatics); the FuelLib-4 numbers and the R-C4 FuelLib correlation are dominated by single fuels and carry no statistical weight.
* Question A: the corrected margin rests on the WebBook set being extrapolated 31-42 K below its range and checked against a single measured point at 437.83 K plus the `thermo` library reference; it is not an independent low-T measurement. The mesitylene lab system pairs mesitylene with 2,2,4,6,6-pentamethylheptane, whose own constants are listed as a needs-removal candidate in `data/antoine_sources.csv` (CAS/SMILES mismatch), so that system's error is not attributable to mesitylene alone. The CRC 647 r3 flash was not rerun with the corrected mesitylene (42 of 57 reconstructions contain it at 0.05-1.4 percent mass); only the lab filter was.
* Question B: the validity range is taken from the stored source, which for 49 compounds is just the 273.15-373.15 K anchoring window (a data window, not a certified limit), so flags there mean "outside the window the constants were anchored on". Predicted T is the ideal-model flash T of the existing v4/r3 outputs; the UNIFAC T is in the records file.
* Nothing was written to `configs/`, no v4 output was changed, nothing was committed.
