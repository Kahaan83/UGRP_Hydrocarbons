# v4 Task 9 - flash point: ideal vs original-UNIFAC (zero-parameter physics)

Reproduce: `.venv/Scripts/python.exe scripts/v4/eval_flash_v4.py` (~30 s). Code: `src/fuelblend/baselines/flash_unifac.py`
(`flash_point_unifac`, `flash_point_ideal`, `flash_point_unifac_detail`), solver hook `gamma_fn` in `baselines/flash.py`.

Model: Le Chatelier H(T)=sum x_i*gamma_i*Psat_i/(P*LFL_i(T)) = 1, brentq (150-650 K), Zabetakis LFL correction, units kPa/K/degC.
gamma = **original (Hansen 1991) UNIFAC** from `thermo.unifac` (UFSG/UFIP, version=0), NOT Dortmund; hydrocarbon groups
assigned by an RDKit fragmenter (CH3/CH2/CH/C, alkene 5-8,70, ACH/AC, ACCH3/ACCH2/ACCH).
Fallback (documented, warned, reported): a molecule with no groups (alkyne, fluorene, ACC(quaternary), triphenylmethane, ...) gets gamma=1;
UNIFAC runs on the rest renormalised.

**Nothing is fitted.** Ideal and UNIFAC have no trained parameter, so "leave-one-system-out" = each system is scored by a model that never
saw it; no threshold/margin was chosen on any scored row. Ideal-model numbers reproduce `results/v3/R2` "Raw Zero-Shot Physical Raoult" exactly
(3.74 / 4.55 / -3.30 / +10.60); the 2.28 C CLM-12 figure is the other (learned) R2 row, not this baseline.

## Lab mixtures (126 records, 31 systems; 123 binary/30 systems, 3 ternary/1 system)
| model | N | MAE C | RMSE C | bias C | over>5 C | worst over C |
|---|---|---|---|---|---|---|
| ideal  | 126 | 3.74 | 4.55 | -3.30 | 1/126 (0.8 %) | +10.60 |
| UNIFAC | 126 | 4.20 | 4.95 | -3.91 | 1/126 (0.8 %) | +11.10 |

* The single >5 C over-prediction is the same record in both (system 39a914ff0579, measured 57.7 C; C12+ mixture).
* System-cluster paired bootstrap (5000): MAE(UNIFAC)-MAE(ideal) = +0.45 C, 95 % CI [+0.11, +0.78] -> UNIFAC is **worse**, small effect.
  Systems: UNIFAC better in 16, worse in 13, ~equal in 2 (|dMAE|<0.05 C).
* Aromatic-containing records (85 rec, 13 sys): ideal MAE 3.98 -> UNIFAC 4.83 (bias -3.58 -> -4.69; gamma up to 1.43 at the flash point lowers flash).
  Non-aromatic records (41 rec, 18 sys): 3.26 -> 2.89 (gamma slightly <1 raises flash).
  Reading: the gamma>1 direction is physically right for aromatic/alkane, but the ideal model already under-predicts; the residual error is not
  non-ideality (Antoine/LFL/anchoring dominate). Under-prediction is the safe direction; UNIFAC makes it more conservative, not more accurate.
* No lab mixture used a UNIFAC fallback (all 126 fully assignable).

## Real fuels (only fuels with a `raw_file` flash point in `jet_properties_v4.csv`)
CRC 647 (N=57, **ASTM D93**; composition provenance = `reconstructed`, i.e. formula-derived, NOT measured; 48 species, unresolved tail dropped
and renormalised; ~3.0 % mass (1.3-5.2 %) in UNIFAC-fallback species). Only the flash value is raw-file; the composition is derived.
| slice | model | N | MAE | bias | over>5 C | worst over |
|---|---|---|---|---|---|---|
| D93 dev folds | ideal / UNIFAC | 45 | 9.44 / 10.06 | -9.44 / -10.06 | 0 / 0 | -2.53 / -2.89 (all under-predicted) |
| D93 locked | ideal / UNIFAC | 12 | 8.44 / 9.03 | -8.44 / -9.03 | 0 / 0 | -5.68 / -6.17 |
| D93 all (57) | ideal / UNIFAC | 57 | 9.23 / 9.84 | -9.23 / -9.84 | 0 / 0 | -2.53 / -2.89 |
Locked fuels were scored only as a never-tuned slice (nothing was tuned). Ideal D93 bias -9.23 C matches CLM-26 (-9.23). ASTM D56: no CRC fuel has D56 -> N=0.
ASTM D56: one fuel only, POSF-10325 (55.0 C, av2718 p14): **n=1, headline_statistics_allowed=false: no statistics are reported for it.**
Its row in `flash_real_fuels.csv` is for traceability only (ideal 43.1, UNIFAC 42.9; molecules = registry representatives of class x carbon-number
pseudo-components; across 20 random representative draws the ideal prediction spans ~26-45 C, so representative choice dominates).
No diesel/gasoline fuel has a raw-file flash in v4 (SD-02, diesel nil; gasoline blends are octane-only).
Per-family N: jet/CRC647 57 (D93), jet GC 1 (D56, n=1), diesel 0, gasoline 0.

## Pure inputs used (`pure_inputs_used.csv`, 98 species)
* Antoine: `data/processed/v2/compounds_anchored_vp.csv` (CLM-12): 63 anchored_thermo_ambient (chemicals/thermo correlations, **reference_only-type, unreviewed**),
  31 refit_measured, 4 webbook. NOT the Cox-gate registry (uncommitted work in `results/v4/antoine_cox_gate` leaves this file unchanged).
* LFL: 16 measured (Zabetakis 627, OCR-read), 82 Jones 0.55*Cst estimates -> LFL is the largest unverified input.
* Cross-check vs `pure_props_v4.csv` (Psat at 37.8 C, status as in Task 4): 28 of 98 species have a v4 row (9 verified, 14 disputed, 3 unverified,
  2 reference_only). The "disputed" status concerns the hand-LISTED value in the old registry vs chemicals, not the Antoine used here; the anchored
  Psat is within 5 % of the chemicals reference for 26/28; outliers: mesitylene 0.52x, n-hexylbenzene 1.20x. 70 species have no v4 Psat row (no independent check).
  Pure flash status (v4): 8 verified, 12 reference_only, 5 disputed for the species that have one.
* Registry quirk: CRC647 reconstructed compositions include perdeuterated n-alkane SMILES ([2H]); grouped as ordinary alkanes (tested).

## Claims register rows these results would update (not edited)
CLM-12 (ideal anchored zero-shot reproduced: 3.74 MAE; note which R2 row the 2.28 is), CLM-26 (CRC647 D93 bias -9.23 C reproduced on 57; locked 12 / dev 45 split),
CLM-13 (D56/general vs D93 conclusion: only D93 N=57 is supportable; D56 has n=1), CLM-09/CLM-28 (flash margin: ideal worst lab over-prediction +10.6 C
on this exact baseline; margin recalibration must come from v4 data, UNIFAC does not reduce it).
