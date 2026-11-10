# v4 Task 17: design engine re-run (jet, diesel, gasoline; ASTM and BIS)

Engine `src/fuelblend/design/engine_v4.py`; physics wrappers `src/fuelblend/design/models_v4.py`; driver
`scripts/v4/run_designs_v4.py`; independent check `scripts/v4/verify_designs_v4.py`.  Specs come only from
`configs/specs_astm.yaml` / `specs_bis.yaml`; margins only from `configs/uncertainty_v4.yaml` via `applied_margin`.
All numbers below were recomputed from the CSVs in this folder (`designs_v4.csv`, `constraints_v4.csv`,
`robustness_v4.csv`, `comparison_v3_v4.csv`, `verification_v4.csv`).

## Bottom line

**0 of 13 designs pass.  No design is certified.**  (6 specs x 2 pools, plus one extra gasoline pool built to exercise the DVPE minimum.)
In v3, 6 of 6 designs reported `all_constraints_passed = True`; those properties were formula-derived (see the v3
fabricated-data note), so v3 never tested anything.  Re-scored by the v4 engine with the same recipes, none of the 6
v3 recipes passes either (comparison_v3_v4.csv: 0 of 6 `v3_recipe_rescored_passed`).

Status of each constraint over all 13 designs (77 constraint evaluations):

| constraint | pass / by_construction | fail | uncertified |
|---|---|---|---|
| flash_point_C | 5 | 1 | 2 |
| kin_visc_m20C_cSt (jet) | 4 | 0 | 0 |
| nhoc_MJ_kg (jet) | 4 | 0 | 0 |
| density_15C_kg_m3 | 5 | 4 | 4 |
| aromatics_wt_pct / vol_pct / ethanol_vol_pct | 6 | 1 | 1 |
| freezing_point_C | 0 | 2 | 2 |
| smoke_point_mm | 0 | 2 | 2 |
| kin_visc_40C_cSt (diesel) | 0 | 4 | 0 |
| dcn, cloud_point_C, pour_point_C | 0 | 0 | 8 |
| ron, mon | 0 | 2 | 8 |
| dvpe_kPa | 0 | 1 | 4 |
| distillation T50 / E70 | 0 | 0 | 5 |

## Which properties are certified (pass with a finite yaml margin) and which are not

Certifiable today (a pass can be reported, with the margin that was applied):

* **Flash point** (Le Chatelier, ideal): margin = max(lab 3.70 C, D93 0.0 C) = 3.70 C.  Verified against thermo/chemicals
  (all flash comparisons agree within 1 C, far inside the margin).  D56 (Jet A in some specs) has margin `inf`, so a
  D56-method pass is uncertified.
* **Density at 15 C**: margin 22.25 kg/m3 (two-sided), which is wider than most spec windows.  A pass is therefore
  only possible for gasoline (window 720-775, recipes at ~747).  Jet (775-840) fails or is uncertified with this margin.
* **Net heat of combustion**: 0.194 MJ/kg (one-sided, protects the minimum).
* **Kinematic viscosity at -20 C** (jet): recalibrated Refutas, ln-ratio margin 0.3465, valid -40..20 C.
* **Class mass/volume fractions and ethanol volume**: by construction (no model error).

NOT certified, and why:

* **Freezing point** (all jet designs): the only margin for an ideal solid-liquid model is `inf` (8 groups), so freezing is
  never certified.  The Coutinho-library alternate has a finite margin (4.05 C) but is validated on CRC real fuels only.
* **Smoke point**: AJFTD bulk-regression margin is `inf` (6 groups).
* **Diesel 40 C viscosity**: the margin is validated only at -40..20 C.  There are 0 diesel fuels in the validation
  corpus, so no diesel property has ever been validated.  Nothing is invented here: diesel designs are reported as
  uncertified/failing, never as validated.
* **DCN, cloud point, pour point**: excluded from the yaml (`excluded` section); no valid margin or no model.
* **RON / MON**: Task 15 validated octane models are base gasoline plus ONE blendstock, so multi-compound pools have no
  validated octane model.  The linear-by-volume values shown are nominal only.  Several pure RON/MON inputs are unverified (R16).
* **DVPE**: units of the DVPE data are unconfirmed and it is excluded from the yaml, so a DVPE pass is never certified.
  The **minimum is enforced**: a DVPE below the yaml minimum is `fail` (see design `gasoline_astm_d4814__v4pool_dvpe`,
  42.2 kPa against the 45 kPa minimum).  The BIS E10 gasoline spec carries only a maximum for DVPE in the yaml; I did not add an
  unverified minimum.
* **Distillation T50 / E70**: no model.
* **Pools using decalin_* or isopentane**: pure constants are `disputed`/`unverified`, so the model raises `ModelUnavailable`
  and the constraint is uncertified.  The v3 pool designs therefore collapse to degenerate recipes (the optimiser excludes
  the unavailable constraints; `n_excluded_from_optimiser` in designs_v4.csv records how many).

## What data would unlock certification

| property | needed |
|---|---|
| freezing, jet | an SLE validation set of real or measured multi-component blends (>= 20 independent groups) so the ideal SvL margin is finite |
| smoke point | >= 20 independent fuels with measured smoke point and composition |
| diesel anything | diesel fuels with composition and measured viscosity at 40 C, density, flash, cetane, cloud/pour |
| viscosity at 40 C | CRC-style fuels measured at 40 C, to extend the validated T range |
| RON/MON | measured multi-component blends (not base+one blendstock) |
| DVPE | confirmed units, then a conformal margin |
| density margin of 22 kg/m3 | a larger validation set; the current margin is dominated by few groups |
| D56 flash | more than one D56 flash group (n = 1 now) |

## Robustness (robustness_v4.csv, 108 rows)

Each constraint is evaluated with every model form (primary and alternates), each form judged with its own margin.  Forms
without a margin are uncertified; none borrows another form's margin.  Primary vs alternates disagree: for jet freezing the
Coutinho library passes on 2 designs while the ideal SvL and ideal-mole forms fail or are uncertified, so the freezing result
is not robust to model form.  Viscosity: the 7-feature ridge was moved to an alternate and now refuses to predict outside its CRC
feature envelope (it gave 0.13 cSt for pure decane when allowed to extrapolate).

## Independent verification (verification_v4.csv)

thermo/chemicals recompute Psat, Le Chatelier flash, ideal SvL liquidus and ideal-volume density (no fuelblend import,
enforced by `tests/test_v4_verify_imports.py`).  17 comparisons agree; 8 designs/properties were skipped (engine value not
available, so nothing to compare).  0 disagreements exceed the margin.  Where the engine margin is `inf`/absent (freezing) the
tolerance is a fixed 5 C stated in the script (`tolerance_is_engine_margin = False`); the freezing differences are 3.5-3.7 C.
This verifies the arithmetic only.  It says nothing about the validity of the ideal-mixture assumption.

## v3 vs v4, side by side (comparison_v3_v4.csv)

| design | v3 | v4 re-score of the same recipe | v4 redesign (same pool) | why it changed |
|---|---|---|---|---|
| jet ASTM / BIS | pass | not passed: 3-4 constraints uncertified (decalin density/dHfus disputed) | degenerate recipe (100% n-decane) with 3-4 constraints uncertified | v3 used formula-derived properties; v4 refuses unverified constants |
| diesel ASTM / BIS | pass | fail: 1.81 cSt at 40 C below the diesel minimum and out of the margin's validated T range; 4 uncertified | degenerate (100% n-hexadecane) | v3 viscosity at 40 C was not validated |
| gasoline ASTM | pass | not passed: RON/MON/DVPE/T50 uncertified (isopentane RON/MON unverified, isooctane Psat disputed) | still not passed | no validated octane/DVPE for this pool |
| gasoline BIS E10 | pass | not passed: same, only density and ethanol pass | still not passed | same |

With a pool of verified compounds ("v4pool"): jet fails density (22.25 kg/m3 margin), freezing (nominal -30.5 C and -39.3 C against
the limit) and smoke point; diesel fails viscosity, density and flash; gasoline passes density only.
