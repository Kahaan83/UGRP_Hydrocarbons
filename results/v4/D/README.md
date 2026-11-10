# v4 Task 17: design engine re-run (jet, diesel, gasoline; ASTM and BIS)

Engine `src/fuelblend/design/engine_v4.py`; physics wrappers `src/fuelblend/design/models_v4.py`; driver
`scripts/v4/run_designs_v4.py`; independent check `scripts/v4/verify_designs_v4.py`.  Specs come only from
`configs/specs_astm.yaml` / `specs_bis.yaml`; margins only from `configs/uncertainty_v4.yaml` via `applied_margin`.
All numbers below were recomputed from the CSVs in this folder (`designs_v4.csv`, `constraints_v4.csv`,
`robustness_v4.csv`, `comparison_v3_v4.csv`, `verification_v4.csv`).

## Bottom line

13 designs (6 specs x 2 pools, plus one extra gasoline pool built to exercise the DVPE minimum).  `category` in
`designs_v4.csv` separates what the numbers can and cannot say:

| category | designs | meaning |
|---|---|---|
| (a) A_fails_at_point_estimate | 5 | at least one constraint is violated by the point estimate itself (margins ignored) |
| (b) B_meets_at_point_estimate_not_certifiable | 0 | every property has a point estimate and meets its limit, but a margin is inf or out of scope |
| (c) C_certified | **0** | every constraint pass/by_construction, no extrapolation |
| (d) D_undetermined | 8 | no violation among the properties that have a point estimate, but 1-4 properties per design have none (disputed or unverified constants, no model); cannot be called a pass or a fail |

So "0 certified" does NOT mean "all 13 physically fail": 5 fail at the point estimate and 8 cannot be assessed.
The v3 comparison (6 of 6 reported pass in v3, 0 of 6 now) must be read the same way.  v3 properties were formula-derived (see
the v3 fabricated-data note), so v3 never tested anything.  Re-scoring the v3 recipes with the v4 engine, they are blocked by
uncertified constraints, and for the v3-pool redesigns also by the exclusion of constraints that involve disputed or unverified
constants (decalin_*, isopentane) from the optimiser.  That is a data-availability cause, not evidence that the recipes
physically fail.  The exception is the two v3 diesel recipes: their nominal 40 C viscosity (1.81 cSt) is below the limit.

Constraint evaluations over the 13 designs (77 rows in constraints_v4.csv):

| constraint | pass / by_construction | pass_jet_derived_margin_only (not certified) | fail | uncertified |
|---|---|---|---|---|
| flash_point_C | 4 | 1 (diesel_bis_is1460 v4pool) | 1 | 2 |
| kin_visc_m20C_cSt (jet) | 4 | 0 | 0 | 0 |
| nhoc_MJ_kg (jet) | 4 | 0 | 0 | 0 |
| density_15C_kg_m3 | 0 | 5 (all gasoline) | 4 | 4 |
| aromatics_wt / vol, ethanol_vol | 6 | 0 | 1 | 1 |
| freezing_point_C | 0 | 0 | 2 | 2 |
| smoke_point_mm | 0 | 0 | 2 | 2 |
| kin_visc_40C_cSt (diesel) | 0 | 0 | 4 | 0 |
| dcn, cloud_point_C, pour_point_C | 0 | 0 | 0 | 8 |
| ron, mon | 0 | 0 | 2 | 8 |
| dvpe_kPa | 0 | 0 | 1 | 4 |
| distillation T50 / E70 | 0 | 0 | 0 | 5 |

Certification scope: every yaml margin was validated on jet fuels (CRC 647 / AJFTD / jet-range lab blends), so a pass with a
margin counts as certified only on a jet spec.  The same pass on a diesel or gasoline spec is reported as
`pass_jet_derived_margin_only` and is excluded from the certified count (6 evaluations: the diesel BIS flash and the 5
gasoline densities).  Only class-fraction and ethanol quantities fixed by the recipe are by_construction for any fuel kind.

## Which properties are certified (jet specs only) and which are not

Individual constraint passes on jet specs (no complete design is certified):

* **Flash point** (Le Chatelier, ideal): margin = max(lab 3.70 C, D93 0.0 C) = 3.70 C.  Matches thermo/chemicals within 1 C
  (verification_v4.csv).  D56 has margin `inf` (n = 1), so a D56-method pass is uncertified.
* **Net heat of combustion**: 0.194 MJ/kg (one-sided, protects the minimum).
* **Kinematic viscosity at -20 C**: recalibrated Refutas, ln-ratio margin 0.3465, valid -40..20 C.
* **Density at 15 C**: margin 22.25 kg/m3 (two-sided) narrows the 775-840 window to about 797-818.  The v4 jet pools reach
  790.9 and 786.7, so jet density fails after the margin.  The gasoline density passes (747.5 against a 720-775 window) use this
  jet-derived margin and are NOT certified.
* **Class mass/volume fractions and ethanol volume**: by construction (no model error).

NOT certified, and why:

* **Diesel anything**: the corpus has 0 diesel fuels, so nothing about diesel is validated.  The single diesel flash "pass"
  (diesel_bis_is1460 v4pool) is `pass_jet_derived_margin_only`.  Diesel 40 C viscosity is also outside the margin's
  validated -40..20 C range.
* **Freezing point** (all jet designs): the margin for an ideal solid-liquid model is `inf` (8 groups).  The Coutinho-library
  alternate has a finite margin (4.05 C) but is validated on CRC real fuels only.
* **Smoke point**: AJFTD bulk-regression margin is `inf` (6 groups).
* **DCN, cloud point, pour point**: excluded from the yaml (`excluded` section); no valid margin or no model.
* **RON / MON**: Task 15 validated octane models are base gasoline plus ONE blendstock, so multi-compound pools have no
  validated octane model.  Several pure RON/MON inputs are unverified (R16).
* **DVPE**: units of the DVPE data are unconfirmed and it is excluded from the yaml, so a DVPE pass is never certified.
  The **minimum is enforced**: a DVPE below the yaml minimum is `fail` (design `gasoline_astm_d4814__v4pool_dvpe`, 42.2 kPa
  against the 45 kPa minimum).
* **Distillation T50 / E70**: no model.
* **Unverified spec limits**: a limit whose spec entry has `verification` other than "verified" (case-insensitive) or
  `unverified: true` can never yield a pass (tested; the check was dead code before this fix).
* **Pools using decalin_* or isopentane**: pure constants are `disputed`/`unverified`, so the model raises `ModelUnavailable`;
  those constraints were excluded from the optimiser, which is why the v3-pool designs collapse to degenerate recipes
  (`n_excluded_from_optimiser` in designs_v4.csv).

### R4 deviation (DVPE minimum in the specs yaml)

Ruling R4 asked for a DVPE min/max per gasoline class from IS 2796 / ASTM D4814 text where confirmable, otherwise an
`unverified: true` value.  **R4 not applied for BIS E10** because I could not cite a source for an IS 2796 DVPE minimum (the
yaml entry has only a 67 kPa maximum), and no number is invented.  The ASTM D4814 entry already has min 45 / max 70 and that
is enforced.  Whoever can cite IS 2796 should add the minimum to `configs/specs_bis.yaml` with `unverified: true` and a source
comment; the engine will then enforce it as a nominal limit and never certify it.

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
without a margin are uncertified; none borrows another form's margin.  For jet freezing the Coutinho library passes on 2
designs while the ideal SvL and ideal-mole forms fail or are uncertified, so the freezing result is not robust to model form.

### Viscosity model choice (disclosure)

The Task 11 winner was the structured ridge (margin 0.162 ln) against recalibrated Refutas (0.3465 ln).  The primary here is
recalibrated Refutas on purpose: it is the conservative choice (wider margin) and, as a 2-parameter physics model, does not
extrapolate.  The ridge is linear in 7 class-fraction features fitted on CRC jet fuels.  It is kept as an alternate behind an
envelope guard that raises `ModelUnavailable` when any feature lies outside the dev-fuel min/max.  That is why the ridge row is
NaN/uncertified for all 8 designs that have a viscosity constraint (jet and diesel; gasoline has none) in robustness_v4.csv:
every pool blend has a feature outside the CRC envelope (for example aromatic fraction 0 against a CRC minimum of 0.060, or
0.44-0.49 against the CRC range).  Reproducible check: `viscosity_envelope_check_v4.csv` and
`test_ridge_envelope_guard_fires_for_pure_decane`.  Pure n-decane at -20 C gives 0.126 cSt from the unguarded ridge (the guard
refuses, NaN), against 3.26 cSt recalibrated Refutas and 2.30 cSt raw Refutas.

## Independent verification (verification_v4.csv)

thermo/chemicals recompute Psat, Le Chatelier flash, ideal SvL liquidus and ideal-volume density (no fuelblend import,
enforced by `tests/test_v4_verify_imports.py`).  17 comparisons agree; 8 designs/properties were skipped (engine value not
available, so nothing to compare).  0 disagreements exceed the margin.  Where the engine margin is `inf`/absent (freezing) the
tolerance is a fixed 5 C stated in the script (`tolerance_is_engine_margin = False`); the freezing differences are 3.5-3.7 C.
This verifies the arithmetic only.  It says nothing about the validity of the ideal-mixture assumption.

## v3 vs v4, side by side (comparison_v3_v4.csv)

v3-pool rows: constraints involving disputed or unverified constants (decalin_*, isopentane) were EXCLUDED from the optimiser,
so those redesigns are degenerate and say nothing about whether a good blend from the v3 pool exists.

| design | v3 | v4 re-score of the same recipe | v4 redesign (same pool) | why it changed |
|---|---|---|---|---|
| jet ASTM / BIS | pass | not certified: 3-4 constraints uncertified (decalin density/dHfus disputed); no violation among computed properties | degenerate (100% n-decane), 3-4 constraints uncertified | v3 used formula-derived properties; v4 refuses unverified constants |
| diesel ASTM / BIS | pass | fail: nominal 1.81 cSt at 40 C is below the limit; 4 uncertified | degenerate (100% n-hexadecane) | v3 viscosity at 40 C was not validated |
| gasoline ASTM | pass | not certified: RON/MON/DVPE/T50 uncertified; density pass is jet-derived only | same | no validated octane/DVPE for this pool |
| gasoline BIS E10 | pass | not certified: same | same | same |

With a pool of verified compounds ("v4pool"): jet fails at the point estimate on freezing (-30.5 and -39.3 C against the
limits) and smoke point, and fails density only after the margin.  Diesel fails at the point estimate on 40 C viscosity (and on
aromatics for ASTM); its density and flash fail only after the margin.  Gasoline has no limit that can be tested
(octane and DVPE have no certifiable model).
