# delta_t2_residual (Task 13)

Declared 2026-11-12 before `delta_t2_residual.py` exists; this commit contains no result file. Frozen reference `5a62e20`. `label = exploratory`.

Model: physics baseline + ridge residual on 22 composition features (`fuelblend.mixture.t2_residual.FEATURE_NAMES`; the plan said 20 names but listed 22,
all 22 are used). Features standardised on training rows only; alpha in (0.1, 1, 10, 100) by inner leave-one-group-out on training rows only.
Targets: residual (measured - designated physics; ln for viscosity) of flash, density at 15 C, ln viscosity at 40, -20 and -40 C, heat of combustion (ga_comp_H),
freezing (S5 F1). Validation: train AV-33-22 -> test AJFTD GCxGC, and the reverse; FuelLib test-only (trained on both). Training rows whose double-count
group appears in the test set are removed. One fuel per group.

Reading per property (fixed now): "the residual model beats raw physics in both directions with the bootstrap interval above 0 and does not raise the unsafe rate".
Scored with the paired group bootstrap (raw physics as the comparator); fewer than 8 test groups -> n_too_small, which fails the reading for that direction.

## Result 2026-11-12

Run with `scripts/v5/base/delta_t2_residual.py`; per-direction rows in `delta.csv`, verdicts in `readings.json`. Exploratory; nothing adopted.
* PASS (both directions beat raw physics with the bootstrap interval above 0, unsafe rate not raised): ln viscosity at -20 C and at -40 C.
  At -40 C AV-33-22 MAE 0.637 -> 0.112 and AJFTD 0.698 -> 0.425 (10 groups); unsafe rate 91 % -> 0 % and 100 % -> 30 %. This is mostly the intercept, i.e. the same
  cold offset K1 already provides; the composition terms add little (compare `delta_cold_viscosity`).
* FAIL: flash (no skill either way; MAE worse, 1.50 -> 1.91 and 3.30 -> 4.46 C), density at 15 C (worse on AV-33-22, 2.47 -> 7.86 kg/m3), viscosity at 40 C
  (better on AV-33-22, worse on AJFTD), heat (better on AJFTD only), freezing (worse on AV-33-22).
* FuelLib (3-4 groups) is always n_too_small.
Reading: the learned residual transfers between sources only where the physics error is a near-constant offset (cold viscosity). With 35 and 16-35 fuels per source it does not
learn composition-dependent corrections that survive the source change.
