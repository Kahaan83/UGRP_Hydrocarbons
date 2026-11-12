# delta_t2_residual (Task 13)

Declared 2026-11-12 before `delta_t2_residual.py` exists; this commit contains no result file. Frozen reference `5a62e20`. `label = exploratory`.

Model: physics baseline + ridge residual on 22 composition features (`fuelblend.mixture.t2_residual.FEATURE_NAMES`; the plan said 20 names but listed 22,
all 22 are used). Features standardised on training rows only; alpha in (0.1, 1, 10, 100) by inner leave-one-group-out on training rows only.
Targets: residual (measured - designated physics; ln for viscosity) of flash, density at 15 C, ln viscosity at 40, -20 and -40 C, heat of combustion (ga_comp_H),
freezing (S5 F1). Validation: train AV-33-22 -> test AJFTD GCxGC, and the reverse; FuelLib test-only (trained on both). Training rows whose double-count
group appears in the test set are removed. One fuel per group.

Reading per property (fixed now): "the residual model beats raw physics in both directions with the bootstrap interval above 0 and does not raise the unsafe rate".
Scored with the paired group bootstrap (raw physics as the comparator); fewer than 8 test groups -> n_too_small, which fails the reading for that direction.
