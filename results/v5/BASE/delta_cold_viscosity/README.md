# delta_cold_viscosity (Task 11)

Declared 2026-11-12 before `delta_cold_viscosity.py` exists; this commit contains no result file.
Frozen reference: scoreboard commit `5a62e20`. `scoreboard.csv` is never edited. User decision U1: adopt the K1 cold offset below 263 K
as a variant; this delta only measures whether it can carry a margin.

Method. The K1 offset is the per-temperature mean ln(meas/pred) of the other source, refit with the scored double-count group removed
(AV-33-22 scored with an AJFTD-trained offset, AJFTD with an AV-33-22-trained one; FuelLib with AV-33-22). Residual = ln(pred_K1) - ln(meas), one fuel
per group at each temperature. The two-sided margin is `two_sided_margin(|residual|, groups, alpha=0.10)` at 233.15 K (-40 C) and 253.15 K (-20 C),
pooled over the scored sources and per scored source. `inf` when fewer than 9 groups. Written to `proposed_margins.yaml` in this folder, not to `configs/`.

Reading (fixed now): "the corrected-model margin at -20 C is below the raw S1 margin 0.4154 and 90 % coverage holds on the other source",
i.e. the margin set on AV-33-22 residuals covers at least 90 % of AJFTD groups at -20 C and the margin set on AJFTD residuals covers at least 90 % of AV-33-22 groups.
The -40 C margin is reported with the same two checks; where AJFTD has fewer than 9 groups its margin is `inf` and that check is `not evaluable`.

Cause check (audit proposal 2): for every AV-33-22 representative with more than 2 % mean mapped mass, record whether the pure-viscosity provider is inside its
measured fit range at 253.15 K and 233.15 K (`fit` within the 25 K tolerance of its range, otherwise `ml_extrapolation`) and its ln error against Boehm 2022
pure data where present. `pure_cold_extrapolation.csv` holds this. If more than half of the mapped mass is on `ml_extrapolation` at a temperature, the README result note
says the fix belongs in the pure provider and names Walther extrapolation from each compound's own measured range as the next delta.
