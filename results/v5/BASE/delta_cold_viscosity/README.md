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

## Result 2026-11-12

Run with `scripts/v5/base/delta_cold_viscosity.py`; numbers in `delta.csv`, `proposed_margins.yaml`, `pure_cold_extrapolation.csv`, `readings.json`.
Group counts are those of the frozen scoreboard (AV-33-22 35 at both temperatures; AJFTD 28 at -20 C and 10 at -40 C; the plan's 19 and 8 were stale).
* Reading, first half PASSES: the K1-corrected two-sided margin at -20 C is 0.172 (pooled, 65 groups), below the raw S1 margin 0.4154. At -40 C it is 0.293 (48 groups).
* Reading, second half FAILS: the margin set on AV-33-22 residuals (0.152 at -20 C) covers only 21 of 28 AJFTD groups (75 %), under the 90 % target.
  The margin set on AJFTD residuals (0.286) covers all 35 AV-33-22 groups. At -40 C the AV-33-22 margin (0.293) covers 9 of 10 AJFTD groups (90 %) and the AJFTD margin (0.434) covers all 35.
  So the pooled 0.172 is NOT shown to transfer to AJFTD; the AJFTD-set 0.286 is the transferable value at -20 C. FuelLib (4 groups) has margin `inf`.
  `proposed_margins.yaml` holds all of these; nothing is in `configs/`, and the margin is not adopted.
* Cause check: the pure provider is on `ml_extrapolation` (outside its own measured range by more than 25 K, or no measured fit) for 50.7 % of mapped mass at 253 K and 59.9 % at 233 K,
  so by the declared rule the cold-viscosity fix belongs in the pure provider. Next delta: Walther extrapolation from each compound's own measured range.
  Only 2 of 17 representatives have a Boehm pure value at 253 K (ln errors -0.075, -0.177) and 1 at 233 K (-0.460), so the cause is not confirmed by pure data for most of the mass; the 50.7 % at 253 K is only just over half.
