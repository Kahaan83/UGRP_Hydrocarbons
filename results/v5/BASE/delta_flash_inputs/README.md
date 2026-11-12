# delta_flash_inputs (Task 10)

Declared 2026-11-12 before `delta_flash_inputs.py` exists; this commit contains no result file.
Frozen reference: scoreboard commit `5a62e20`. `scoreboard.csv` is never edited.

User decision U1 (2026-11-12): adopt the LFL overlay; skip the corrected mesitylene Antoine set. Mesitylene variants are still
scored, as rejected side-by-side rows, so the decision stays checkable.

Method: registry overrides through `RouteCtx` (in memory only). Lab mixtures (126 rows, 31 systems) are re-predicted per row for four
variants (designated, lfl_overlay, mesitylene_fixed, lfl_overlay+mesitylene_fixed). The lab flash margin is recomputed with
`one_sided_margin(direction="over", alpha=0.10)` for each variant, over all 31 systems and with the mesitylene system left out; the v4 value is 3.70 C.
CRC 647 r3 and T2 rows come from the committed `L1/F_r3/flash_real_fuels.csv` and `rows_t2_flash.csv`, compared with the frozen designated rows.

Readings (fixed now):
1. With the LFL overlay, no set's unsafe rate (over-prediction > 5 C) rises against designated: lab, CRC 647 dev, CRC 647 fold_b, T2 pooled and each T2 source.
2. The leave-mesitylene-system-out lab margin is within 0.5 C of the full margin (checked for designated and for lfl_overlay).
3. The overlay changes the lab margin by less than 0.5 C versus 3.70 C; otherwise the new value is reported and adoption is the user's call again.

If readings 1 and 2 hold, a follow-up commit writes the overlay flash margin once into a new `configs/uncertainty_v5.yaml` (v4 yaml untouched).
If either fails, nothing is written to `configs/` and the failure is reported.

## Result 2026-11-12

Run with `scripts/v5/base/delta_flash_inputs.py`; numbers in `delta.csv`, `lab_margin.csv`, `readings.json`.
* Reading 1 PASSES: with the overlay no unsafe rate rises. Lab 0.8 % -> 0.8 %; CRC 647 dev 13.3 % -> 4.4 %; fold_b 25 % -> 8.3 %; T2 pooled 5.7 % -> 5.7 % and unchanged in every source and family.
* Reading 2 FAILS. The lab margin 3.70 C is the group score of the mesitylene system itself (it is the 29th of 31 ranked group scores). Leaving that
  system out gives 0.31 C (30 systems), a 3.39 C difference, for designated and overlay alike. With the mesitylene Antoine set the all-31 margin is also 0.31 C.
  So the v4 lab flash margin is set by one system with a suspect Antoine input, not by the overlay.
* Reading 3 PASSES trivially: the overlay leaves the lab margin at 3.70 C (the overlay does not change the tail group scores).
* Overlay effect on bias/MAE: lab bias -3.30 -> -3.25 C, MAE 3.74 -> 3.69; CRC dev bias +0.44 -> -0.73 C; T2 pooled bias -0.73 -> -0.69 C, MAE unchanged (2.58). The mesitylene set worsens
  T2 pooled MAE to 2.81 C and lab MAE to 3.89 C, consistent with U1 (not adopted).

Consequence: `configs/uncertainty_v5.yaml` is NOT written. The overlay does not move the margin, so writing a v5 file would only restate 3.70 C, and the margin's
dependence on one system is unresolved. If you want the margin settled, decide between keeping 3.70 C (conservative, rests on the suspect mesitylene system)
or the 0.31 C leave-one-system-out value (not safe to adopt: it rests on 30 mostly binary alkane systems and shows how fragile the 29th-of-31 rank is).
