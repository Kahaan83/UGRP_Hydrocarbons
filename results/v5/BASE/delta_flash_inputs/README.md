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
