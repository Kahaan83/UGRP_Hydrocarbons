# v5 BASE: full-baseline scoreboard

Status: DECLARED 2026-11-12, not scored.
Start commit: see `START_COMMIT.txt`. Plan: `docs/superpowers/plans/2026-11-12-v5-full-baseline.md`.

Every T2 fuel was read before this declaration; the baseline is descriptive.
The CRC 647 locked 12 fuels have been read five times; they are reported as `fold_b`, not as a locked test.

## Scoreboard

Rows = property x tier x set x model x variant. Columns:
`property, tier, set, model, variant, status, n_groups, n_rows, bias, mae, rmse, unsafe_rate, trivial_model, mae_trivial, diff_trivial_minus_model, ci_lo, ci_hi, verdict, margin90, coverage90, n_excluded_censored, n_excluded_nonfinite, n_low_coverage, evidence_path, label`

- `tier`: `T1_lab | T2_measured | CRC647_r3 | T3_ajftd | diesel | designs`
- `status`: `recomputed_v5 | aggregated_v5 | carried_v4 | no_evaluation`
- `label`: `descriptive | dev | fold_b | already_read | exploratory`
- `verdict`: `beats_trivial | no_skill_shown | worse_than_trivial | n_too_small | not_applicable`.
  Rule: paired group bootstrap (5000 resamples, seed 0), 95 % interval of `MAE(trivial) - MAE(model)`;
  above 0 = beats, below 0 = worse, else no skill shown; fewer than 8 groups = `n_too_small`.
- `trivial_model`: leave-one-group-out constant. For flash and heat of combustion the constant is fitted within method family.

## Definitions

- Error = predicted - measured. Viscosity error = ln(pred) - ln(meas).
- Unsafe side: flash = fraction over-predicted by more than 5 C. Viscosity at T < 283 K = fraction under-predicted
  by more than the S1 margin (0.4154 below 263 K, 0.6084 at 263-283 K).
- Group rule: one fuel per double-count group, groups from `results/v5/S3/corrected_2026-11-12_groups.csv`
  (source order AV-33-22, FuelLib, AJFTD).
- Censoring rule: values with a non-null `qualifier` (e.g. "<") are bounds, excluded and counted.
- Finite rule: a prediction with status `ok` must be finite; others are excluded and counted in `n_excluded_nonfinite`.
- `low_coverage`: a fuel whose composition maps under 90 % of its mass; scored apart, never in the headline, never dropped silently.
- Several measured values for one (fuel, property, T): one median target; method family recorded.
- Designated model per tier: T1/T2 = API `composition_kind="known"`; CRC 647 r3 = `composition_kind="reconstructed"`;
  T3 = no shipped model (hgb, ridge vs constant); freezing, smoke, pour, cetane = named baseline in the `model` column.
- Locked stages are never re-invoked; locked cells are aggregated from committed CSVs and labelled `already_read`.

## Cells that must exist

| Property | T1 lab | T2 measured | CRC 647 r3 | T3 AJFTD | Diesel |
|---|---|---|---|---|---|
| Viscosity by T bin | aggregate S1 | recompute, 9 T | carry v4 R3 | aggregate S4 | aggregate S6 |
| Density by T | v4 L + ThermoML 2020 new systems | recompute, 9 T | carry | no_evaluation | carry |
| Flash point | v4 F, 126 rows | recompute | carry | aggregate | carry |
| Freezing point | SLE2 tiers | aggregate S5 | carry | aggregate | no_evaluation |
| Net heat of combustion | no_evaluation | new | carry | S4 + S4b | no_evaluation |
| Smoke point | no_evaluation | new (YSI) | no_evaluation | carry v4 SP | no_evaluation |
| Distillation T10/T50/T90 | no_evaluation | new | no_evaluation (input) | no_evaluation (input) | no_evaluation |
| Pour point | no_evaluation | no_evaluation | bulk baseline n=57 | no_evaluation | no_evaluation |
| Cetane / DCN | no_evaluation | no_evaluation | no_evaluation | no_evaluation | S6 not evaluable |
| RON / MON, DVPE | carry v4 O | - | - | - | - |
| Designs certified | aggregate D4 | | | | |

## Variants

`designated`, `lfl_overlay`, `mesitylene_fixed`, `lfl_overlay+mesitylene_fixed`, `k1_cold_offset_loso`,
`sle2_two_compilations`, `sle2_single_compilation`. The `designated` variant stays the v4 routing until the user decides.

## Readings (declared before any BASE result exists)

- R-BASE-1: every cell in the table exists with a status.
- R-BASE-2: each `recomputed_v5` cell that has a committed counterpart reproduces it within tolerance
  (flash 0.05 C, ln viscosity 0.005, density 0.5 kg/m3, freezing 0.05 C) or the difference is listed with its cause.
- R-BASE-3: per property, name the tier where the designated model beats the trivial model and the tier where it does not.
- R-BASE-4: per property, name the largest unsafe-side rate and the set it occurs on.
- R-BASE-5: list every cell whose verdict is `n_too_small` or `no_evaluation`; this is the data-gap list.

## Note 2026-11-12 (Task 4, committed before the density cell is scored)

Reading for the one new T1 cell (ideal-volume density on ThermoML 2020 systems absent from the unified table and from the
v2 ThermoML extract, selected by system_id and by DOI): **absolute bias <= 2 kg/m3 and 90 % of systems inside the v4 lab density
margin** (4.5866 kg/m3 over, 3.2499 kg/m3 under). These systems have never been scored for density, so the cell is labelled `dev`
("first read").

## Note 2026-11-12 (Task 6 readings, heat / smoke / distillation)

Readings (text fixed by the plan, Task 6):
* heat: group additivity beats a within-method constant on measured rows (n >= 8 groups).
* smoke and distillation: reported, no skill claim; the verdict comes from the rule.

Provenance, stated plainly: the plan required this note to be committed before scoring. It was not. `rows_t2_bulk.csv`
was generated earlier on 2026-11-12, and the heat cells were printed in the session before this commit. This reading is
therefore NOT a pre-registration. B6 is not fixed for these three cells. The reading text is the plan's, unchanged.

## Note 2026-11-12 (SCORED)

All Part B cells exist (`cells_t1_lab`, `cells_t2_flash`, `cells_t2_bulk`, `cells_t3_bulk`, `cells_carried`) and are assembled by
`scripts/v5/base/build_scoreboard.py` into `scoreboard.csv`, `data_gaps.csv` and `SCOREBOARD.md`. The declaration text above is unchanged.
R-BASE-1 to R-BASE-5 are answered in `SCOREBOARD.md` from the data.

What was and was not run on this machine:
* `reproduce_base.py` (default) rebuilds from the committed `rows_*.csv`: `cells_t2_bulk.py --from-rows` refits only the K1 offsets, so the
  T2 viscosity, density, heat, smoke and distillation predictions in `rows_t2_bulk.csv` were generated once and not regenerated in the default run.
  `--full` re-runs every prediction stage; the T2 prediction pass alone took over 80 minutes, so `--full` has NOT been run end to end.
* R-BASE-2: the AJFTD GCxGC flash headline (31 fuels, one per double-count group) differs from G1 (37 rows, bias -0.426, MAE 3.536) by 0.21 C in bias.
  The cause is the group rule, shown by re-scoring all 37 rows from `rows_t2_flash.csv`, which reproduces G1 exactly.
* Every T2 cell is `descriptive` (all T2 fuels have been read). CRC 647 `fold_b` was read five times before. Neither is a held-out number.
* The readings for Task 6 cells 4-6 (heat, smoke, distillation) were committed after those results were seen (see the Task 6 note above).
* Smoke point (0 of 42 fuels) and distillation (1 of 78 fuels) on T2 are `no_evaluation`: mass coverage of the YSI table and of the Antoine boiling points is under 0.90.

## Note 2026-11-12 (frozen reference)

The frozen scoreboard is the one committed in `5a62e20`. Part C deltas are measured against it; `scoreboard.csv` is never edited by a delta.
