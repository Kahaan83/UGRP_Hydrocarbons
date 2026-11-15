# delta_adopt_rescore: BASE re-score under the adopted defaults (decisions 2026-11-15)

Declared 2026-11-15 BEFORE the long run (`reproduce_base.py --full` has not been started for this folder). Base: v5 at commit `fc618b5`.
Plan: `docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md`, Task 10. Spec of the decisions: `docs/decisions_2026-11-15.md`.
Frozen reference: `results/v5/BASE/**` (excluding every `delta_*` folder) and `SHA256SUMS_v5_base.txt`, both from commit `5a62e20`.
Gates, margins, spec limits, the scoring rules of `README.md` (parent) and `configs/uncertainty_v5.yaml` are not changed by this delta.

## What this delta measures

Scope = every file under `results/v5/BASE/**` (non-`delta_*`) that `python scripts/v5/base/reproduce_base.py --full` regenerates with a
SHA-256 different from `SHA256SUMS_v5_base.txt`, and, inside `scoreboard.csv` / `cells_*.csv` / `SCOREBOARD.md`, every cell whose value differs.
For each changed cell the delta lists the old value, the new value and the decision that caused it. The frozen BASE is the **old** side;
the full re-run of the current code is the **new** side. Readings are reported, not judged against gates (no gate is moved).

The parent `README.md` R-BASE-2 tolerances (flash 0.05 C, ln viscosity 0.005, density 0.5 kg/m3, freezing 0.05 C) are used only to label a
difference "inside / outside the frozen reproduction tolerance"; a cell is "changed" at any difference > 1e-6 (absolute, or relative for
values above 1), because the CSVs are written with 8 significant digits. A hash difference with no cell difference is listed as `format_only`.

## Decisions that can move numbers

| # | decision | commit | where it acts |
|---|---|---|---|
| 1 | designated freezing route = UNIFAC Schroeder-van Laar (`DESIGNATED_FREEZING_MODEL = "svl_unifac"`) | `bcab8b4`, `9fbff24` | `composition_route.freezing_C`, `predict_fuel["freezing_point_C"]`, design judge |
| 3 | flash-registry additions overlay (bicyclohexyl, cis-decalin Antoine + LFL) ON by default | `5986a6b` | `flash_unifac._default_registry()` -> `_lookup`, `RouteCtx.flash_reg`, `SubReps(vp_lookup=_lookup)` |
| 4 | bicyclohexyl SvL dHfus = 13.94 kJ/mol (fusion + II->I) | `3ec5aa1` | `tm_dhfus_verified.csv` -> `build_tmap_verified` -> freezing only |
| 5 | exact aromatic matching in the Coutinho table (`build_tmap`) | `0b248e6` | `build_tmap` -> freezing only |
| 7 | cis-decalin MW verified (`pure_props_v5.csv`) | `3ec5aa1` | `design/library_route.py` only (design library) |

Decisions 2 (D86 stays Engler), 6 (`strict_constraints=True`, design only) and the margin decisions do not act on any BASE stage.

A sixth class exists that is not a 2026-11-15 decision: **PRE**, the code changed between the frozen commit `5a62e20` and the decisions
(for example `delta_api_route` moved `predict()` onto the BASE route; `delta_ml_viscosity` made `ml_pair_residual` the designated in-domain viscosity;
`delta_flash_margin_q2`, `delta_smoke_margin`). The T2 `api_as_shipped` rows call `predict()`, so PRE can move them. A change that survives with
decisions 1 and 3 reverted is attributed to PRE; this is stated in `cause`, not hidden.

## Attribution rule (declared now; a full re-run per switch takes hours and is not done)

Cause is assigned in this order:

1. **By code path (`by_code_path`)**: which switch can reach the cell, read from the code. Declared reading of the code before the run:
   * Stage `build_t2_targets`, `cells_t2_flash`, `cells_t2_bulk` build `SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))`; `_lookup` reads
     `_default_registry()`, i.e. the overlay (#3). Fuels whose GCxGC bins have bicyclohexyl or cis-decalin as the representative, or contain them, can change
     in **every** property (representative choice, mapped coverage, `low_coverage`, flash, Refutas, ideal volume, heat, smoke, D86). `RouteCtx.flash_reg`
     is a copy of the same registry (flash predictions).
   * `cells_t1_lab`, `cells_t3_bulk`, `cells_carried`, `build_scoreboard` read committed CSVs/YAML (`results/v4/**`, `results/v5/S1,S2,S4,S4b,S5,S6,SLE2,L1,D4,G1`);
     their only recomputed part is `cells_t1_lab` density on ThermoML 2020 systems (`PurePropertyProvider` + `CompoundRegistry`, no overlay, no freezing).
     **No decision in the table above reaches them.**
   * **Every BASE freezing cell is an aggregation of a committed CSV** (`results/v5/S5/s5_predictions.csv`, `results/v5/SLE2/summary_by_tier.csv`,
     `results/v4/R3/B/freezing_v4.csv`, S4, S6). No BASE script calls `freezing_C`, `ideal_sle_C` with a tmap, or reads `predict_fuel()["freezing_point_C"]`
     (`RouteCtx(reg, prov)` is built with an empty `tmap`). So decisions 1, 4 and 5 **cannot** change a freezing cell; decision 7 is design-only.
   * The one indirect path for 1, 4, 5: `predict()` calls `predict_fuel`, which now also computes freezing. `cells_t2_flash` / `cells_t2_bulk` wrap `predict()`
     in `try/except` and store NaN with `api_method = "error:..."` on any exception. If a freezing/UNIFAC/tmap error turned a viscosity or density API row into
     an `error:` row, that row would be a #1/#4/#5 effect. The delta checks `api_method` / `dmethod` / `api_error` for `error` in the new `rows_t2_bulk.csv` and
     `rows_t2_flash.csv` and reports the count (expected 0).
2. **Targeted recompute for ambiguous cells ONLY.** A cell is ambiguous if, by code path, more than one of {#3, PRE, (#1 via the error path)} can reach it.
   Its fuels/rows are recomputed in a scratch process (outputs in the scratchpad, never in `results/`) with exactly one switch reverted:
   * `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` (plus `_default_registry.cache_clear()` / a fresh process) reverts #3;
   * `DESIGNATED_FREEZING_MODEL = "ideal_svl"` reverts #1 (and `predict_route._tmap_freezing.cache_clear()`);
   * #4, #5, #7 have no revert switch and no declared path into a BASE cell: attributed `by_code_path` only, labelled so.
   Cause labels in `drift.csv`: `D3_overlay`, `D1_freezing_route`, `D4_bicyclohexyl_dHfus`, `D5_coutinho_exact`, `D7_decalin_MW`, `PRE_<commit-area>`,
   `format_only`, `unattributed`. Suffix `:by_code_path` or `:by_recompute`. A cell attributed to D1/D4/D5/D7 is **unexpected** under the reading above and is
   flagged as such. `unattributed` is allowed and is reported as unattributed.
3. A cell that is unchanged is not listed (it appears as a count per `cells_*.csv`).

## Declared expectations (reported, not gates)

* E1: every changed cell traces to D3 (overlay) or PRE; none to D1, D4, D5, D7.
* E2: all freezing cells (T1, T2, T3, CRC 647, diesel) are bit-identical to the frozen scoreboard.
* E3: cells that move are those of fuels whose mapped composition touches bicyclohexyl or cis-decalin.

Honest disclosure, as this is not a blind pre-registration: before this declaration, `tests/test_predict_route.py::test_heat_matches_bulk_rows` and
`tests/test_base_composition_route.py::test_reproduces_s3_corrected_M_predictions` were run on the current code and **fail**, both on AJFTD:2
(heat_net_MJ_kg 43.5287022 live vs 43.508441 in frozen `rows_t2_bulk.csv`; density 793.08801 vs 791.99300 kg/m3 and viscosity 4.764561 vs 4.60734 cSt at 233.2 K,
3.374664 vs 3.27281 at 253.2 K in the S3 corrected predictions). 110 other tests across 13 BASE-related test files pass. The delta therefore starts with the knowledge
that AJFTD:2 moves; E1-E3 were set after that.

## Rules

* **The frozen BASE files and `SHA256SUMS_v5_base.txt` are NOT overwritten in this task.** `--full` writes in place into `results/v5/BASE/`, so a copy of the
  frozen files is taken first (`BASE_frozen`, 23 files, 22 of them listed in the manifest, SHA-256 verified against the manifest) and is used to restore them.
* Before restoring, the regenerated `SCOREBOARD.md` and all regenerated source CSVs/JSON are copied into `delta_adopt_rescore/rerun/` and committed, so the
  new complete numbers are in git. `rerun/` is **not** a frozen BASE and is not hashed into the manifest.
* `--update-manifest` is never passed. Adopting the re-run as the new frozen BASE is gated on the user (question U4), not decided here.
* After restoring, `reproduce_base.py` (default, fast) must show 0 changed and 0 missing. One `NEW` line is **pre-existing and not drift from this delta**:
  `results/v5/BASE/data_thermoml_sle/README.md` (added in `23e4aa7`, never in the manifest; checked read-only with the manifest functions before this declaration).
* Outputs of this delta: `README.md` (this declaration, then results appended below the line), `drift.csv` (`path, cell, old, new, cause`), `rerun/`.

## Pinned tests and the legacy switches

Tests that pin frozen BASE/S3 values against the live route and fail (or can fail) under the adopted defaults will be pointed at the **legacy switch with a
fixture**; no expected value is edited. Known now: `test_predict_route.py::test_heat_matches_bulk_rows` and
`test_base_composition_route.py::test_reproduces_s3_corrected_M_predictions` (fixture: `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` via `monkeypatch.setenv` with
`flash_unifac._default_registry.cache_clear()` and `predict_route.get_ctx.cache_clear()` before and after, as `tests/test_delta_pool_bicyclohexyl.py` does).
The final list comes from the re-run (tests that read the regenerated files in `rerun/` are not touched; tests that read the frozen files keep reading them).
A new test will assert the live-vs-frozen difference is confined to the cells listed in `drift.csv`. Test edits are made in a later commit, after `drift.csv` exists.

## Long run (for the controller)

`.venv\Scripts\python.exe scripts\v5\base\reproduce_base.py --full` in the background, `FUELBLEND_FLASH_REGISTRY_ADDITIONS` unset, no `--update-manifest`.
It runs, in order: `build_t2_targets.py`, `cells_t1_lab.py`, `cells_t2_flash.py`, `cells_t2_bulk.py` (full prediction pass, > 80 min on the development machine
per the parent README), `cells_t3_bulk.py`, `cells_carried.py`, `build_scoreboard.py`, then `pytest tests/test_base_*.py`, then the manifest check.
A pytest failure ends the script with exit 1 BEFORE the manifest check; the files are already written, so drift is computed against `BASE_frozen` directly.
`predict_route` reads `results/v5/BASE/rows_t2_bulk.csv` and `rows_t2_flash.csv` (K1 offsets, indicative intervals), so the run must not overlap any process that
loads them (the `delta_adopt_designs` runs).

## Results

(to be appended after the run; the declaration above is not edited)

Appended 2026-11-15 after the run. The declaration above is unchanged.

### What was run

* `reproduce_base.py --full` ran to the end, with the adopted defaults and `FUELBLEND_FLASH_REGISTRY_ADDITIONS` unset. It exited 1 at the pytest stage, with the 4 pinned-value failures listed below, so the manifest check did not run. Log: `.superpowers/sdd/2026-11-15-apply-decisions-and-audit-fixes/reproduce_base_full.log`.
* The regenerated files were copied into `rerun/` **before** anything was restored:
  * `SCOREBOARD.md`, `scoreboard.csv`, `data_gaps.csv`, `cells_*.csv`, `rows_*.csv`, `t2_*.csv/json`, `*_error_vs_T.csv` (21 files);
  * the 11 tracked `delta_*` files that the pytest stage rewrote with live values, under `rerun/delta_files/<folder>/`.

  `rerun/` is not a frozen BASE and is not in the manifest.
* **Targeted recompute (step 2 of the attribution rule).**
  * What was run: `build_t2_targets`, `cells_t2_flash` and `cells_t2_bulk` (full prediction pass, 7.4 min), in a scratch process with exactly one switch reverted: `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0`, set before any import.
  * Where it read and wrote: the module `OUT` was redirected to the scratchpad. `predict_route` read the restored frozen `rows_t2_*.csv`.
  * Scoreboard: the overlay-off scoreboard (`X`) was assembled the same way, with `build_scoreboard.BASE` pointed at a scratch copy of the frozen files plus the recomputed T2 files.
  * Scope: the recompute covered every T2 row, not only the ambiguous cells. It costs minutes, not the hours the declaration assumed.
  * Result: it reproduces the frozen designated rows exactly. Over 1,736 `rows_t2_bulk.csv` designated, K1, heat, smoke and D86 rows, the maximum relative difference is 9.4e-8. All 79 `rows_t2_flash.csv` `pred_*` values match, and so does `t2_fuels.n_components`. So decision 3 is the only cause on the designated route.
* The pure-property grid matched its pinned SHA-256 in the recompute (`grid_sha_status` True).

### Files and cells that changed

* **Files.** 11 of the 21 regenerated files have a SHA-256 different from `SHA256SUMS_v5_base.txt`: `SCOREBOARD.md`, `scoreboard.csv`, `cells_t2_flash.csv`, `cells_t2_bulk.csv`, `rows_t2_flash.csv`, `rows_t2_flash_api_differences.csv`, `rows_t2_flash_outliers.csv`, `rows_t2_bulk.csv`, `t2_fuels.csv`, `viscosity_error_vs_T.csv` and `density_error_vs_T.csv`.
  * Bit-identical: `cells_t1_lab`, `cells_t3_bulk`, `cells_carried`, `data_gaps`, `rows_pour_point`, `rows_t1_density_thermoml2020`, `t2_targets`, `t2_censored_excluded` and `t2_build_log`.
  * `format_only`: none.
* **scoreboard.csv.** 2,021 cells changed in 294 of 781 rows. The diff key is property/tier/set/model/variant/label/T_K **plus** cells_file/status/evidence_path. Without the last three, four duplicate-key pairs (T3 density, T3 flash hgb/ridge) get paired by row order and show as false drift.
* **drift.csv.** One row per changed cell in `scoreboard.csv`, `cells_t2_flash.csv`, `cells_t2_bulk.csv`, `t2_fuels.csv` and both `*_error_vs_T.csv` files. The `rows_*` files are summarised as one row per (file, property/model/variant, column, cause), with `n_changed` and `max_abs_diff`. For viscosity `pred`/`lo`/`hi` that difference is ln(new/old). `SCOREBOARD.md` gets one `derived` row (348 changed lines).
  * Extra columns: `kind`, `designated`, `n_changed`, `max_abs_diff`, `tol_rbase2`, `outside_frozen_tol`, `overlay_off` (the X value) and `note`.

Cause of the 2,021 changed scoreboard cells:

| cause | cells | meaning |
|---|---|---|
| `D3_overlay:by_code_path` | 949 | designated-route cells (Refutas, ideal volume, Le Chatelier and its input variants, K1, heat); overlay-off recompute reproduces the frozen value |
| `PRE_api_route+D3_overlay:by_recompute` | 818 | `api_as_shipped` cells: old -> X is the API change (delta_api_route), X -> new is the overlay |
| `PRE_api_route:by_recompute` | 246 | `api_as_shipped` cells where the overlay-off value already equals the new one (verdicts, notes, method counts, unsafe rates) |
| `PRE_flash_margin_q2:by_recompute` | 4 | `api_as_shipped` flash `margin90` NaN -> 2.2069 (delta_flash_margin_q2) |
| `D3_overlay:by_recompute` | 4 | `api_as_shipped` reconstructed-route viscosity verdicts (av3322 and pooled at 273.15 K and in the 263-283 K bin, `beats_trivial` -> `no_skill_shown`) that the overlay alone flips |
| D1 / D4 / D5 / D7 | 0 | |
| `unattributed` | 0 | |

How combined labels are read: `A+B:by_recompute` means the overlay-off value X lies strictly between old and new (to 1e-5 relative). The step old -> X is A, and X -> new is B. The X value is in the `overlay_off` column.

Outside the scoreboard, 6 cells of `density_error_vs_T.csv` carry `PRE_from_rows_precision:by_recompute`:
* The difference is at most 1.3e-5 relative, for example `-0.15446` -> `-0.154462`.
* X equals the new value.
* The frozen file was scored `--from-rows`, from the 8-significant-digit CSV and written at 6 digits. `--full` scores the in-memory predictions.

This label is not in the declared vocabulary. It is named here rather than hidden under `unattributed`.

### Declared expectations

* **E1 holds.** Every changed cell traces to D3 or PRE. None traces to D1, D4, D5 or D7.
* **E2 holds.** Every freezing cell, and every T1, T3, CRC 647 and carried cell, is bit-identical to the frozen scoreboard.
* **E3 holds, with one qualification.**
  * Using the mapped composition under the overlay, 144 of 155 T2 fuels contain bicyclohexyl or cis-decalin. The median mass fraction is 0.088 and the maximum is 0.31. `n_components` rises by 1 in 133 fuels.
  * Every scored fuel that contains them moved (74 of 74). None of the 70 unscored ones did.
  * 6 fuels without them moved **only** in the `refutas_raw+k1` variant: AJFTD:161, 174, 189, 93, 97 and FuelLib posf11498. The reason is that the K1 offset is fitted on the other source's designated rows, and those rows moved.
* **Declared `error:` check: 0.** No `error` appears in `api_method` in either `rows_t2_bulk.csv` or `rows_t2_flash.csv`, and neither file has an `api_error` column. So decisions 1, 4 and 5 have no path into a BASE cell.
  * Separately, two API rows became `non_finite` with `api_method = "none"`: AV-33-22:CRC_AV3322_28 viscosity at 373.15 K, `known` and `reconstructed_route`.
  * These are not exceptions. The API now follows the route, and the designated Refutas row of that fuel and temperature is `non_finite` in both the frozen and the new file. Label: `PRE_api_route`.

### Five largest movements

Ranked by |new - old| divided by the R-BASE-2 tolerance: flash 0.05 C, ln viscosity 0.005, density 0.5 kg/m3. The ranking covers bias and MAE of scoreboard cells, one entry per cell.

| # | cell | metric | old -> new | x tol | cause |
|---|---|---|---|---|---|
| 1 | flash T2 fuellib, api_as_shipped | bias | -24.79 -> 2.643 C | 549 | PRE_api_route+D3_overlay |
| 2 | flash T2 ajftd_gcxgc, api_as_shipped | bias | -23.42 -> 0.0019 C | 468 | PRE_api_route+D3_overlay |
| 3 | flash T2 pooled, api_as_shipped | bias | -21.66 -> 0.0961 C (X = -0.692) | 435 | PRE_api_route+D3_overlay |
| 4 | flash T2 av3322, api_as_shipped | bias | -19.74 -> -0.112 C | 393 | PRE_api_route+D3_overlay |
| 5 | ln viscosity T2 av3322 233.15 K, api_as_shipped reconstructed_route | bias | -0.604 -> 0.0319 (X = -0.0315) | 127 | PRE_api_route+D3_overlay |

The frozen API flash rows ran `gbm_mix` and the frozen API viscosity and density rows ran `gn_fallback`: that was `predict()` before delta_api_route moved it onto the route. Those movements are mostly PRE. On the designated route alone (D3 only), the largest movements are:
* flash `family_other`: MAE 3.791 -> 2.570 C, 24 x tol;
* flash `family_D56`: bias -1.065 -> 0.090 C, 23 x tol;
* flash `av3322`: bias -1.062 -> -0.028 C, 21 x tol;
* flash `pooled`: bias -0.734 -> 0.059 C, 16 x tol;
* density `av3322` 313.15 K: bias -2.84 -> -0.70 kg/m3, 4.3 x tol.

### Readings from the new scoreboard (R-BASE-2/3/4 style; reported only, no gate is moved)

**R-BASE-2** (rows of the re-run `SCOREBOARD.md` against their committed counterparts):

| cell | frozen BASE | re-run | committed counterpart | within tol (re-run) |
|---|---|---|---|---|
| flash T2 AJFTD GCxGC, all 37 rows: bias | -0.426 | 0.141 | -0.426 (G1) | no |
| flash T2 AJFTD GCxGC, all 37 rows: MAE | 3.54 | 3.09 | 3.54 (G1) | no |
| flash T2 AJFTD GCxGC, headline one-per-group (31): bias | -0.634 | -0.152 | -0.426 (G1) | no (already no in the frozen BASE: group rule) |
| ln viscosity AV-33-22 233.15 K: bias | -0.637 | -0.574 | -0.637 (C1) | no |
| ln viscosity AV-33-22 313.15 K: bias | 0.0656 | 0.109 | 0.066 (C1) | no |
| flash CRC 647 r3 dev (45): bias | 0.436 | 0.436 | 0.436 (L1) | yes |
| freezing T2 pooled F1: MAE | 3.46 | 3.46 | 3.46 (S5) | yes |

The four T2 rows move outside tolerance because of decision 3 (`D3_overlay:by_code_path`). The re-run's own "cause" text for the viscosity rows ("C1 used S3's mapping route") is generic: it is printed whenever the difference exceeds 0.005. CRC 647 and freezing are unchanged.

**R-BASE-3** (designated model vs the leave-one-group-out constant):
* The per-property table is unchanged except for one cell: flash_point_C "worse" goes from `T2_measured` to `-`. The cause is the `api_as_shipped` designated rows (pooled, ajftd_gcxgc, av3322), which go from `worse_than_trivial` to `beats_trivial`. That is PRE_api_route. The overlay moves only the size of the errors.
* Designated-route verdict changes:
  * flash ajftd_gcxgc (31): `no_skill_shown` -> `beats_trivial` (D3).
  * ln viscosity refutas_raw, av3322 at 253.15 K: `worse_than_trivial` -> `no_skill_shown`.
  * ln viscosity refutas_raw at 313/333/353 K (av3322 and pooled): `beats_trivial` -> `no_skill_shown`. The warm bias grows, for example av3322 313.15 K goes from +0.066 to +0.109 ln (D3).
* Otherwise unchanged: density beats in every set; the cold-viscosity bins stay `worse_than_trivial`; heat `ga_comp_H` stays `no_skill_shown`.

**R-BASE-4** (largest unsafe-side rate, designated rows): unchanged.
* Flash: 0.25 on CRC 647 `fold_b` (12 groups).
* Viscosity: 1.0 on T2 ajftd_gcxgc 233.15 K, refutas_raw (10 groups).
* Within T2, unsafe rates drop or stay the same under the overlay:
  * viscosity av3322 233.15 K: 0.914 -> 0.886;
  * pooled 253.15 K: 0.092 -> 0.031;
  * ajftd_gcxgc 253.15 K: 0.179 -> 0.071;
  * flash T2 pooled: 0.057 (unchanged).

Changed designated cells (scoreboard rows with `variant = designated`; T_K `nan` = aggregate bin; flash and viscosity errors in C and ln; density in kg/m3; heat in MJ/kg):

| property | set | T_K | model | n | bias old -> new | MAE old -> new | unsafe old -> new | verdict old -> new | cause |
|---|---|---|---|---|---|---|---|---|---|
| density_kg_m3 | ajftd_gcxgc | 233.15 | ideal_volume | 3 | -7.495 -> -6.35 | 7.495 -> 6.35 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc | 253.15 | ideal_volume | 3 | -7.726 -> -6.559 | 7.726 -> 6.559 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc | 273.15 | ideal_volume | 3 | -8.13 -> -6.941 | 8.13 -> 6.941 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc | 288.15 | ideal_volume | 27 | -5.952 -> -5.092 | 8.15 -> 7.295 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc | 293.15 | ideal_volume | 3 | -8.706 -> -7.496 | 8.706 -> 7.496 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc | 313.15 | ideal_volume | 14 | -6.343 -> -5.774 | 7.228 -> 6.66 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc | 363.15 | ideal_volume | 7 | -5.063 -> -4.866 | 6.157 -> 6.079 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc_T_263_283K |  | ideal_volume | 3 | -8.13 -> -6.941 | 8.13 -> 6.941 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc_T_ge_283K |  | ideal_volume | 32 | -5.305 -> -4.565 | 7.546 -> 6.811 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | ajftd_gcxgc_T_lt_263K |  | ideal_volume | 3 | -7.495 -> -6.35 | 7.495 -> 6.35 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 233.15 | ideal_volume | 35 | -0.8739 -> 1.18 | 1.848 -> 1.986 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 253.15 | ideal_volume | 35 | -1.289 -> 0.7868 | 1.925 -> 1.822 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 273.15 | ideal_volume | 35 | -1.867 -> 0.2304 | 2.22 -> 1.673 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 288.15 | ideal_volume | 35 | -2.217 -> -0.1034 | 2.47 -> 1.549 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 293.15 | ideal_volume | 35 | -2.369 -> -0.2498 | 2.584 -> 1.588 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 313.15 | ideal_volume | 35 | -2.84 -> -0.6999 | 2.962 -> 1.673 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 333.15 | ideal_volume | 35 | -2.955 -> -0.7947 | 3.038 -> 1.667 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 353.15 | ideal_volume | 35 | -3.314 -> -1.134 | 3.346 -> 1.766 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322 | 373.15 | ideal_volume | 35 | -2.747 -> -0.5483 | 2.882 -> 1.568 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322_T_263_283K |  | ideal_volume | 35 | -1.867 -> 0.2304 | 2.22 -> 1.673 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322_T_ge_283K |  | ideal_volume | 35 | -2.217 -> -0.1034 | 2.47 -> 1.549 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | av3322_T_lt_263K |  | ideal_volume | 35 | -0.8739 -> 1.18 | 1.848 -> 1.986 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib | 233.15 | ideal_volume | 2 | -7.739 -> -6.934 | 7.739 -> 6.934 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib | 253.15 | ideal_volume | 2 | -8.722 -> -7.909 | 8.722 -> 7.909 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib | 273.15 | ideal_volume | 3 | -9.565 -> -8.108 | 9.565 -> 8.108 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib | 293.15 | ideal_volume | 3 | -10.57 -> -9.085 | 10.57 -> 9.085 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib | 313.15 | ideal_volume | 3 | -11.39 -> -9.877 | 11.39 -> 9.877 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib_T_263_283K |  | ideal_volume | 4 | -7.212 -> -6.12 | 7.212 -> 6.12 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib_T_ge_283K |  | ideal_volume | 4 | -8.175 -> -7.063 | 8.175 -> 7.063 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | fuellib_T_lt_263K |  | ideal_volume | 2 | -7.739 -> -6.934 | 7.739 -> 6.934 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 233.15 | ideal_volume | 39 | -1.536 -> 0.4021 | 2.41 -> 2.438 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 253.15 | ideal_volume | 39 | -1.97 -> -0.01045 | 2.541 -> 2.351 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 273.15 | ideal_volume | 40 | -2.729 -> -0.7297 | 3.038 -> 2.395 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 288.15 | ideal_volume | 63 | -3.798 -> -2.256 | 4.881 -> 4.002 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 293.15 | ideal_volume | 40 | -3.273 -> -1.251 | 3.46 -> 2.422 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 313.15 | ideal_volume | 50 | -4.072 -> -2.41 | 4.406 -> 3.339 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 333.15 | ideal_volume | 35 | -2.955 -> -0.7947 | 3.038 -> 1.667 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 353.15 | ideal_volume | 35 | -3.314 -> -1.134 | 3.346 -> 1.766 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 363.15 | ideal_volume | 7 | -5.063 -> -4.866 | 6.157 -> 6.079 |  ->  | n_too_small | D3_overlay (by_code_path) |
| density_kg_m3 | pooled | 373.15 | ideal_volume | 35 | -2.747 -> -0.5483 | 2.882 -> 1.568 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled_T_263_283K |  | ideal_volume | 41 | -2.667 -> -0.7157 | 2.967 -> 2.341 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled_T_ge_283K |  | ideal_volume | 69 | -3.802 -> -2.378 | 4.97 -> 4.153 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| density_kg_m3 | pooled_T_lt_263K |  | ideal_volume | 39 | -1.536 -> 0.4021 | 2.41 -> 2.438 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| flash_point_C | ajftd_gcxgc |  | api_as_shipped | 31 | -23.42 -> 0.001942 | 24 -> 2.925 | 0.03226 -> 0.09677 | worse_than_trivial -> beats_trivial | PRE_api_route + PRE_api_route+D3_overlay + PRE_flash_margin_q2 (by_recompute) |
| flash_point_C | av3322 |  | api_as_shipped | 35 | -19.74 -> -0.1115 | 19.74 -> 1.239 | 0 -> 0 | worse_than_trivial -> beats_trivial | PRE_api_route + PRE_api_route+D3_overlay + PRE_flash_margin_q2 (by_recompute) |
| flash_point_C | fuellib |  | api_as_shipped | 4 | -24.79 -> 2.643 | 24.79 -> 4.784 | 0 -> 0.25 | n_too_small | PRE_api_route + PRE_api_route+D3_overlay + PRE_flash_margin_q2 (by_recompute) |
| flash_point_C | pooled |  | api_as_shipped | 70 | -21.66 -> 0.09613 | 21.92 -> 2.188 | 0.01429 -> 0.05714 | worse_than_trivial -> beats_trivial | PRE_api_route + PRE_api_route+D3_overlay + PRE_flash_margin_q2 (by_recompute) |
| flash_point_C | ajftd_gcxgc |  | ideal_le_chatelier | 31 | -0.6343 -> -0.152 | 3.381 -> 2.995 | 0.09677 -> 0.09677 | no_skill_shown -> beats_trivial | D3_overlay (by_code_path) |
| flash_point_C | av3322 |  | ideal_le_chatelier | 35 | -1.062 -> -0.02755 | 1.499 -> 1.198 | 0 -> 0 | beats_trivial | D3_overlay (by_code_path) |
| flash_point_C | family_D56 |  | ideal_le_chatelier | 28 | -1.065 -> 0.09036 | 1.555 -> 1.286 | 0 -> 0 | beats_trivial | D3_overlay (by_code_path) |
| flash_point_C | family_D93_or_IP170 |  | ideal_le_chatelier | 32 | -0.3867 -> -0.04463 | 3.096 -> 2.876 | 0.09375 -> 0.09375 | beats_trivial | D3_overlay (by_code_path) |
| flash_point_C | family_other |  | ideal_le_chatelier | 10 | -0.9214 -> 0.2996 | 3.791 -> 2.57 | 0.1 -> 0.1 | no_skill_shown | D3_overlay (by_code_path) |
| flash_point_C | fuellib |  | ideal_le_chatelier | 4 | 1.359 -> 2.443 | 5.815 -> 4.731 | 0.25 -> 0.25 | n_too_small | D3_overlay (by_code_path) |
| flash_point_C | pooled |  | ideal_le_chatelier | 70 | -0.7344 -> 0.05854 | 2.579 -> 2.196 | 0.05714 -> 0.05714 | beats_trivial | D3_overlay (by_code_path) |
| heat_of_combustion_net_MJ_kg | ajftd_gcxgc_measured |  | ga_comp_H | 28 | 0.3481 -> 0.3576 | 0.3501 -> 0.3596 |  ->  | no_skill_shown | D3_overlay (by_code_path) |
| heat_of_combustion_net_MJ_kg | av3322 |  | ga_comp_H | 35 | 0.09376 -> 0.1314 | 0.09708 -> 0.1314 |  ->  | no_skill_shown | D3_overlay (by_code_path) |
| heat_of_combustion_net_MJ_kg | measured_pooled |  | ga_comp_H | 28 | 0.3481 -> 0.3576 | 0.3501 -> 0.3596 |  ->  | no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 233.15 | refutas_raw | 10 | -0.6979 -> -0.6866 | 0.6979 -> 0.6866 | 1 -> 1 | worse_than_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 253.15 | refutas_raw | 28 | -0.2911 -> -0.2727 | 0.3442 -> 0.3258 | 0.1786 -> 0.07143 | no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 273.15 | refutas_raw | 3 | -0.07527 -> -0.04211 | 0.07527 -> 0.04211 | 0 -> 0 | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 288.15 | refutas_raw | 3 | -0.00338 -> 0.02758 | 0.009087 -> 0.02758 |  ->  | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 293.15 | refutas_raw | 4 | 0.002334 -> 0.04123 | 0.02064 -> 0.04123 |  ->  | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 313.15 | refutas_raw | 16 | -0.006955 -> -0.000297 | 0.05411 -> 0.05938 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc | 363.15 | refutas_raw | 8 | -0.03667 -> -0.03405 | 0.03837 -> 0.03579 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc_T_263_283K |  | refutas_raw | 3 | -0.07527 -> -0.04211 | 0.07527 -> 0.04211 | 0 -> 0 | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc_T_ge_283K |  | refutas_raw | 17 | -0.01802 -> -0.00696 | 0.04412 -> 0.04921 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | ajftd_gcxgc_T_lt_263K |  | refutas_raw | 28 | -0.4051 -> -0.3864 | 0.4581 -> 0.4394 | 0.4643 -> 0.3571 | worse_than_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 233.15 | refutas_raw | 35 | -0.6373 -> -0.574 | 0.6373 -> 0.574 | 0.9143 -> 0.8857 | worse_than_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 253.15 | refutas_raw | 35 | -0.2643 -> -0.2066 | 0.2643 -> 0.2066 | 0.02857 -> 0 | worse_than_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 273.15 | refutas_raw | 35 | -0.05539 -> -0.00283 | 0.06196 -> 0.03185 | 0 -> 0 | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 288.15 | refutas_raw | 35 | 0.01912 -> 0.06802 | 0.03095 -> 0.06802 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 293.15 | refutas_raw | 35 | 0.03549 -> 0.08318 | 0.04009 -> 0.08318 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 313.15 | refutas_raw | 35 | 0.06564 -> 0.1086 | 0.06564 -> 0.1086 |  ->  | beats_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 333.15 | refutas_raw | 35 | 0.05817 -> 0.09657 | 0.05817 -> 0.09657 |  ->  | beats_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 353.15 | refutas_raw | 35 | 0.035 -> 0.06886 | 0.03525 -> 0.06886 |  ->  | beats_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322 | 373.15 | refutas_raw | 34 | -0.00697 -> 0.02197 | 0.01587 -> 0.02225 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322_T_263_283K |  | refutas_raw | 35 | -0.05539 -> -0.00283 | 0.06196 -> 0.03185 | 0 -> 0 | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322_T_ge_283K |  | refutas_raw | 35 | 0.01912 -> 0.06802 | 0.03095 -> 0.06802 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | av3322_T_lt_263K |  | refutas_raw | 35 | -0.6373 -> -0.574 | 0.6373 -> 0.574 | 0.9143 -> 0.8857 | worse_than_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | fuellib | 233.15 | refutas_raw | 4 | -0.6645 -> -0.6284 | 0.6645 -> 0.6284 | 1 -> 1 | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | fuellib | 253.15 | refutas_raw | 4 | -0.3123 -> -0.2793 | 0.3123 -> 0.2793 | 0.25 -> 0 | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | fuellib | 313.15 | refutas_raw | 4 | 0.0325 -> 0.0573 | 0.0325 -> 0.0573 |  ->  | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | fuellib | 373.15 | refutas_raw | 3 | -0.01966 -> 0.003418 | 0.01966 -> 0.01437 |  ->  | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | fuellib_T_ge_283K |  | refutas_raw | 4 | 0.03188 -> 0.05668 | 0.03188 -> 0.05668 |  ->  | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | fuellib_T_lt_263K |  | refutas_raw | 4 | -0.6645 -> -0.6284 | 0.6645 -> 0.6284 | 1 -> 1 | n_too_small | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 233.15 | refutas_raw | 48 | -0.6503 -> -0.5995 | 0.6503 -> 0.5995 | 0.9375 -> 0.9167 | worse_than_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 253.15 | refutas_raw | 65 | -0.2754 -> -0.2362 | 0.2983 -> 0.259 | 0.09231 -> 0.03077 | worse_than_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 273.15 | refutas_raw | 38 | -0.05696 -> -0.005931 | 0.06301 -> 0.03266 | 0 -> 0 | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 288.15 | refutas_raw | 38 | 0.01735 -> 0.06483 | 0.02922 -> 0.06483 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 293.15 | refutas_raw | 39 | 0.03209 -> 0.07888 | 0.03809 -> 0.07888 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 313.15 | refutas_raw | 54 | 0.04203 -> 0.07328 | 0.06012 -> 0.09097 |  ->  | beats_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 333.15 | refutas_raw | 35 | 0.05817 -> 0.09657 | 0.05817 -> 0.09657 |  ->  | beats_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 353.15 | refutas_raw | 35 | 0.035 -> 0.06886 | 0.03525 -> 0.06886 |  ->  | beats_trivial -> no_skill_shown | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 363.15 | refutas_raw | 8 | -0.03667 -> -0.03405 | 0.03837 -> 0.03579 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled | 373.15 | refutas_raw | 37 | -0.007999 -> 0.02046 | 0.01618 -> 0.02161 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled_T_263_283K |  | refutas_raw | 38 | -0.05696 -> -0.005931 | 0.06301 -> 0.03266 | 0 -> 0 | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled_T_ge_283K |  | refutas_raw | 54 | 0.009763 -> 0.04496 | 0.03563 -> 0.06264 |  ->  | beats_trivial | D3_overlay (by_code_path) |
| viscosity_kin_cSt | pooled_T_lt_263K |  | refutas_raw | 65 | -0.547 -> -0.5045 | 0.5699 -> 0.5273 | 0.7385 -> 0.6923 | worse_than_trivial | D3_overlay (by_code_path) |

### Tests

**The 4 tests that failed in the `--full` pytest stage.** Every frozen value pin is unchanged.

| test | why it failed | now |
|---|---|---|
| `test_base_cells.py::test_cells_t2_flash` | reads the regenerated `rows_t2_flash.csv` / `cells_t2_flash.csv` (headline bias -0.152, pin -0.634) | passes on the restored frozen files; no edit |
| `test_base_cells.py::test_cells_t2_bulk` | reads the regenerated `viscosity_error_vs_T.csv` (AV-33-22 233.15 K bias -0.574, pin -0.637) | passes on the restored frozen files; no edit |
| `test_base_scoreboard.py::test_designated_columns_untouched_by_variants` | reads the regenerated `scoreboard.csv` (flash AJFTD bias -0.152, pin -0.634) | passes on the restored frozen files; no edit |
| `test_base_composition_route.py::test_reproduces_s3_corrected_M_predictions` | recomputes live and compares with the frozen S3 file. Under the overlay, 125 AJFTD rows differ, for example AJFTD:2 density 791.993 -> 793.088 kg/m3 and viscosity at 233.2 K 4.6073 -> 4.7646 cSt | **legacy-pointer fixture**; passes |

**Tests that got a legacy-pointer fixture.** No expected value was edited, and `DESIGNATED_FREEZING_MODEL` was not needed.

* `test_base_composition_route.py::test_reproduces_s3_corrected_M_predictions`
  * Uses the new module-scoped fixture `legacy_parts_ctx`.
  * The fixture runs inside `pytest.MonkeyPatch.context()`: it sets `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0`, calls `flash_unifac._default_registry.cache_clear()`, then builds `SubReps` and `RouteCtx`. Both read the registry at construction (`has_extra`, `flash_reg`). It clears the cache again on exit.
  * The other tests in the file keep the default `parts` / `route_ctx`.
* `test_predict_route.py::test_heat_matches_bulk_rows` (AJFTD:2 `heat_net_MJ_kg`: live 43.5287022 under the overlay vs frozen 43.508441)
  * Uses the new fixture `legacy_flash_registry`.
  * The fixture calls `monkeypatch.setenv(..., "0")` and clears `_default_registry`, `predict_route.get_ctx` and `predict_route._overlay_ctx` before and after. `_overlay_ctx` holds `get_ctx().reps`.
  * The re-run's `rows_t2_bulk.csv` holds 43.528702 for AJFTD:2, the same value the live API gives under the adopted default.

**New test** `tests/test_delta_adopt_rescore.py` (3 tests, reads committed files only):
* every scoreboard cell that differs between the frozen BASE and `rerun/` is listed in `drift.csv`, and nothing else is;
* no `D1`/`D4`/`D5`/`D7` cause appears;
* freezing and all non-T2 cells are bit-identical (E2);
* no `error` rows appear in `rerun/`.

**Run on the restored tree:** the 4 tests above + `test_heat_matches_bulk_rows` + the 3 new tests: 8 passed.

**Side effect of the pytest stage of `reproduce_base.py`.** `test_base_delta_representatives` and `test_base_delta_vapor_pressure` regenerate their `delta_*` files live, now under the overlay. They pass, but they rewrite 7 tracked files. Those files were restored with `git checkout` after each run. The `--full` versions (11 files) are kept in `rerun/delta_files/` for inspection; they are not adopted.

### Restore and final check

* **How the files were restored.** The frozen BASE files were restored by **copying back `BASE_frozen`**, not with `git checkout -- <paths>`:
  * `.gitattributes` has `*.csv text eol=lf`, so a checkout writes LF.
  * The manifest hashed the CRLF files that pandas writes on Windows.
  * After `git checkout`, `scoreboard.csv` hashed `a15d...` instead of the manifest's `941b...`.
* **Byte-level check.** All 23 frozen files (22 in the manifest, plus `data_thermoml_sle/README.md`) are byte-identical to `BASE_frozen`.
* **What git shows.** `git status` still lists 11 of them as ` M`. That is line endings only: `git diff` is empty and `git hash-object` equals the index blob (HEAD content).
* `.venv/Scripts/python.exe scripts/v5/base/reproduce_base.py` (default, `--from-rows`), run after the fixtures were added:
  * pytest stage: 47 passed, 1 skipped;
  * `[MANIFEST] DRIFT: 0 changed, 0 missing, 1 new (23 files hashed)`. The one `NEW` is the declared, pre-existing `results/v5/BASE/data_thermoml_sle/README.md`, which is why the script returns 1.
  * Its `--from-rows` pass rewrote `rows_t2_bulk.csv`, cells and scoreboard byte-identically, as the cmp check after the run confirms.
* No untracked file was left in `results/v5/BASE` outside `delta_adopt_rescore/`. `results/v5/BASE/delta_adopt_designs` was not touched.

### Status

**The frozen BASE (`results/v5/BASE/**`, non-`delta_*`) and `SHA256SUMS_v5_base.txt` are NOT overwritten.** `--update-manifest` was never passed.

The re-scored numbers live only in `rerun/` and `drift.csv`. Adopting the re-run as the new frozen BASE is question **U4** for the user. If it is adopted, these updates go with it:
* the 3 frozen-value pins in `test_base_cells.py` and `test_base_scoreboard.py`;
* the 2 legacy fixtures;
* the parent README's R-BASE-2 statement.

None of these is done here.
