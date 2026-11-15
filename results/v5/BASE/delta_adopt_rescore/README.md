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
