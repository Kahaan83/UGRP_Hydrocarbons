# delta_api_route: declared readings (written before any result exists)

## Purpose

Check that the re-routed prediction API (decision D1 to D4 of the v5 API-route plan) reproduces the scoring route that BASE scored, instead of the old shipped API path. This README is committed before any delta result is generated. It declares what will be read and what counts as agreement; it reports no results.

## Baseline

`api_as_shipped` on the T2 flash set: bias -21.66 C, MAE 21.92 C, over 70 groups (`results/v5/BASE/cells_t2_flash.csv`, pooled row, `api_as_shipped`/`designated`). Values were verified against the CSV before writing (bias -21.656749, MAE 21.916079, n_groups 70, n_rows 70). `rows_t2_flash.csv` has 79 rows.

## Declared readings

- **R-API-1.** On all 79 T2 rows, the API's designated flash equals the route's `pred_lfl_overlay`, and `variants['v4_designated']` equals `pred_designated`, both to 1e-6 C (columns of `results/v5/BASE/rows_t2_flash.csv`).
- **R-API-2.** API kinematic viscosity and density equal the route's values on all T2 rows at the 9 BASE temperatures, to 1e-9 relative.
- **R-API-3.** Flash, heat and aromatics are bit-identical across T_K in {233.15, 288.15, 313.15}.
- **R-API-4.** Both D4 v5pool recipes re-predict through `predict()` within 0.01 C for flash, 1e-6 relative for viscosity, and exactly for heat. Density and freezing differences are reported with their cause.

## Grid SHA-256 and disclosure

`data/cache/pure_grid.npz` is an unpinned, gitignored binary built 2026-11-11. Its SHA-256 is pinned here and in `tests/test_delta_api_route.py`:

`GRID_SHA256 = b931e969926698dff52d69ff51dd99a990c24ff12eca545aa98e970afd3e7daf`

As reported by the BASE work, the grid's values differ from today's `_fetch_via_api` by up to 3.25 in ln eta at 230 K. This is disclosed, not fixed here. The grid holds 1042 compounds, 114 with a full finite grid, starting at T = 230 K.

## Out of scope

Grid rebuild; U1 promotions (K1, F1, LFL overlay); distillation.
