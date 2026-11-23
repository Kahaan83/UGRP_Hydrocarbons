# v10 predeclaration (fixed before any v10 value is read)

Date: 2026-11-23. Plan: `docs/superpowers/plans/2026-11-23-v10-pdf-transcription.md`.

## Decisions
- Inputs are PDFs on disk in `data/raw/v10_pdfs/` (never committed). No web fetching.
- Gemini 3.8 Flash drafts; Sonnet checks every row that can move a reading against the rendered page image.
- v10-new means a group with no row in `data/unified_v9/measurements_dedup.csv`.
- Lock: `md5(group_id) mod 4 == 0` over v10-new groups, written once to `results/v10/lock.json`, scored once, 8-fuel floor.
- Commit per task, one push at the end to `origin feat/external-datasets-ingest`. Explicit `git add` paths only.

## Readings
| reading | pass rule |
|---|---|
| R10-1 real-fuel freezing | Model: v5 designated UNIFAC Schroeder-van Laar route (`unifac_desig`, as in `scripts/v8/g2_freezing_check.py`). Constant: leave-one-group-out mean freezing point of dev fuels. Skill shown only if the 90 % paired group bootstrap CI (5000 resamples, seed 42) of (constant MAE minus model MAE) excludes 0 on the v10 lock, with at least 8 v10-new fuels with a closed carbon-split composition. Below the floor: descriptive, status "no skill shown". A dev-only gain cannot reopen it. |
| R10-2 real-fuel cold viscosity at 253.15 K and 233.15 K | Model: v5 Refutas route on carbon-split composition. Constant: leave-one-group-out mean ln(cSt) of dev fuels at that temperature. Same bootstrap, rule and floor, per temperature. |
| R10-3 agreement with existing rows | Descriptive: count and size of differences beyond tolerance (freezing 1 degC, flash 2 degC, else 5 %) for fuels already in v5 or v9. |
| R10-4 design candidates with measured constants | Descriptive: candidates with measured LFL, viscosity, density and a converged Antoine refit; and the subset whose refit T range overlaps 313 to 373 K. No design re-run in v10. |

## Phase 0 observations (headers and page counts only)
See `results/v10/pdf_register.csv`. AIAA 2017-0146 not provided. `sources.txt` absent, titles blank.
