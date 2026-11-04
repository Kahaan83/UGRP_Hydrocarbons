# WP0 — Housekeeping & source check (status)

Date: 2026-10-30 · seed fixed: 20261030 · model seeds: [11, 23, 37, 41, 59]

## What exists and was verified
- Repo: `C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\fuelblend_data` (no git repo yet).
- All 18 test suites pass (`python tests/test_*.py`, pytest installed in venv for `test_nist_density`).
- `audit.py --bundle`: PASS — 0 FAIL, 0 MISSING, 7 WARN (OneDrive location, 1 macOS junk file, 2 empty/read-error files, Antoine coverage, fake DTIC PDFs).
- ThermoML re-extract clean: 46,747 rows, 100 % with T_K (fixed parser).
- Unified build: 71,129 raw rows → 47,531 dedup records, 700 conflicts flagged, 1,038 compounds.
- Leenhouts benchmark + splits downloaded: `data/interim/leenhouts_benchmark/{datasets,splits}`.
- Baselines previously verified: Refutas 4.8 % MAPE on binary viscosity; Le Chatelier flash model runs on 147 predictable records (unified report said 99 — needs reconciliation).

## Data facts (measurements_dedup.csv, 47,531 × 26)
- kind: binary 31,199 · pure 9,665 · real_fuel 4,869 · multi 1,798.
- properties: viscosity_dyn_mPas 18,932 · density_kg_m3 15,907 · vapor_pressure_kPa 6,500 · flash_point_C 989 · dcn 426 · mon 373 · viscosity_kin_cSt 361 (+ others).
- sources: thermoml 26,484 · nist_tde 5,499 · leenhouts 3,649 · crc647 3,511 · av2718 869 · chew2024 574 · flashpoint_li 398 · fortin2012 252 · ajftd 234 · pure_freezing_nist 25 · saldana11 13 · av1109 3.
- quality: 700 conflicts, ~290 rows with n_outlier > 0, 6,020 rows with ≥2 sources.
- sha256 (dedup) = `74df03826348dee88f4252281a5eaac33bebc3c14226b67308fdd6b0073078f9`; measurements `e78b7a0d…`; compounds `8333ea70…`.
- 112 distinct real fuel_id values present.

## Harness artifacts (WP2, built on WP0)
- `data/processed/model_table.csv` — 47,531 rows, HC-only, quality_weight = 1/sqrt(n_obs), source tiers A/B/C, domain tags.
- `src/fuelblend/splits/split_generator.py` — S1 28,729 / S2 8,050 / S3 826 (molecule holdouts, seed 20261030); locked real-fuel rows and locked multi system_ids excluded from S1/S2 (locks generated first).
- `src/fuelblend/splits/locked_sets.py` → `splits/locked.json`: 24 locked real fuel_ids + 12 locked multi system_ids (20 % per source bucket), dev/final halves.
- `src/fuelblend/eval/metrics.py` — MAE/RMSE/medianAPE/frac_within, flash over-prediction rate, cluster bootstrap (by system_id, 2,000 resamples).
- Reports: `DATA_CARD.md`, `GATES.md` (gate classification per handoff 17.4), `NOISE_FLOOR.md`, `LOCKED_ACCESS.log`.
- Leakage test result: LOCKED-SET LEAKAGE TEST PASSED (after excluding locked/real rows from S1/S2).

## Gate classification (conflict-free rows, dedup per system_id)
- G2: viscosity_dyn_mPas (379 bt), density_kg_m3 (235 bt / 94 real), smoke_point_mm (18 real).
- G3: viscosity_kin_cSt, flash_point_C, pour_point_C, mon, dcn, net_heat_content_btu_lb.

## Open items (needs approval or waiting)
1. ANALYSIS_PLAN.md approval (metrics, locks, thresholds, flash asymmetric loss weights {1,3,5,10}).
2. T00 SOURCE_CHECK against ASTM quotes (needs the exact ASTM precision statements).
3. Leenhouts baseline reproduction (paper + splits.zip) — subagent attempt failed on harness/network error (`getaddrinfo ENOTFOUND`); must be re-run.
4. `git init` + move project out of OneDrive (≈800 files sync badly; one file write was observed to land 0-byte/vanish).
5. Reconcile the 147 vs 99 flash-predictable-records discrepancy between baselines.py and UNIFIED_REPORT.md.

## Honest limitations
- No numbers here are invented; all counts come from the unified CSVs on 2026-10-30.
- Noise-floor table is an upper-bound proxy (observed spread within bins), not a calibration.
- Locked AJFTD/fuel_ids are excluded from training splits but the real-fuel G3 gate means no modelling target there — baseline-only.

## Decisions made
- Gates follow the strict handoff 17.4 rule verbatim (G1 needs ≥1,000 binary/ternary AND ≥100 multi AND ≥20 real rows; nothing currently qualifies as G1).
- Locked sets generated before splitting; real_fuel rows never enter S1/S2.
- Failed/absent items listed above as TODO, not silently skipped.
