# UGRP Hydrocarbons

Predict fuel blend properties (Jet A/A-1, diesel, RP-1) from composition, then inverse-design blends. Authors: Kahaan Shah, Namyaa Parmar, Arijit Ganguli, Jayendra Bhalodiya, A. S. K. Sinha (Ahmedabad University).
Docs (read only when needed): `docs/PROJECT_HANDOFF.md`, `docs/GUIDE.md`, `docs/final_numbers_v5.md` (source of truth for numbers).

## Model (v5)
- Viscosity: Refutas blending + optional LightGBM pair residual.
- Flash point: Le Chatelier + Raoult + Antoine + LFL overlay.
- Density: ideal volume mixing + optional LightGBM residual.
- Freezing point: Schroeder-van Laar + UNIFAC.
- Uncertainty: one-sided 90% conformal margins from held-out lab systems.
- Public API: `src/fuelblend/api/predict_route.py` (`predict_fuel()`).

## Commands (Windows, PowerShell, Python 3.12 in `.venv`)
- Tests: `.venv\Scripts\python.exe -m pytest -q` (run one file/test, full suite only at the end). Baseline (2026-11-23, full run after v11c hybrid): 911 passed, 2 skipped, 1 failed. `test_committed_splits_match_rebuild` fails (known, AJFTD input hash); `test_fit_ensemble_logger_none_bitwise` is load-sensitive (passed in this run, failed once under load).
- Setup: AJFTD tools need `beautifulsoup4` and `playwright`; after install run `python -m playwright install chromium`.
- Keep the project out of OneDrive if possible (~800 small files, ~600 MB raw data sync badly; deep zip paths hit the Windows 260-char limit).
- v5 baseline check: `python scripts/v5/base/reproduce_base.py`
- v4 pipeline: `python scripts/v4/reproduce_all.py` (do not rerun unless asked)

## BASE freeze (v5)
- Frozen BASE: `results/v5/BASE/`, manifest `SHA256SUMS_v5_base.txt`. Adopted re-score (U4, 2026-11-15) is adoption commit `91f2380`; the pre-adoption BASE is recoverable at its parent `09e1366` (hashes after the history rewrite, see `results/v5/BASE/README.md`).
- Do NOT `git checkout` BASE CSVs or normalise their line endings: the manifest hashes exact bytes, and `.gitattributes` marks `results/v5/BASE/**` as `-text`. Verify with `reproduce_base.py`; never regenerate BASE numbers.

## Data
- Canonical units: cSt, mPa.s, kg/m3, degC, kPa, MJ/kg, mm. Temperatures in K inside equations.
- Bases: `w` mass (viscosity, density inputs), `x` mole (flash, VLE), `v` volume (pour point).
- Antoine: log10(P/kPa) = A - B/(T_K + C). NIST bar to kPa: A_kPa = A_bar + 2.
- Always `GroupKFold` on `system_id` (md5[:12] of the sorted InChIKey set; `fuel:<id>` for real fuels). Train and evaluate from `measurements_dedup.csv`, never raw rows.
- Dedup on (property, system_id, composition, T, P): drop exact repeats (same source+ref+value), take the median, flag (never silently delete) outliers and conflicts. Bins: T to 0.5 K; P <=150 kPa or missing -> 101, else rounded to 100 kPa; x to 3 decimals.
- `measurements.csv` columns: `source` (leenhouts, nist_tde, thermoml, ajftd), `property`, `value`, `T_K`, `P_kPa`, `kind` (pure, binary, multi, real_fuel), `n_comp`, `is_hc`, `system_id`, `comp_inchikeys`, `x_mole`, `w_mass`, provenance (`fuel_id`, `ref`, `raw_value`, `raw_unit`). Dedup adds `vmin`, `vmax`, `q25`, `q75`, `n_obs`, `n_refs`, `n_sources`, `sources`, `n_outlier` (>25 % from median, needs >=3 obs), `conflict`.
- Conflict tolerances (IQR, or range if 2 obs): abs 2 degC flash, 3 degC pour/cloud/DCN/MON, 1 degC freezing; 5 % relative otherwise.
- Hydrocarbon-only (`is_hc`) is the default filter. All unit conversions happen in loaders only. Every new source needs a loader to this schema, a unit test with a fixture mirroring the real file structure, and an audit check.
- Spec limits (real fuel targets): flash 38 degC min (Jet A, D56); viscosity 8 cSt max at -20 degC; density 775-840 kg/m3 (AJFTD spec limits).
- Never commit `ajftd_session.json`, raw AJFTD pages, or licensed data; `audit_bundle.zip` holds no raw datasets by design.
- Production data: `data/unified/`. `data/quarantine/` is formula-derived, never train on it.

## Traps
1. Viscosity log units: NIST ln(Pa.s) vs Leenhouts ln(mPa.s); offset ln(1000).
2. ThermoML: constraint T is in `Constraint > ConstraintID > ConstraintType`; liquid phase only; skip excess properties.
3. Never evaluate un-refit Antoine constants for C10+ alkanes in the 40-100 degC flash window.
4. Pandas 3: no `.loc[k, [cols]] = [...]` with mixed types; assign columns one by one.
5. Skip `__MACOSX/._*` and `~$*` files in parsers.
6. Flash point: over-prediction is the dangerous direction. Report one-sided over-prediction error, not only MAE.
7. Circular labelling: when a result looks too good, check for synthetic data or the blending rule used as the label (`DatasetsFuelProperties.xlsx` is 50,000 synthetic rows; never ground truth).
8. Historic bugs, do not reintroduce:
   - ThermoML: native units are Pa.s and m2/s, so viscosity needs x1e3 (mPa.s) or x1e6 (cSt); constraint T/P live under `Constraint > ConstraintID > ConstraintType` (61 % of T lost otherwise); keep liquid phase only, whitelist property names (no excess viscosity), ignore vapour-phase composition variables, drop P > 1000 kPa.
   - Leenhouts flash unit string is `Celcius` (sic).
   - DTIC returns an HTML consent page saved as a "PDF"; validate PDF headers.
   - `compounds.csv` name with an embedded comma (`1,2,4-trimethylbenzene`) must be quoted.
   - AJFTD scraper: every fuel was named with the site banner, no resume, no pagination cross-check; check "Showing x to y of N entries" against the fetched count. WebBook `page_name` also returned the banner.
   - Antoine refit crashed on falling/flat pressure data (bad starts); keep clipped starts and per-compound `failed:` status.
   - Synthetic-test lesson: a test built on a guessed XML structure passed while real data failed. Keep real-record regression tests (`tests/data/je8006138.xml`).

## Rules
1. Never invent or alter numbers; every figure must trace to `docs/final_numbers_v5.md` or committed results.
2. State "status: no skill shown" where a model does not beat a constant (real-fuel freezing point, cold viscosity).
3. Design certification is model-level only; no lab measurements exist.
4. Keep `data/unified/`, `scripts/v5/`, `src/fuelblend/` unchanged unless the task requires it.
5. Academic, direct tone. No em dashes in manuscripts or decks. No AI attribution in commits or PRs.
6. Keep personal names, emails, and project details out of external services.

## Usage limits
- Never bulk-read `data/`, `results/` CSV/JSON, `.venv/`, `.docx`, `.pdf`. Use Grep, `head`, `Read` with `limit`, or short Python summaries.
- NO subagents, `/code-review max|ultra`, or web fetching unless asked.
- Cap command output (`Select-Object -First 40`); send long logs to `scratch/`.
- Sonnet  and haiku for edits, Opus only for hard decisions, never `fable`.
- Terse replies, no diff recaps, one recommendation instead of option surveys. New session per task.
