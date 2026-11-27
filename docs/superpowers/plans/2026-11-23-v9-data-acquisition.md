# v9: data acquisition for real-fuel freezing, mixtures, pure constants and cold viscosity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (native). No subagents. Three sequential phases with fixed models; see "Model assignment" and run only your phase. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add measured data where v5 to v8 showed data, not models, is the limit: real-fuel freezing with detailed composition, multicomponent and cold mixture systems, measured pure constants for design candidates, and real-fuel cold viscosity; then re-read the affected models on dev and score a new lock once.

**Architecture:** Each source becomes a raw CSV plus `SOURCES.md`, a loader to the project schema, a fixture test and an audit check. New rows go into a new layer `data/unified_v9/` (v5 unified rows plus new sources, re-deduplicated); `data/unified/` and `results/v5/BASE/` are untouched, so v5 reproduces. A deterministic hash rule locks a quarter of the NEW systems and fuels before any model sees them.

**Tech Stack:** Python 3.12 (`.venv`), pandas, RDKit, pdfplumber/pypdf/pypdfium2, requests, beautifulsoup4, playwright (only where a site needs it), existing `fuelblend.data` loaders and `fuelblend.v7.guard`.

**Spec:** this plan; the evidence it argues from: `docs/final_numbers_v5.md` (sections A, F), `results/v8/substitution/READING.md`, `.superpowers/sdd/2026-11-12-v5-api-route-and-design/ingest-common.md` (ingest rules, followed verbatim).

## User decisions (recorded 2026-11-23, before any v9 work)

* Targets: all four (real-fuel freezing, mixture systems, pure constants, real-fuel cold viscosity).
* Licences: open and government sources plus numeric facts transcribed with citation from journal supporting information (Vozka 2018 practice). Raw PDFs and licensed files are never committed.
* Web access: fetching pages and downloading papers, SI and the ThermoML archive is allowed without asking.
* Re-evaluation: lock a share of the new data before looking, re-read models on dev, score the new lock once.
* Git: commit after each task; push to `origin feat/external-datasets-ingest` once, after Task 9.
* Execution: native, no subagents; Gemini 3.8 Flash for gathering and drafting, Sonnet for verification, locks, readings and push (user decision 2026-11-23).

## Model assignment (Gemini 3.8 Flash for gathering and drafting, Sonnet for verification and anything that can spend a lock)

Three sessions, run one after another by the user (no subagents, no parallel sessions). Each phase starts by reading this plan and the previous phase's handoff file.

| phase | model | tasks / steps | ends with |
|---|---|---|---|
| S0 | Sonnet | Task 0 (precondition, predeclaration, lock rule) | commit; `results/v9/HANDOFF_S0.md` |
| G | Gemini 3.8 Flash | Task 1 (inventory); Tasks 2 and 3 Steps 1 to 6 (search, fetch, draft transcription, loaders, tests); Task 4 Steps 1 to 6 (fetch archive, fixture test, loader, extraction); Task 5 Steps 1 to 5 (collect constants, refit script, tests); Task 7 Step 1 (write `audit_v9.py`, do not interpret) | local commits only; `results/v9/HANDOFF_gemini.md` |
| S1 | Sonnet | verify every Gemini transcription and review every Gemini loader/test (rules below); Task 5 Step 6 (run refit); Task 6; Task 7 Steps 2 and 3; Task 8; Task 9 | push |

**Gemini rules (phase G).**
* Writes transcriptions to `data/raw/<folder>/<source>_*_draft.csv`, never to the final file name. Every row gets `page`, `table`, and `confidence in {clear, unsure}`; a value read from a garbled text layer, a rotated table or a scanned page is `unsure`.
* For every source it saves the exact PDF page numbers used and, for scanned or garbled pages, renders them to PNG under `extras/staging/lit/v9/pages/` (gitignored) for the verifier.
* May write loaders, tests, fetch scripts, `inventory.py`, `audit_v9.py`, and fill `data/raw/v9_source_log.csv`. Loader tests may read `_draft.csv` until phase S1 renames them.
* Must not: create or edit `data/unified_v9/`, `results/v9/lock.json`, `readings.py`, `make_lock.py`, anything under `results/v5/`, `data/unified/`, `splits/`, `configs/`, `src/fuelblend/data/model_table.py`, `src/fuelblend/splits/`; run `reproduce_base.py` writes or any model reading; push.
* Commits locally with prefix `data(v9-draft):` or `feat(v9-draft):`.
* `HANDOFF_gemini.md` lists: sources ingested / not found per target, every draft file with row count and number of `unsure` rows, every judgement call, tests added and their status, and anything left half done.

**Sonnet rules (phase S1).**
* Verification of transcriptions: every `unsure` row, and every row of real-fuel freezing and cold-viscosity drafts (small and decisive), is re-read against the page image. For larger tables (more than 200 rows from one table), every `unsure` row plus a random 20 % of `clear` rows (seed 42); if any sampled `clear` row is wrong, that whole table is re-read. Corrections are recorded in `SOURCES.md`. Then the draft is renamed to its final name and loader tests are pointed at it.
* Code review: each Gemini loader is checked against the Global Constraints (unit conversions in loaders only, ThermoML traps, PDF header check, duplicate handling, `qualifier` for censored values, `blend_basis_unconfirmed`). Fixes get a failing test first.
* Sonnet alone does the identity map, `unified_v9`, the lock, all readings, final numbers and the push.
* Commits use the normal prefixes; `results/v9/verification.md` records rows checked, rows corrected and tables re-read per source.

## Unattended execution rules

* Each phase runs start to finish without asking the user.
* A source that cannot be found, is paywalled with no open copy, or has no usable numbers is a reading, not a stop: add a row to `data/raw/v9_source_log.csv` (`source, target, status in {ingested, not_found, paywalled, plot_only, duplicate, no_composition}, note`) and continue.
* A site that refuses scripts (403, "Access Denied", consent page) is tried once, logged, not retried. CAPTCHAs are never attempted.
* Stop only on a broken invariant: `reproduce_base.py` fails, a write into `data/unified/`, `results/v5/BASE/` or `splits/`, or a test that cannot pass without changing a gate. Then commit what is done, write `results/v9/STOPPED.md`, do not push.
* Long downloads and model runs go in the background with resumable outputs; logs go to `scratch/`.
* Search budget: at most 10 search queries and 15 page fetches per candidate source in Tasks 2 and 3, and at most 150 fetches per task. When the budget runs out, the source is logged `not_found` and the run moves on.
* No personal names, emails or project details in any request (AGENTS.md rule 6): generic client name in User-Agent; skip APIs that require an email (Unpaywall, Crossref polite pool).
* Precondition: v8 is closed (pushed, or `results/v8/STOPPED.md` exists) before Task 0 starts.

## Global Constraints

* Never invent or estimate a value; plot-only data is skipped (`plot_only`). Every number is checked against the page image when the text layer is unsure (minus signs, subscripts).
* Every new source: raw CSV in `data/raw/<folder>/`, `SOURCES.md` (citation, doi, licence, pages/tables, transcription method, every judgement call, no quoted paper text), loader `src/fuelblend/data/loaders_<folder>.py`, test `tests/test_loaders_<folder>.py` with a fixture mirroring the real file structure, and an audit check (Task 7).
* Canonical units only after the loader: cSt, mPa.s, kg/m3, degC, kPa. Conversions happen in loaders only.
* `system_id` = md5[:12] of the sorted InChIKey set; real fuels `fuel:<source>:<sample>`. Hydrocarbon-only (`is_hc`) default.
* Dedup on (property, system_id, composition, T, P) with the AGENTS.md bins and tolerances (flash 2 C, freezing 1 C, pour/cloud/DCN/MON 3 C, 5 % otherwise); flag outliers and conflicts, never silently delete.
* ThermoML traps: constraint T/P under `Constraint > ConstraintID > ConstraintType`; liquid phase only; whitelist property names (no excess properties); Pa.s x1e3, m2/s x1e6; drop P > 1000 kPa; regression test on `tests/data/je8006138.xml` stays green.
* Validate PDF headers (`%PDF-`); a consent page saved as PDF is logged `not_found`. Skip `__MACOSX/._*` and `~$*`.
* Never commit `ajftd_session.json`, raw AJFTD pages, PDFs, page images or licensed files. `data/raw/**/*.pdf` stays gitignored.
* Real fuels already in AJFTD, CRC AV-33-22, CRC 647, FuelLib or Vozka 2018 are matched by label (POSF number, sample id) through the identity map of Task 6. Their new rows are kept but rewritten onto the existing canonical `system_id` (so locks, folds and the guard see one fuel), and logged `duplicate` in the source log.
* No changes to `data/unified/`, `scripts/v5/`, `configs/uncertainty_v5.yaml`, `results/v5/BASE/**`, `splits/**`, existing locks.
* Tests: `.venv\Scripts\python.exe -m pytest -q <file>`; full suite only in Task 9 (baseline 741 passed, 11 skipped, plus the known `test_committed_splits_match_rebuild` failure reported in v8).

## Review Focus

1. The same fuel reached through two sources (for example an AFRL POSF fuel also in AJFTD): rewritten onto one canonical `system_id` by the Task 6 identity map, never double-counted in a fold or a lock, and still caught by `splits_locked` / `freezing_holdout`.
2. Censored or qualified values (`<-60`, "below", "pass"): stored with `qualifier`, excluded from MAE and margins unless a reading says otherwise.
3. Blend basis not stated (vol vs mass): `qualifier=blend_basis_unconfirmed`, never silently treated as mass.
4. A new fuel with composition by class only (no carbon-number split): stored, but flagged `no_carbon_split` and excluded from freezing readings (freezing needs the n-paraffin tail).
5. A new system that overlaps an existing v6/v7 lock: `fuelblend.v7.guard.assert_dev_only` must raise before any dev reading uses it.

---

### Task 0: Predeclaration, source log and lock rule

**Files:**
- Create: `docs/v9_predeclaration.md`, `data/raw/v9_source_log.csv` (header only), `src/fuelblend/v9/__init__.py`, `src/fuelblend/v9/lock.py`
- Test: `tests/v9/test_lock.py`

**Interfaces:**
- Produces: `is_v9_locked(group_id: str) -> bool` (True when `int(md5(group_id).hexdigest(), 16) % 4 == 0`); `v9_lock_ids(new_ids: Iterable[str]) -> frozenset[str]`; lock file `results/v9/lock.json` (`{"rule": "md5 mod 4 == 0", "ids": [...], "created": ...}`), written once by Task 8 and refused if it exists.

- [ ] **Step 0: Precondition check.** `git log origin/feat/external-datasets-ingest..HEAD` is empty for v8 commits or `results/v8/STOPPED.md` exists; `git status --short` shows only the known modified `results/v5/BASE/**` CRLF files and untracked files present at plan time. `splits/locked.json` reappears as modified after test runs with a line-ending-only change: if `git diff --ignore-cr-at-eol --quiet splits/locked.json` succeeds, `git restore splits/locked.json` and continue; any content change stops the run. Otherwise write `results/v9/STOPPED.md` ("v8 not closed") and stop.
- [ ] **Step 1: Write `docs/v9_predeclaration.md`:** the User decisions section verbatim; the lock rule; the readings R9-1 to R9-4 with the pass rules in the table below; the statement that only data new in v9 can enter the v9 lock, so no v5/v6/v7 lock is re-read.

| reading | pass rule (statistics: paired group bootstrap over fuels or systems, 90 % CI) |
|---|---|
| R9-1 real-fuel freezing | skill shown only if the CI of the gain (constant MAE minus model MAE) excludes 0 **on the v9 lock**. A dev-only gain cannot reopen "status: no skill shown". If the v9 lock holds fewer than 8 freezing fuels with carbon split (8 = the v5 locked freezing set, `docs/final_numbers_v5.md` section A), the lock reading is descriptive only and the status stays "no skill shown" |
| R9-2 viscosity residual with more data | retrained model vs the v5 `ml_pair_residual`, both out-of-fold on the v5 dev systems: "more data helps" only if the CI of the MAE gain excludes 0; then confirmed on the v9 lock with the same rule |
| R9-3 real-fuel cold viscosity (253.15 K, 233.15 K) | as R9-1, per temperature, same lock floor of 8 fuels |
| R9-4 v5 routes on new mixture systems | descriptive: MAE and bias per property vs the v5 lab figures; no verdict change |
- [ ] **Step 2: Failing tests**

```python
def test_lock_rule_is_deterministic():
    assert is_v9_locked("fuel:ajftd:2") == is_v9_locked("fuel:ajftd:2")

def test_lock_share_near_quarter():
    ids = [f"fuel:x:{i}" for i in range(4000)]
    assert 0.22 < sum(map(is_v9_locked, ids)) / 4000 < 0.28
```

- [ ] **Step 3: Run; expect FAIL. Step 4: implement. Step 5: run; expect PASS.**
- [ ] **Step 6: Commit** `docs(v9): predeclare data acquisition, readings and lock rule`

### Task 1: Inventory of what is already on disk

**Files:**
- Create: `scripts/v9/inventory.py`
- Output: `results/v9/inventory/{before.csv, README.md}`

- [ ] **Step 1: Write the script:** counts per target from committed files only. Real-fuel freezing: fuels with measured freezing AND class x carbon-number composition (AJFTD GCxGC: 26 have freezing; CRC AV-33-22: 38; plus `data/raw/freezing_lit/`). Mixture systems: systems per property and `n_comp` in `data/unified/measurements_dedup.csv`, split by T < 263 K. Pure constants: `has_antoine/has_lfl/has_viscosity/has_density` for the 7 candidates in `results/v5/BASE/delta_representatives/candidate_gaps.csv`. Cold viscosity: real fuels with viscosity at 253.15 K and 233.15 K.
- [ ] **Step 2: Run** `.venv\Scripts\python.exe scripts/v9/inventory.py`; write `README.md` with the table.
- [ ] **Step 3: Commit** `exp(v9): data inventory before acquisition`

### Task 2: Real-fuel freezing with composition (target 1)

**Files:**
- Create: `data/raw/freezing_v9/{SOURCES.md, <source>_composition.csv, <source>_properties.csv}`, `src/fuelblend/data/loaders_freezing_v9.py`
- Test: `tests/test_loaders_freezing_v9.py`

**Interfaces:**
- Consumes: the class map and `vozka_as_ajftd_composition` shape in `src/fuelblend/data/loaders_freezing_lit.py` (class x carbon-number bins in the AJFTD GCxGC format).
- Produces: `load_freezing_v9() -> tuple[pd.DataFrame, pd.DataFrame]` = (properties in `measurements_dedup` schema, compositions in `ajftd_gcxgc_composition.csv` format) with `source_id` per row.

Candidate sources (search each, log each in `v9_source_log.csv`; this list is where to look, not a claim that each holds usable data):
1. AFRL / NJFCP reference jet fuels (Edwards, "Reference jet fuels for combustion testing", A-1/A-2/A-3, C-1 to C-9): match against AJFTD first.
2. Purdue (Vozka, Kilaz) follow-up papers with GCxGC tables and measured freezing points, beyond Vozka 2018.
3. Heyne group (UDRI, WSU) papers and ASCENT project reports with per-fuel GCxGC and freezing point (ASCENT 065A already ingested).
4. CRC aviation reports with per-sample composition and freezing point beyond AV-33-22 and Report 647.
5. ASTM D4054 qualification research reports (HEFA, ATJ, SIP, CHJ) that are public.
6. Literature on n-alkane distribution vs freezing point of jet fuels (Cookson group and successors), only where per-fuel n-alkane carbon-number data and freezing points are both tabulated.

- [ ] **Step 1: Search and fetch** each candidate; save open PDFs/SI to `extras/staging/lit/v9/` (gitignored), check `%PDF-`; log the status.
- [ ] **Step 2: Transcribe** tables to the raw CSVs (columns per `ingest-common.md`), one source per commit-sized batch; record every judgement in `SOURCES.md`.
- [ ] **Step 3: Failing loader tests:** row counts per source; two spot values per source with page in a comment; every SMILES parses; composition sums to 100 within 0.5 per fuel; `no_carbon_split` flag set where bins are class-only; duplicates against AJFTD/CRC/FuelLib/Vozka resolved to one `system_id`.
- [ ] **Step 4: Run; expect FAIL. Step 5: implement the loader. Step 6: run; expect PASS.**
- [ ] **Step 7: Commit** `data(v9): real-fuel freezing with composition`

### Task 3: Real-fuel cold viscosity (target 4)

**Files:**
- Create: `data/raw/cold_visc_v9/{SOURCES.md, ...}`, `src/fuelblend/data/loaders_cold_visc_v9.py`
- Test: `tests/test_loaders_cold_visc_v9.py`

Sources: the same documents as Task 2 first (most report viscosity at -20 C and -40 C next to freezing); then public reports with low-temperature kinematic viscosity of jet fuels with composition. Fuels without composition are logged `no_composition`.

- [ ] **Steps 1 to 7 as Task 2**, with these test assertions: viscosity stored in cSt after the loader; `T_K` in {253.15, 233.15} or the stated T; a value given in mPa.s is converted only with a measured density at the same T, else `qualifier=dynamic_no_density` and not used.
- [ ] **Commit** `data(v9): real-fuel cold viscosity`

### Task 4: ThermoML archive update for mixture systems (target 2)

**Files:**
- Create: `scripts/v9/fetch_thermoml.py`, `src/fuelblend/data/loaders_thermoml_v9.py`
- Modify: none in `src/fuelblend/data/thermoml_extract.py` unless a test proves a parse gap (then a new test with a real record goes in first)
- Test: `tests/test_loaders_thermoml_v9.py`
- Output: `data/raw/thermoml_v9/` (archive files gitignored), `data/raw/thermoml_v9/SOURCES.md`, `data/raw/thermoml_v9/extracted.csv`

**Interfaces:**
- Consumes: `thermoml_extract.parse_file(root, fname, wanted, max_p, stats)` and `parse_file_sle` unchanged.
- Produces: `load_thermoml_v9() -> pd.DataFrame` in `measurements_dedup` schema, only files not already in `thermoml_2020` / `thermoml_archive_2020` (matched by file name and doi).

- [ ] **Step 1: Fetch** the NIST TRC ThermoML archive journal files published after the 2020 archive (J. Chem. Eng. Data, J. Chem. Thermodyn., Fluid Phase Equilib., Int. J. Thermophys., Thermochim. Acta); record archive date, URL and sha256 in `SOURCES.md`.
- [ ] **Step 2: Failing tests:** one real new record saved as a fixture under `tests/data/` (structure mirrored, not synthetic) parses to the expected T, phase and value; `je8006138.xml` regression stays green; no excess properties; no vapour-phase composition; P <= 1000 kPa.
- [ ] **Step 3: Run; expect FAIL. Step 4: implement. Step 5: run; expect PASS.**
- [ ] **Step 6: Run extraction** for viscosity, density, flash and SLE; write counts by `n_comp` and T < 263 K into `SOURCES.md`.
- [ ] **Step 7: Commit** `data(v9): ThermoML archive update`

### Task 5: Measured pure constants for design candidates (target 3)

**Files:**
- Create: `data/raw/pure_v9/{SOURCES.md, vapour_pressure.csv, lfl.csv, viscosity.csv, density.csv}`, `src/fuelblend/data/loaders_pure_v9.py`, `scripts/v9/refit_antoine_v9.py`
- Test: `tests/test_loaders_pure_v9.py`, `tests/v9/test_refit_antoine_v9.py`
- Output: `results/v9/antoine_refit/refits.csv`

**Interfaces:**
- Consumes: the 7 candidates in `results/v5/BASE/delta_representatives/candidate_gaps.csv` and its `what_would_fill_it` column; Task 4 extraction for vapour pressure; the existing Antoine refit with clipped starts and per-compound `failed:` status.
- Produces: `refits.csv` columns `name, smiles, A_kPa, B, C, T_min_K, T_max_K, n_points, rms_log10, status` (`refit_measured` or `failed:<reason>`); LFL rows with `lfl_source in {measured, tabulated}` (Jones-rule estimates are refused).

- [ ] **Step 1: Collect** vapour pressure from ThermoML (Task 4) and NIST WebBook (Antoine bar to kPa: A_kPa = A_bar + 2); LFL from open tabulations named in `what_would_fill_it`; viscosity and density from ThermoML or Boehm 2022.
- [ ] **Step 2: Failing tests:** a refit with fewer than 5 points or a range missing the 40 to 100 C flash window gets `failed:` (trap 3); falling or flat pressure data does not crash; a WebBook constant converts with +2.
- [ ] **Step 3: Run; expect FAIL. Step 4: implement. Step 5: run; expect PASS.**
- [ ] **Step 6: Run the refit; record which candidates now have measured Antoine + LFL + viscosity + density** (the P4 condition from v8 Task 5).
- [ ] **Step 7: Commit** `data(v9): measured pure constants and Antoine refits for design candidates`

### Task 6: Build `data/unified_v9/`

**Files:**
- Create: `scripts/v9/build_unified_v9.py`
- Test: `tests/v9/test_unified_v9.py`
- Output: `data/unified_v9/{measurements.csv, measurements_dedup.csv, real_fuel_composition_v9.csv, fuel_identity_map.csv, compounds_v9_added.csv, build_log.json}`

**Interfaces:**
- Consumes: `data/unified/measurements.csv` (read only), the real-fuel composition files under `data/real_fuels/v5/` (read only) and the loaders of Tasks 2 to 5.
- Produces:
  - `fuel_identity_map.csv` columns `source_id, source_sample, canonical_system_id, match_key (posf|label|none)`; `canonical_id(source_id: str, sample: str) -> str` in `scripts/v9/build_unified_v9.py`, matching POSF numbers and sample labels against AJFTD, CRC AV-33-22, CRC 647, FuelLib and Vozka 2018.
  - `measurements_dedup.csv` with the dedup columns of AGENTS.md (`vmin, vmax, q25, q75, n_obs, n_refs, n_sources, sources, n_outlier, conflict`) plus `v9_new` (True for groups with no v5 row). New-source rows of a known fuel use its canonical id.
  - `real_fuel_composition_v9.csv`: v5 GCxGC compositions plus Task 2/3 compositions, in the `ajftd_gcxgc_composition.csv` format, keyed by canonical id, with `no_carbon_split`.

- [ ] **Step 1: Failing tests:** `canonical_id("afrl", "POSF-10325") == ` the AJFTD fuel 2 id; a new freezing value for a known AJFTD fuel lands on that fuel's id and is not `v9_new` at fuel level; every v5 dedup row appears unchanged when no new row shares its key; `v9_new` groups have no v5 row; conflicts flagged with the AGENTS.md tolerances; `1,2,4-trimethylbenzene` stays quoted in the compounds file; `data/unified/` mtime and hash unchanged after the build.
- [ ] **Step 2: Run; expect FAIL. Step 3: implement. Step 4: run; expect PASS.**
- [ ] **Step 5: Run the build; then** `.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py`; expect pass.
- [ ] **Step 6: Commit** `data(v9): unified_v9 layer with dedup`

### Task 7: Audit checks and after-inventory

**Files:**
- Create: `scripts/v9/audit_v9.py`
- Output: `results/v9/inventory/{after.csv, delta.csv}`, `results/v9/audit/report.md`

- [ ] **Step 1: Audit per new source:** units in canonical ranges (sanity bounds for catching unit errors, not results: density 600 to 1100 kg/m3, freezing -120 to 20 C, viscosity at -40 C 1 to 50 cSt, flash 0 to 150 C), no `data/quarantine/` or `DatasetsFuelProperties.xlsx` rows, provenance columns filled, every `SOURCES.md` has a licence line, no PDFs staged in git.
- [ ] **Step 2: Re-run `inventory.py` on `data/unified_v9/`**; write `delta.csv` (before, after, added per target).
- [ ] **Step 3: Commit** `exp(v9): audit and inventory after acquisition`

### Task 8: Lock the new data, dev readings, then the lock (once)

**Files:**
- Create: `scripts/v9/make_lock.py`, `scripts/v9/readings.py`
- Modify: `src/fuelblend/data/model_table.py`, `src/fuelblend/splits/locked_sets.py` (read `FUELBLEND_DEDUP_PATH`, default unchanged)
- Output: `results/v9/lock.json`, `results/v9/readings/{R9-1..R9-4}.csv`, `results/v9/readings/READING.md`

**Interfaces:**
- Consumes: `is_v9_locked` (Task 0), `fuelblend.v7.guard.assert_dev_only` with every existing lock named, the v5 designated routes, the v5 `ml_pair_residual` training code, the WS6 F1x real-fuel freezing reading.

**Data path.** `src/fuelblend/data/model_table.py` and `src/fuelblend/splits/locked_sets.py` hard-code `data/unified/measurements_dedup.csv`. Add an environment variable `FUELBLEND_DEDUP_PATH` read by both (default unchanged, so BASE reproduces); `readings.py` sets it to `data/unified_v9/measurements_dedup.csv` and reads real-fuel composition from `data/unified_v9/real_fuel_composition_v9.csv`.

Readings (pass rules fixed in Task 0, same statistics as v5: paired group bootstrap over systems or fuels, 90 % CI):
* R9-1 real-fuel freezing (dev, old plus new fuels with carbon split): UNIFAC SvL and F1x vs leave-one-group-out constant; skill shown only if the CI of the gain excludes 0.
* R9-2 viscosity pair residual retrained on v5 plus new dev systems vs the v5 model, on the v5 dev systems (does more data help).
* R9-3 real-fuel cold viscosity at 253.15 K and 233.15 K: Refutas and the WS6 constant-offset model vs constant.
* R9-4 density and flash routes on new dev mixture systems (transfer of v5 routes to new data; descriptive).

- [ ] **Step 1: Run** `make_lock.py` once: lock = `v9_new` groups with `is_v9_locked`; write `lock.json`; refuse if it exists.
- [ ] **Step 2: Tests** `tests/v9/test_readings_guard.py`: `readings.py --set dev` raises `LockedRowError` if any v9 or earlier locked id is passed; with `FUELBLEND_DEDUP_PATH` unset, `model_table` reads `data/unified/measurements_dedup.csv` (BASE path); the R9-2 training set's system count equals the v9 dev count, not the v5 one; then run `reproduce_base.py`, expect pass.
- [ ] **Step 3: Run dev readings; write `READING.md` with each verdict using the Task 0 rule verbatim; commit** `exp(v9): dev readings on expanded data`.
- [ ] **Step 4: Run `readings.py --set v9_lock` once; append to `READING.md`; commit** `exp(v9): v9 lock readings (scored once)`.

### Task 9: Numbers, docs, full suite, push

**Files:**
- Create: `scripts/v9/build_final_numbers_v9.py`, `docs/final_numbers_v9.md` (generated)
- Modify: `docs/methodology_v7.md` (section "v9 data"), `docs/GUIDE.md` (status, data sources), `docs/v5_negative_results.md` (failed readings)

- [ ] **Step 1: Generator reads only committed v9 CSVs:** inventory delta, source log summary, R9-1 to R9-4 dev and lock tables; status lines "status: no skill shown" where a reading fails.
- [ ] **Step 2: Run it; commit** `docs(v9): final numbers and methodology`.
- [ ] **Step 3: Full suite** `.venv\Scripts\python.exe -m pytest -q > scratch/v9_full.log`; expect 741 + new passed, 11 skipped, only the known `test_committed_splits_match_rebuild` failure.
- [ ] **Step 4:** `git status --short`; confirm no PDF, page image, `ajftd_session.json` or licensed file is staged.
- [ ] **Step 5:** `git push origin feat/external-datasets-ingest` (skipped if `results/v9/STOPPED.md` exists).
- [ ] **Step 6: Report:** fuels and systems added per target, sources ingested vs not found, R9-1 to R9-4 dev and lock verdicts, and whether real-fuel freezing now shows skill (the certification bottleneck).
