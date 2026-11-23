# v10 PDF Transcription Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans (native, no subagents). Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** turn source PDFs that the user downloads into verified rows in a new data layer `data/unified_v10/`, then score the v10-new data once on a fresh lock. Targets are real-fuel freezing with composition, real-fuel cold viscosity, and measured pure constants for the 7 design candidates.

**Architecture:** v9 failed for one reason: the drafts could not be checked against a primary file. This happened for Edwards 2017 and for the pure constants, because DTIC and ROSAP served HTML in place of PDFs. In v10 the user puts the PDFs on disk first. Gemini then drafts a page index and transcriptions, quoting each value's page and table. Sonnet checks every value that can move a reading against the page itself (text layer or rendered image), then builds `unified_v10` on top of `unified_v9`, locks the v10-new groups and scores them once.

**Tech Stack:** Python 3.12 (`.venv`), pandas, pdfplumber (text layer, `page.to_image`), pypdf, pytest. All are already installed.

**Spec:** user decisions below; the v9 plan `docs/superpowers/plans/2026-11-23-v9-data-acquisition.md`; the verification record `results/v9/verification.md`.

## User decisions (recorded 2026-11-23, before any v10 work)

* **PDFs:** the user downloads these into `data/raw/v10_pdfs/` before phase S0 starts:
  * AFRL-TR-2020-0017 (Edwards, Jet Fuel Properties, DTIC AD1093317).
  * Edwards AIAA 2017-0146 (Reference Jet Fuels for Combustion Testing).
  * The primary papers for the 7 design candidates' constants.
* **Who does what:** Gemini 3.8 Flash drafts, Sonnet checks every row that matters against the page.
* **Scoring:** new lock on v10-new groups only, `md5(group_id) mod 4 == 0`, scored once, with an 8-fuel floor.
* **Carried over from v9:**
  * Commit after each task and push once at the end to `origin feat/external-datasets-ingest`.
  * Run natively, unattended, with no subagents.
  * Raw PDFs are never committed.
  * No web fetching in v10. The only inputs are the PDFs on disk.

## Phase U: the user's step (before any session)

1. Create `data/raw/v10_pdfs/` and save the PDFs there. Use the browser's own "Save as" so the file is the real PDF, not a viewer page.
2. Add `data/raw/v10_pdfs/sources.txt`, one line per file: `filename | short title | where it came from (DTIC, AIAA, ACS ...)`.
3. For the pure-constant papers, a useful list to look for is the citations in `data/raw/pure_v9/*_draft.csv` (column `doi` or `ref`). Any paper with measured vapour pressure, viscosity, density or flammability limit for these candidates counts:
   * 1-methylnaphthalene
   * 2,6-dimethylnaphthalene
   * perhydrophenanthrene
   * perhydrofluorene
   * 1-methyltetralin
   * n-butyltetralin
   * 2-methyldecalin

The plan never asks the user anything after this. A missing PDF is logged as `not provided` and the run continues.

## Model assignment

| Phase | Model | Tasks | May write | May never |
|---|---|---|---|---|
| S0 | Sonnet | Task 0 | `docs/v10_predeclaration.md`, `src/fuelblend/v10/lock.py`, `results/v10/pdf_register.csv`, tests | read any PDF content beyond headers and page counts |
| G | Gemini 3.8 Flash | Tasks 1 to 4 | `data/raw/*_v10/*_draft.csv`, loaders under `src/fuelblend/data/loaders_*_v10.py`, their tests and fixtures, `results/v10/page_index.csv`, `results/v10/HANDOFF_gemini.md`; local commits prefixed `data(v10-draft):` or `feat(v10-draft):` | touch `results/v10/lock.json`, readings, `data/unified_v10/`, `data/unified*/`, BASE, `splits/`, `configs/`, `model_table.py`; invent a `system_id`; push |
| S1 | Sonnet | Tasks 5 to 9 | everything else in this plan | alter a number to fit a check |

Gemini ends phase G by writing `results/v10/HANDOFF_gemini.md`. It covers:
* per PDF, which tables it transcribed and their row counts;
* the count of `unsure` rows;
* every place where the text layer and the page image disagreed;
* what it left undone.

## Unattended rules

* **Failed checks:** a failed check or gate is recorded and the run continues.
* **Stopping:** stop only on a broken invariant:
  * a BASE hash mismatch,
  * a lock file that already exists,
  * a content change in `splits/locked.json`,
  * a raw PDF staged for commit.

  Stopping writes `results/v10/STOPPED.md` and does not push.
* **`splits/locked.json`:** a CRLF-only diff is restored with `git checkout -- splits/locked.json` and noted in the handoff.
* **Staging:** use `git add <explicit paths>` only, never `git add -A` or `git add .`. v9 swept unrelated `.docx` abstracts into a commit this way.

## Global Constraints

* Canonical units: cSt, mPa.s, kg/m3, degC, kPa. T is in K inside equations. Conversions happen in loaders only.
* Antoine form: log10(P/kPa) = A - B/(T_K + C). NIST bar to kPa: A_kPa = A_bar + 2.
* Carbon split: R10-1 and R10-2 need composition by class x carbon number (the AJFTD GCxGC layout that `vozka_as_ajftd_composition()` emits). A class-only table (D2425 hydrocarbon types) is kept, flagged `class_only`, and does not count toward a floor.
* Real-fuel ids come only from `canonical_id(source_id, sample)` in `scripts/v9/build_unified_v9.py`, which matches POSF numbers against the AJFTD label table. A fuel with no match gets `fuel:v10:<source_id>:<sample>`, nothing else.
* Conflict tolerances against v5 or v9 rows for the same fuel:
  * Freezing: 1 degC.
  * Flash: 2 degC.
  * Otherwise: 5 % relative.

  Conflicts are flagged, never deleted.
* Flammability limits: measured or tabulated values only, with page and table. No Jones-rule or homologue estimates.
* Never evaluate un-refit Antoine constants for C10+ alkanes in the 40 to 100 degC window.
* Every new source needs a loader, a unit test with a fixture that mirrors the real file layout, and an audit check.
* Never alter numbers. Every figure in docs traces to a committed results file.
* Do not checkout or normalise `results/v5/BASE/**`.
* No AI attribution in commits.

## Review Focus

1. **A DTIC consent page saved as `.pdf`:** Task 0 must reject a file that does not start with `%PDF-` (test `test_register_rejects_html`).
2. **A scanned PDF with no text layer:** Gemini must mark every row from it `unsure` and render the page. Sonnet must then read every row from the image (test `test_scanned_rows_forced_unsure` in Task 4).
3. **Negative signs lost in extraction:** freezing points like `-47` read as `47` or with an en dash. The loader must accept `−` and `–`, and the audit must reject jet freezing points above 0 degC (test `test_minus_sign_variants`).
4. **A fuel that is in AJFTD under a POSF label:** it must map to the existing `fuel:<id>` and stay dev, not enter the v10 lock as new (test `test_posf_maps_to_ajftd`).
5. **Composition tables that do not close, or have no carbon split:** a fuel whose class sums are outside 100 ± 0.5 mass % is kept but flagged `composition_open`. A class-only table is flagged `class_only`. Both are excluded from R10-1 and R10-2 (tests `test_open_composition_flagged`, `test_class_only_flagged`).

---

### Task 0: Precondition, PDF register, predeclaration, lock rule (S0, Sonnet)

**Files:**
- Create: `docs/v10_predeclaration.md`
- Create: `src/fuelblend/v10/__init__.py`, `src/fuelblend/v10/lock.py`
- Create: `scripts/v10/register_pdfs.py`, which writes `results/v10/pdf_register.csv`
- Test: `tests/v10/test_lock.py`, `tests/v10/test_register.py`

**Interfaces:**
- Produces:
  - `is_v10_locked(group_id: str) -> bool`: `int(md5(group_id).hexdigest(), 16) % 4 == 0`.
  - `v10_lock_ids(new_ids: Iterable[str], existing_ids: Iterable[str]) -> frozenset[str]`: the locked subset of `new_ids` minus `existing_ids`.
  - `register(pdf_dir: Path) -> pd.DataFrame` with columns `filename, sha256, is_pdf, n_pages, has_text_layer, title_from_sources_txt, status`. `status` is `ok`, `not_pdf` or `not provided`.

- [ ] **Step 0:** confirm v9 is closed. `git fetch origin` then `git merge-base --is-ancestor e46bb3c origin/feat/external-datasets-ingest` must exit 0. Otherwise write `STOPPED.md`.
- [ ] **Step 1:** write the failing tests.
  * `test_register_rejects_html`: a fixture file beginning `<!DOCTYPE` gives `is_pdf False, status not_pdf`.
  * `test_register_hash`: sha256 matches `hashlib` on the fixture.
  * `test_lock_rule_matches_v9_form`: `is_v10_locked` equals the v9 formula on 20 ids.
  * `test_lock_only_new`: `v10_lock_ids(new_ids, existing_ids)` never returns an id in `existing_ids` (in the run, the `system_id` set of `data/unified_v9/measurements_dedup.csv`).
- [ ] **Step 2:** implement both modules. `has_text_layer` is True when pdfplumber extracts at least 200 characters from the first 3 pages.
- [ ] **Step 3:** run `register_pdfs.py`. Expected outcomes:
  * Each file the user named is `ok`.
  * Expected but absent files are listed as `not provided`. These are AFRL-TR-2020-0017, AIAA 2017-0146, and each pure paper named in `sources.txt`.
  * A `not_pdf` file is logged and skipped. It is not a stop.
- [ ] **Step 4:** write `docs/v10_predeclaration.md`. It copies the user decisions, the lock rule (v10-new means a group with no row in `unified_v9`), and this readings table, fixed before any value is read:

| reading | pass rule |
|---|---|
| R10-1 real-fuel freezing | model: the v5 designated UNIFAC Schroeder-van Laar route (`unifac_desig`, as in `scripts/v8/g2_freezing_check.py`). Constant: the leave-one-group-out mean freezing point of the dev fuels (column `C` there). Skill shown only if the 90 % paired group bootstrap CI (5000 resamples, seed 42) of (constant MAE minus model MAE) excludes 0 on the v10 lock, with at least 8 v10-new fuels that have a carbon-split composition that closes. Below the floor the reading is descriptive and the status stays "no skill shown". A dev-only gain cannot reopen it |
| R10-2 real-fuel cold viscosity at 253.15 K and 233.15 K | model: the v5 Refutas route on the carbon-split composition. Constant: the leave-one-group-out mean ln(cSt) of dev fuels at that temperature. Same bootstrap, rule and floor as R10-1, per temperature |
| R10-3 agreement with existing rows | descriptive: for fuels already in v5 or v9, the count and size of differences beyond the conflict tolerance |
| R10-4 design candidates with measured constants | descriptive, two counts: candidates with measured LFL, viscosity, density and a converged Antoine refit; and the subset whose refit T range overlaps the flash window 313 to 373 K (only these meet the P4 condition). No design is re-run in v10 |

- [ ] **Step 5:** run `pytest tests/v10 -q` and confirm it passes. Then commit with `git add` on the explicit paths and the message `chore(v10): PDF register, predeclaration, lock rule`. Write `results/v10/HANDOFF_S0.md` with the register summary.

### Task 1: Page index (G, Gemini)

**Files:** create `results/v10/page_index.csv` with columns `filename, page, table_label, content, target, extraction`, where:
* `target` is one of `freezing`, `cold_visc`, `composition`, `pure_vp`, `pure_visc`, `pure_density`, `pure_lfl`, `other`;
* `extraction` is `text_layer` or `image_only`.

- [ ] **Step 1:** for every `ok` PDF, list each table or figure-with-numbers that carries a target property. Figures with no tabulated numbers are `other`. Values are never read off plots.
- [ ] **Step 2:** render every indexed page to `scratch/v10_pages/<filename>_p<page>.png` at 200 dpi with `pdfplumber` `page.to_image(resolution=200)`. `scratch/` is never committed.
- [ ] **Step 3:** commit `page_index.csv` with `data(v10-draft): page index`.

### Task 2: AFRL-TR-2020-0017 and AIAA 2017-0146 real-fuel drafts (G, Gemini)

**Files:**
- Create: `data/raw/fuels_v10/<source_id>_properties_draft.csv`, `data/raw/fuels_v10/<source_id>_composition_draft.csv`, and `data/raw/fuels_v10/SOURCES.md`. The source ids are `afrl2020` and `edwards2017`.
- Create: `src/fuelblend/data/loaders_fuels_v10.py`
- Test: `tests/test_loaders_fuels_v10.py`, with a fixture `tests/data/fuels_v10_fixture.csv` that copies 5 real rows including one en-dash negative.

**Interfaces:**
- Draft columns are those of the v9 drafts: `source_id, table_or_page, sample_id, kind, property, value, unit, T_C, method, qualifier, note, page, table, confidence`. Add `raw_text` (the cell exactly as printed) and `extraction` (`text_layer` or `image`).
- Produces:
  - `load_fuel_properties_v10(path) -> pd.DataFrame` in the `measurements.csv` schema.
  - `load_fuel_composition_v10(path) -> pd.DataFrame` in the AJFTD GCxGC long layout used by `vozka_as_ajftd_composition()`.

- [ ] **Step 1:** transcribe every freezing, cold-viscosity, density and flash value for every fuel in the tables indexed in Task 1, and every composition table (class x carbon number, or hydrocarbon type).
  * `raw_text` holds the printed cell.
  * `confidence` is `unsure` when any of these holds:
    * the page is `image_only`;
    * text and image disagree;
    * the sign or decimal is unclear;
    * a cell spans rows.
- [ ] **Step 2:** write the tests:
  * `test_minus_sign_variants`: `−47`, `–47` and `-47` all load as -47.0.
  * `test_open_composition_flagged`: a fixture fuel summing to 97 % gets `composition_open True`.
  * `test_class_only_flagged`: a hydrocarbon-type table with no carbon numbers gets `class_only True`.
  * `test_posf_maps_to_ajftd`: `POSF-10325` maps through `canonical_id` to the AJFTD id. The test reads the id from the label table; it is not hard-coded.
  * `test_units`: cSt kept, mm2/s maps to cSt, degF converted.
- [ ] **Step 3:** implement the loaders, run the tests, and commit with `data(v10-draft): AFRL 2020 and Edwards 2017 drafts`.

### Task 3: Pure-constant drafts for the 7 candidates (G, Gemini)

**Files:**
- Create: `data/raw/pure_v10/{vapour_pressure,viscosity,density,lfl}_draft.csv` and `data/raw/pure_v10/SOURCES.md`
- Create: `src/fuelblend/data/loaders_pure_v10.py`
- Test: `tests/test_loaders_pure_v10.py`

- [ ] **Step 1:** transcribe only from the downloaded papers. Use the same columns as Task 2 plus `compound, inchikey`. Take the InChIKey from `data/unified/compounds.csv` by name, or leave it blank.
  * Pressures in mmHg or bar are converted in the loader, not in the CSV.
  * A value cited by a paper from another work is kept only with `qualifier=secondary` and the original reference in `note`.
- [ ] **Step 2:** write the tests:
  * `test_vp_units`: mmHg, bar and Pa all map to kPa.
  * `test_lfl_measured_only`: rows whose `method` contains `estimate`, `Jones` or `homologue` are rejected.
  * `test_scanned_rows_forced_unsure`: rows with `extraction=image` and `confidence=clear` load with `confidence` forced to `unsure`.
- [ ] **Step 3:** implement the loaders, run the tests, and commit with `data(v10-draft): pure constants drafts`.

### Task 4: Gemini handoff (G, Gemini)

- [ ] **Step 1:** write `results/v10/HANDOFF_gemini.md` with the content listed under Model assignment.
  * Report row counts from the CSVs by a script, not from memory. v9's summary disagreed with its own CSV.
  * By the same script, report how many fuels are v10-new (no POSF match in the AJFTD label table) and have both a freezing point and a closed carbon-split composition. R10-1 needs about 32 such fuels for 8 to fall in a 25 % lock; AJFTD is AFRL's own data, so many AFRL-TR-2020-0017 fuels will map to existing ids and stay dev.
  * End phase G. Do not push.

### Task 5: Verification against the page (S1, Sonnet)

**Files:** create `results/v10/verification.md` and `results/v10/verification_rows.csv` with columns `file, row, value_draft, value_page, match, checked_by`, where `checked_by` is `text_layer` or `image`.

- [ ] **Step 1:** check every row with property freezing, cold viscosity (T at or below 253.15 K) or any pure constant against the rendered PNG of its page. The text layer may be used as a second read, never as the only one, because Gemini may have drafted from the same text layer: it misses column shifts (a value on the wrong fuel) and dropped minus signs. For each table, also confirm that the column headers and the row order match the image.
  * Record each check in `verification_rows.csv`.
- [ ] **Step 2:** check composition and other property tables.
  * Check every `unsure` row plus a random 20 % of `clear` rows (`random_state=42`).
  * If any sampled row is wrong, re-check the whole table.
  * Compare composition sums against 100 ± 0.5.
- [ ] **Step 3:** cross-check sample ids. Check each POSF and sample id against the report's own fuel list and against Boehm 2022, where the same POSF appears. v9 found a draft C-1 labelled POSF-11498 against Boehm's POSF 13718.
- [ ] **Step 4:** finalise the drafts.
  * Correct each wrong value to the page value and log it.
  * Rename the verified drafts without the `_draft` suffix.
  * Unverifiable rows stay in `_draft` files and are logged `unverified`.
  * Review Gemini's loaders for invented ids and silent unit changes.
- [ ] **Step 5:** commit with `data(v10): verified transcriptions`.

### Task 6: Antoine refits for candidates with verified vapour pressure (S1)

- [ ] **Step 1:** run `scripts/v9/refit_antoine_v9.py` on the verified `pure_v10/vapour_pressure.csv`, with output to `results/v10/antoine_refit/refits.csv`. If the script hard-codes its input path, add a `--input` and `--out` argument and a test that the default is unchanged.
  * Each candidate gets `ok` or `failed:<reason>`.
  * A fit needs at least 4 points spanning at least 40 K. Record `T_min_K` and `T_max_K` of each fit and a column `in_flash_window` (True when the range overlaps 313 to 373 K). Trap 3: a fit outside that window is not used for flash.
- [ ] **Step 2:** commit with `exp(v10): Antoine refits`.

### Task 7: Build `data/unified_v10/` (S1)

**Files:**
- Create: `scripts/v10/build_unified_v10.py`, which reuses `canonical_id` and the dedup rules from `scripts/v9/build_unified_v9.py` by import, not by copying.
- Test: `tests/v10/test_unified_v10.py`

- [ ] **Step 1:** write the tests:
  * `test_v9_rows_unchanged`: every v9 dedup row is present with the same value.
  * `test_v10_new_flag`: `v10_new` is True only for groups with no v9 row.
  * `test_conflicts_flagged_not_dropped`: an injected 3 degC freezing disagreement on an AJFTD fuel stays, with `conflict True`.
  * `test_no_unverified_rows`: no row comes from a `_draft` file.
- [ ] **Step 2:** implement and build. Write `data/unified_v10/measurements_dedup.csv`, `fuel_identity_map_v10.csv` and `real_fuel_composition_v10.csv`. Force-add only those three, plus the build log, since `data/*` is gitignored. Run `python scripts/v5/base/reproduce_base.py` and confirm it passes.
- [ ] **Step 3:** commit with `data(v10): unified_v10 layer`.

### Task 8: Audit, lock, readings (S1)

- [ ] **Step 1:** run the audit. Copy `scripts/v9/audit_v9.py` to `scripts/v10/audit_v10.py` and add these checks:
  * no file under `data/raw/v10_pdfs/` is tracked;
  * jet freezing values lie between -80 and 0 degC;
  * cold viscosity values lie between 1 and 30 cSt;
  * every row has a page and table.

  Write the report to `results/v10/audit/report.md`. Run `scripts/v9/inventory.py` against `unified_v10`, before versus after.
- [ ] **Step 2:** write the lock. Write `results/v10/lock.json` once, using `v10_lock_ids(new_ids, v9_ids)` over v10-new groups. Refuse to run if the file exists.
- [ ] **Step 3a:** before any real reading, write `tests/v10/test_readings.py::test_r10_1_runs_on_fixture`: a fixture of 10 synthetic fuels with carbon-split compositions and freezing points runs the R10-1 path end to end (route prediction, constant, bootstrap) and returns a finite CI and `n=10`. In v9 every R9 reading had n=0, so this branch has never run. Fixture fuels are never written to any data file.
- [ ] **Step 3:** run the readings. Create `scripts/v10/readings.py`, adapted from `scripts/v9/readings.py`, with `--set dev|v10_lock` and `FUELBLEND_DEDUP_PATH=data/unified_v10/measurements_dedup.csv`.
  * Exclude fuels locked in v5 to v9 from dev, as v9 did.
  * Run the dev set first, then the lock once.
  * Apply the predeclared rules unchanged.
  * Write `results/v10/readings/READING.md`.
- [ ] **Step 4:** commit with `exp(v10): lock and readings (scored once)`.

### Task 9: Numbers, docs, full suite, push (S1)

- [ ] **Step 1:** write the numbers and docs.
  * Create `scripts/v10/build_final_numbers_v10.py`, which writes `docs/final_numbers_v10.md` from committed files only.
  * Add a v10 section to the methodology doc and to `docs/v5_negative_results.md` if a reading is "no skill shown".
  * If R10-4 is at least 1, add one line saying the design re-run (v8 Tasks 5 and 6) can now run for those candidates.
- [ ] **Step 2:** run the full suite with `.venv\Scripts\python.exe -m pytest -q`. Expected: everything passes except the known `test_committed_splits_match_rebuild`.
- [ ] **Step 3:** clean up and push.
  * Restore any CRLF-only change to `splits/locked.json`.
  * Run `git status` and confirm no PDF, page image, `.docx` or session file is staged.
  * Commit with `docs(v10): final numbers and methodology`.
  * Push to `origin feat/external-datasets-ingest`.
