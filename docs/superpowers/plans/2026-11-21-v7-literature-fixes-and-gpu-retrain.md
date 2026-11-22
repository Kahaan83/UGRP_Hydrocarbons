# v7: literature fixes and instrumented GPU retrain, implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (Native; AGENTS.md forbids subagents unless asked) to implement this plan task by task. Steps use checkbox (`- [ ]`) syntax for tracking. Each phase is runnable in its own session; start each session by reading this file and the Interfaces blocks of earlier tasks.

**Goal:** Turn the two v6 literature reviews into gated experiments, and retrain every trainable model with per-epoch learning curves, sample-size learning curves, longer training and nested hyperparameter search, the torch models on the RTX 3050 GPU.

**Architecture:** A shared curve logger (CSV, plus optional TensorBoard) is hooked into the four torch training loops and into LightGBM, with defaults unchanged. A lock guard stops v7 code from scoring any row already scored in v6. Cheap diagnostics (epoch curves of the v1 config, a descriptor GBM comparator, a Chemprop parity check, sample-size curves) run first and set the HPO budget. Then come GPU and CPU retrains with Optuna inside inner folds, and the literature methods (HANNA gamma for flash, freezing by family, equation-parameter heads, family-conditional conformal). All outputs go to `results/v7/`, and nothing is adopted into BASE or into `predict_fuel()` without a declared gate and the user's yes.

**Tech Stack:** Python 3.12, torch (CUDA 12.8 wheel in a separate GPU venv), RDKit, LightGBM, scikit-learn, Optuna, TensorBoard, Chemprop v2 (parity only).

**Spec:** `docs/v6_failure_analysis_and_ml_options.md` (F1 to F4, sections 2 to 4b) and `docs/literature_gaps_v6.md` (gaps 1 to 8, fixes 3.1 to 3.6). Numbers: `docs/v6_summary.md`, `docs/final_numbers_v5.md`.

## Global Constraints

- Never invent or alter numbers. Every figure traces to a committed results file. State no expected gains; write "status: no skill shown" where a model does not beat its comparator.
- BASE (`results/v5/BASE/`, `SHA256SUMS_v5_base.txt`) is read-only. Never `git checkout` it, normalise its line endings or regenerate it.
- Never overwrite the production model files that `predict_fuel()` loads (`src/fuelblend/ml/*` `MODEL_DIR/model.txt`, the `fusion_estimate` model).
- Every v6 holdout is spent: the YSI holdout, the freezing holdout (12 fuels), the mixture `lock.json` heads, the 51-row replication and the flash oxygenate `lock.json`. v7 numbers are dev-only and exploratory unless they come from a new lock declared and committed before fitting.
- `GroupKFold` on `system_id` for mixtures and fuels, and on generic Murcko scaffold for pure compounds (as in v6 Track B). Hyperparameters are chosen in inner folds only; each outer fold is scored once with its chosen config.
- CPU and GPU numbers are never mixed in one comparison (CUDA `index_add_` is nondeterministic). Each results file records `device`, `torch_version` and `cuda_version`.
- Flash point: report one-sided over-prediction alongside MAE (Trap 6).
- Multiplicity: every lock declares its family of comparisons and applies Holm at alpha 0.10 (one-sided) across that family.
- The GPU venv, the TensorBoard run directories and the Optuna storage live outside OneDrive under `C:\fuelblend-gpu\` (env var `FUELBLEND_RUNS`, default `C:\fuelblend-gpu\runs`). Commit only CSV curve summaries and PNGs.
- Tests: `.venv\Scripts\python.exe -m pytest -q <file>` per task. Run the full suite at phase ends; the baseline is 741 passed, 11 skipped.
- Academic, direct tone, no em dashes in docs. No AI attribution in commits.

## Review Focus

1. **A v7 script scores a v6-locked row.** Expect a `LockedRowError` naming the lock file. Pinned in Task 1 (`test_guard_refuses_each_lock_file`).
2. **The instrumented loop changes v6 results when logging is off.** With `logger=None` on CPU at a fixed seed, predictions must be bitwise equal to the current code. Pinned in Task 2 (`test_fit_ensemble_logger_none_bitwise`, `test_molpool_defaults_bitwise`).
3. **CUDA is missing or out of memory on the 4 GB card.** The script must fail with a clear message, or halve the batch size once and log it, and never fall back silently to CPU while the files say `cuda`. Pinned in Task 3 (`test_device_resolution_reports_actual`) and Task 10 (`test_oom_retry_halves_bs_once`).
4. **An Optuna study resumed after a crash or a laptop sleep.** Expect it to continue from the SQLite storage without repeating finished trials, and the outer score to be computed once per fold. Pinned in Task 9 (`test_study_resume_no_duplicate_trials`).
5. **The curve CSV is written while a monitor reads it.** Rows must be appended and flushed per epoch, so that a partial file always parses. Pinned in Task 2 (`test_curve_csv_partial_read`).

---

## Merged findings to task map

| Finding (source) | Task |
|---|---|
| Track B compared only to a constant; add RDKit-descriptor RF/GBM (failure F3; gaps doc gap 1). Deferred earlier for ordering, not dropped | 6 |
| Hand-written D-MPNN not checked against Chemprop (gaps doc gap 2) | 7 |
| Data scale is the cross-cutting limit (failure doc section 2.2; gaps doc table) | 8 |
| More epochs and optimisation (user request) | 5, 9 to 13 |
| HANNA gamma vs UNIFAC, hydrocarbon binaries, flash route only, pretrained, no refit (failure doc sections 4.1 and 4b) | 17 |
| Freezing: limiting-species baseline (d) has no number (gaps doc gap 4) | 14 |
| Freezing: leave-one-family-out, heavy-tail GP, novelty flag, interval-only output (F1 fixes 1, 3, 4) | 15 |
| Freezing: joint pure plus fuel MolPool, M-QSPR style (F1 fix 2; gaps doc gap 5) | 16 |
| Equation-parameter heads (Antoine, Andrade) instead of point values (F3; ExPUFFIN) | 18 |
| Family-conditional and shift-weighted conformal, coverage by family (failure doc section 2.4; gaps doc 3.4) | 19 |
| Data scale is the limit, so which data, and where would it come from (follow-up, 2026-11-21) | 0.5 |
| Fresh lock needed for any claim (all holdouts spent) | 20 |
| Multiplicity rule, citation metadata, link docs (gaps doc 3.3, 3.5, 3.6). Fixes 3.1 and 3.2 are already in the uncommitted `docs/v6_summary.md` | 21 |
| Deferred: COSMO-SAC mixture pretraining (F2 fix 2), CheMeleon encoder (gap 3), tabular foundation model (gap 7), GE-GNN flash oxygenate comparator (F4, gap 6) | After Task 17, see "Deferred" |
| Not recommended: further MolPool tuning on the same 23 to 241 rows; any head scored on a set already seen | Enforced by Task 1 |

---

## Phase 0: guardrails

State at plan time (2026-11-21): `SHA256SUMS_v5_base.txt` verifies 23 of 23 with CRs stripped. The 24 modified `results/v5/BASE/delta_*` files and `splits/locked.json` differ from HEAD only in line endings (`git diff --ignore-cr-at-eol` is empty), and the `delta_*` files are not in the manifest. Leave them as they are; do not commit them.

Step 0 (ask the user first): both spec docs are untracked, and `docs/v6_summary.md` carries uncommitted edits (fixes 3.1 and 3.2). Every phase session reads them, so ask whether to commit them, with this plan, before Task 1.

### Task 0.5: Data inventory and data-source research

Purpose: Task 8 can show that more data helps, but not which data. This task answers that question in advance. It has three parts: what we already hold but don't use, where the cited papers got their data, and what else exists online. It produces no model and no number about our models. Ordering: Part A uses `guard` (Task 1) for `n_locked`, so run Part A right after Task 1. Parts B and C have no code dependency and can run first.

**Authorisation (user, 2026-11-21):** web search, web fetching, Consensus searches and the deep-research skill (`academic-research-skills:deep-research`, or `anthropic-skills:deep-research`, which may spawn research subagents) are allowed for this task only. AGENTS.md rule 6 still applies: search queries name properties, fuel types and methods, never our project, names or institution. Downloading anything (datasets, PDFs over a few MB, supplementary files) needs a separate yes from the user. Licensed data is listed, never downloaded or committed.

**Files:**
- Create: `scripts/v7/data/inventory_local.py`, `results/v7/data_inventory/local_counts.csv`, `results/v7/data_inventory/unused_candidates.csv`
- Create: `docs/v7_data_sources.md` (the paper-by-paper and online survey) and `docs/v7_data_plan.md` (a ranked acquisition list)
- Test: `tests/test_v7_inventory.py`

**Part A: local inventory (no network).**
- `inventory_local.py` writes one row per `property, kind` to `local_counts.csv`, with these columns:
  - `n_rows`, `n_systems`, `n_fuels`, `n_compounds`, `n_scaffolds`
  - `n_in_v6_tables`: rows the v6 tables actually use
  - `n_unused`: rows that pass the `is_hc` filter and dedup but no v6 table uses
  - `n_locked`: rows inside a lock, from `guard`
- Questions it must answer, each as a row in `unused_candidates.csv` with an id, the data available and the reason it is unused:
  1. AJFTD fuels with GCxGC class by carbon-number composition (the detail file has 113 fuels) and a measured freezing point that are not among the 64 v6 freezing groups.
  2. AJFTD or other real fuels with composition and a viscosity measured at or below −20 °C.
  3. ThermoML and NIST TDE hydrocarbon mixture systems (`n_comp` ≥ 2) for viscosity, density and flash that are not in `results/v6/mixture/rows.csv` or `results/v6/flash_oxy/rows.csv`. List them by `system_id`, with component classes (n-alkane, iso-alkane, cycloalkane, aromatic) and the temperature range.
  4. Pure compounds in `data/unified/` with any of the 10 Track B targets that are not in `results/v6/gnn/pool.csv`.
  5. Pure compounds with vapour pressure or viscosity at three or more temperatures (inputs for the Task 18 parameter heads).
  6. Oxygenates in the flash rows that have no Antoine or LFL constants.
  7. What `data/interim/external_staging/` adds (rows, systems, overlap with `measurements_dedup.csv`).
- Rules: read with column subsets and short pandas summaries, never bulk-read CSVs (AGENTS.md usage limits), and never read `data/quarantine/` as a source.
- [ ] **Step 1: Write the failing tests.** `test_counts_sum` (`n_in_v6_tables + n_unused + n_locked <= n_rows` for every row) and `test_unused_not_in_v6` (no id in `unused_candidates.csv` appears in the v6 table it is listed against).
- [ ] **Step 2: Run.** `.venv\Scripts\python.exe -m pytest tests/test_v7_inventory.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement and run.** Expected: PASS. Commit with `data(v7): local inventory of unused labelled data`.

**Part B: where the literature gets its data.**
- For every paper cited in `docs/v6_failure_analysis_and_ml_options.md` section 6 and `docs/literature_gaps_v6.md` that reports a model on fuels, mixtures or pure compounds, record one row in `docs/v7_data_sources.md` with these columns:
  - `paper, property, data source as stated by the authors, n (fuels / systems / compounds / points), composition method (GCxGC, GC-MS, PIONA, none), split type (random, group, family), access (open / on request / licensed / unknown), URL or DOI, overlap with our data (by fuel id, InChIKey or system), verified (full text / abstract only)`
- Priority papers:
  - Hall 2021 and Hall's dissertation (freezing: 74 conventional and 50 synthetic fuels)
  - Boehm 2022 and 2024
  - Bell 2025
  - Shi 2017, Shao 2024, Dong 2026
  - Leenhouts 2024
  - Dobbelaere 2022
  - Saldana 2013, Gaudin 2015, Paricaud 2019
  - Specht 2024 (HANNA, built on the Dortmund Data Bank)
  - Sun 2020, Lee 2026 (flash sets of 10k+ molecules)
  - the melting-point sets (Huang 2025, Dubonos 2026)
- Where a paper says "data from X", follow X to its primary database or programme and record that, not the paper. The DLR and AFRL fuel sets, the NJFCP fuels, the CRC reports and the Dortmund Data Bank are likely candidates; confirm each one.
- Rule: no number goes in this file unless it is quoted from a source read in this task. Mark abstract-only claims "abstract only", as the v6 reviews did.

**Part C: online deep research.**
- Invoke the deep-research skill in its literature-review mode with this question: *"What public or obtainable datasets give (a) real jet, diesel and SAF fuels with detailed hydrocarbon composition (GCxGC or equivalent) together with freezing point, low-temperature viscosity, density or flash point; (b) multicomponent hydrocarbon mixture measurements of viscosity, density and flash point for C7 to C20 alkanes, cycloalkanes and aromatics; (c) pure-compound freezing point, enthalpy of fusion, flash point, Antoine and LFL data for hydrocarbon isomers and oxygenates; (d) synthetic or computed mixture data (COSMO-RS, COSMO-SAC) used to pretrain mixture models? For each, give the size, access terms and licence."*
- Sources to cover:
  - Consensus
  - web search
  - data repositories: Zenodo, Figshare, Mendeley Data, NIST (ThermoML archive, TDE, WebBook), the NREL fuel databases, the AFRL/UDRI and DLR fuel programmes, the Co-Optima YSI and fuel-property databases, CRC, the ASTM research reports
  - the supplementary files of the Part B papers
- Results are appended to `docs/v7_data_sources.md` as a second table with the same columns, plus `already ingested (yes / partly / no)`.

**Part D: synthesis.** Write `docs/v7_data_plan.md`.
- For each property, give a ranked list of data to acquire.
- Each entry states: the gap it fills (by Task 0.5 Part A row, or by the Task 8 reading once it exists), the expected size in independent groups (fuels, systems or scaffolds, not rows), the access route and licence, the effort to ingest (a loader to the `measurements.csv` schema, a unit test with a real-structure fixture and an audit check, per AGENTS.md), and whether it could serve as a fresh lock (no overlap with any v6-scored row).
- Put the items Part A found on disk first, before anything that needs a download or a data request.
- [ ] Commit the two docs with `docs(v7): data sources survey and acquisition plan`. Then stop and ask the user which items to acquire, and which, if any, to reserve as a fresh lock (this feeds Task 20).

### Task 1: Lock guard

**Files:**
- Create: `src/fuelblend/v7/__init__.py`, `src/fuelblend/v7/guard.py`
- Test: `tests/test_v7_guard.py`

**Interfaces:**
- Produces: `class LockedRowError(ValueError)`. `LOCKS: dict[str, LockSpec]`, where `LockSpec(name, path, id_kind, scope)`. `id_kind` is one of `inchikey`, `system_id`, `group` or `fuel_id`. `scope` is the property, or the head, that the lock covers: `ysi_holdout` covers the YSI label only; each `mixture_lock` head covers that head's property only; `freezing_holdout` covers freezing; `flash_oxy_lock` and `mixture_holdout51` cover their property; `splits_locked` covers the v4/v5 real-fuel lock (`real_fuel_ids_locked`), all properties.
- `locked_ids(name: str) -> set[str]`. `assert_dev_only(ids, *, locks: tuple[str, ...]) -> None` raises on overlap with the named locks. `mask_locked(df, id_col: str, *, locks: tuple[str, ...], label_cols: list[str] | None = None) -> tuple[pd.DataFrame, int]`: with `label_cols` it sets only those labels to NaN (for example YSI for the YSI lock) and keeps the row for the other targets; without `label_cols` it drops the row. Every call names its locks explicitly; there is no global default.
- Sources: `results/v6/gnn/ysi_holdout/lock.json`, `results/v6/freezing/holdout_lock.json`, `results/v6/mixture/lock.json` (per head `locked`), `results/v6/mixture/holdout51_pred.csv` (its system ids), `results/v6/flash_oxy/lock.json`, `splits/locked.json`. Read each file's real schema before coding; do not guess it.

- [ ] **Step 1: Write the failing tests.**
  - `test_guard_refuses_each_lock_file`: parametrised over the six lock names; one real id from each file raises `LockedRowError` whose message contains the lock name.
  - `test_guard_passes_dev_ids`: one dev id per source passes.
  - `test_mask_label_only`: the YSI lock masks the `ysi` label of a locked compound and keeps its `flash_C` label.
  - `test_drop_counts`: the returned count equals the overlap size.
- [ ] **Step 2: Run.** `.venv\Scripts\python.exe -m pytest tests/test_v7_guard.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement `guard.py`.** Use plain JSON and CSV reads, cached with `functools.lru_cache`.
- [ ] **Step 4: Run.** Expected: PASS.
- [ ] **Step 5: Commit.** `feat(v7): lock guard refusing v6-scored rows`

---

## Phase 1: instrumentation and GPU environment

### Task 2: Curve logger and hooks in the training loops

**Files:**
- Create: `src/fuelblend/ml/curves.py`
- Modify: `src/fuelblend/ml/dmpnn.py:96-128` (`fit_ensemble`), `src/fuelblend/ml/molpool.py:67-97` (`fit_predict`), `src/fuelblend/ml/set_models.py:210-250` (`train_ensemble`)
- Test: `tests/test_v7_curves.py`

**Interfaces:**
- Step 0, before any edit: run the current `fit_ensemble`, `fit_predict` and `train_ensemble` on CPU at seed 0 on a 30-row fixture (30 Track B pool rows; for MolPool, 30 mixture rows from `results/v6/mixture/rows.csv`). Save the predictions as `tests/data/v7_golden_{dmpnn,molpool,set_models}.npy` and commit them.
- Produces: `class CurveLogger(csv_path: Path, *, run: dict[str, str], tb_dir: Path | None = None)`, with `.log(epoch: int, **metrics: float) -> None`, `.child(**run_extra) -> CurveLogger` (same file, extra id columns such as `fold` and `seed`) and `.close()`. CSV columns: `run_id, model, fold, seed, epoch, train_loss, val_loss, lr, sec, gpu_mem_mb`, plus any `val_<task>` keys. Append and flush every row. TensorBoard is used only when `tb_dir` is set and `torch.utils.tensorboard` imports.
- Produces: `lgb_curve_callback(logger: CurveLogger) -> Callable` for `lgb.train(callbacks=[...])`.
- All three loops also take `on_epoch: Callable[[int, float], None] | None = None`, which is called with `(epoch, val_loss)` after each epoch. The Optuna pruner raises `TrialPruned` through it; the exception must propagate.
- `dmpnn.fit_ensemble(..., logger=None, on_epoch=None, scheduler: str = "none", min_epochs: int = 0)`. `scheduler` is one of `none`, `plateau` (ReduceLROnPlateau factor 0.5, patience = patience // 3) or `cosine`. Per-task masked val MAE is logged as `val_<task>`. The best epoch per model is returned as the attribute `net.best_epoch_`.
- `molpool.fit_predict(..., device=None, val_idx=None, patience=None, logger=None)`. When `val_idx` is given, train on the rest, early-stop on val MSE and restore the best state. Defaults reproduce the current behaviour.
- `set_models.train_ensemble(..., curve=None)`. The existing `log` callable is unchanged.

- [ ] **Step 1: Write the failing tests.**
  - `test_fit_ensemble_logger_none_bitwise`: the modified `fit_ensemble` with default kwargs, on CPU at seed 0, gives predictions `np.array_equal` to `tests/data/v7_golden_dmpnn.npy`.
  - `test_molpool_defaults_bitwise` and `test_set_models_defaults_bitwise`: the same comparison against their golden files.
  - `test_on_epoch_prune_propagates`: an `on_epoch` that raises a custom exception at epoch 3 stops the fit with that exception.
  - `test_logger_rows_per_epoch`: with a logger, the CSV has one row per epoch run per model, and `epoch` is monotone within `(fold, seed)`.
  - `test_curve_csv_partial_read`: after 3 `.log` calls without `.close()`, `pd.read_csv` returns 3 rows.
  - `test_lgb_callback_records`: a 50-row synthetic LightGBM fit with `valid_sets` gives 20 rows for `num_boost_round=20`.
- [ ] **Step 2: Run.** `.venv\Scripts\python.exe -m pytest tests/test_v7_curves.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement `curves.py` and the hooks.** No behaviour change when the new kwargs are at their defaults.
- [ ] **Step 4: Run** this file, then `tests/test_v6_gnn.py` and `tests/test_v6_mixture.py`. Expected: PASS.
- [ ] **Step 5: Commit.** `feat(v7): per-epoch curve logging in torch and LightGBM loops`

### Task 3: GPU venv, device resolution and monitoring

**Files:**
- Create: `scripts/v7/gpu/setup_gpu_env.ps1`, `scripts/v7/gpu/check_gpu.py`, `scripts/v7/monitor.ps1`, `scripts/v7/plot_curves.py`
- Modify: `src/fuelblend/ml/dmpnn.py:14-15` (replace the module constant with `resolve_device()`, keep `DEVICE` as an alias), `requirements.txt` (add `optuna`, `tensorboard`; torch stays `>=2.1`). Install `optuna` and `tensorboard` into `.venv` too, so that every test runs in `.venv`. GPU-only scripts are never imported by tests.
- Test: `tests/test_v7_device.py`

**Interfaces:**
- Produces: `dmpnn.resolve_device(req: str | None = None) -> torch.device`, which reads `FUELBLEND_DEVICE`. If `cuda` is asked for and unavailable it raises `RuntimeError("CUDA requested but torch <ver> has no CUDA")`; it never falls back. `device_info() -> dict` returns `device`, `torch_version`, `cuda_version` and `gpu_name`. Every v7 results JSON embeds this dict.
- `setup_gpu_env.ps1` creates `C:\fuelblend-gpu\.venv` (Python 3.12), installs `requirements.txt` without torch, then `torch` from `https://download.pytorch.org/whl/cu128`, then `optuna` and `tensorboard`. Code runs from the repo with `$env:PYTHONPATH="<repo>\src"`. This replaces the `scratch/torch_gpu` PYTHONPATH workaround (torch 2.11+cu128 against the venv's 2.15+cpu).
- `monitor.ps1` runs `tensorboard --logdir $env:FUELBLEND_RUNS --port 6006`. `plot_curves.py <curves.csv> --out <png>` draws median and IQR train/val loss across seeds per fold and marks the best epoch. It also works for a run that is still going.

- [ ] **Step 1: Write the failing tests.** `test_device_resolution_reports_actual` (monkeypatch `torch.cuda.is_available` to False with `FUELBLEND_DEVICE=cuda`: expect `RuntimeError`; unset: expect `cpu`); `test_device_info_keys`.
- [ ] **Step 2: Run.** Expected: FAIL.
- [ ] **Step 3: Implement.** Then ask the user before running `setup_gpu_env.ps1`, because it downloads about 3 GB (torch cu128 plus dependencies). Ask the same way before the Chemprop install (Task 7) and the HANNA weights (Task 17).
- [ ] **Step 4: Verify.** Run the tests (PASS). Then run `C:\fuelblend-gpu\.venv\Scripts\python.exe scripts/v7/gpu/check_gpu.py`, which should print `cuda True NVIDIA GeForce RTX 3050 Laptop GPU` and a 1-epoch D-MPNN timing on CPU and on GPU.
- [ ] **Step 5: Commit.** `feat(v7): GPU venv setup, strict device resolution, curve monitor`

---

## Phase 2: cheap diagnostics (these set the Phase 3 budget)

All four tasks follow the same order: write `results/v7/<area>/PREDECLARATION.md` (comparators, rows, metric, reading rule) and commit it, then run and commit the results. Rows are the v6 dev rows with `guard.mask_locked` applied (locks named per task); the v6 comparator numbers are recomputed on the same reduced rows from v6 `oof.csv`.

### Task 4: v7 runner skeleton (shared by Tasks 5 to 13)

**Files:**
- Create: `src/fuelblend/v7/runner.py`
- Test: `tests/test_v7_runner.py`

**Interfaces:**
- Produces: `load_track_b_pool() -> tuple[pd.DataFrame, np.ndarray, list[str]]`, returning the pool, the folds and the task names. It loads the full v6 pool (`results/v6/gnn/pool.csv`) and computes folds on the full pool first, with the same generic-Murcko `GroupKFold(5)` as `scripts/v6/gnn/train_eval.py`. `oof.csv` has no fold column, so the folds are recomputed. Only then does it call `guard.mask_locked(..., locks=("ysi_holdout",), label_cols=["ysi"])`, so folds are never recomputed on a reduced pool. `scaffold_folds(df, k=5) -> np.ndarray` is exposed for the test. `write_results(out_dir, *, summary: pd.DataFrame, oof: pd.DataFrame, meta: dict)`, which adds `device_info()`, the git commit and `guard` drop counts. `paired_boot_gain(err_a, err_b, groups, n=4000, seed=0) -> tuple[float, float]` (90 % CI, group bootstrap, reusing `octane_blend_v4.cluster_boot`).

- [ ] Steps 1 to 5 as in Task 1. Tests: `test_folds_match_v6` (fold ids equal those recomputed by the v6 code on the full pool), `test_ysi_locked_label_masked`, `test_write_results_meta_keys`. Commit `feat(v7): shared runner for v7 experiments`.

**Order note:** Task 9 (the HPO harness) is built directly after Task 4, because Task 6 already uses it. Execute in the order 4, 9, 5, 6, 7, 8.

### Task 5: Epoch curves of the declared v1 D-MPNN config

**Files:** Create `scripts/v7/track_b/curves_v1.py`; output to `results/v7/track_b/curves_v1/` (`curves.csv`, `best_epoch.csv`, `curves.png`, `meta.json`).

**Interfaces:** Consumes `runner.*`, `fit_ensemble(logger=..., epochs=1000, patience=100)` and both v1 `CFGS` from `train_eval.py:22`.

- [ ] **Step 1:** Predeclare. Reading rule: if the median `best_epoch_` is below 0.5 x the v1 cap (75) and the val loss is flat or rising afterwards, longer training is not the limit. Report the v1-cap and long-cap outer MAE side by side.
- [ ] **Step 2:** Run on the GPU, in the background, with the log in `scratch/v7_curves_v1.log`. Record seconds per epoch per model in `meta.json`. `PHASE2_READINGS.md` turns this figure into projected hours and the `timeout_s` value for each Phase 3 task.
- [ ] **Step 3:** Commit the results and the reading. `exp(v7): epoch curves for v1 D-MPNN (dev only)`

### Task 6: Descriptor GBM and RF comparator for Track B

**Files:** Create `src/fuelblend/features/rdkit_desc.py` (`rdkit_2d(smiles: list[str]) -> pd.DataFrame`: all `Descriptors.descList`, NaN kept, constant columns dropped inside the fold) and `scripts/v7/track_b/descriptor_gbm.py`; output to `results/v7/track_b/descriptor_gbm/`. Test: `tests/test_v7_rdkit_desc.py` (`test_rdkit_2d_shape_and_names` on 3 SMILES, `test_invalid_smiles_row_nan`).

**Interfaces:** One LightGBM and one `RandomForestRegressor(n_estimators=500)` per target. Inner-fold Optuna is used for LightGBM only (30 trials, `num_leaves`, `min_child_samples`, `learning_rate`, `n_estimators` with early stopping on the inner val). Same outer folds as Task 4. The summary columns match `results/v6/gnn/cv_summary.csv`, plus `mae_gbm, mae_rf, gain_gnn_vs_gbm_lo, gain_gnn_vs_gbm_hi`.

- [ ] Predeclare (family: 10 targets, Holm) → test → implement → run on CPU → commit. Reading: "GNN earns its place on target t" only if the Holm-adjusted gain interval over the better of GBM and RF excludes 0.

### Task 7: Chemprop parity

**Files:** Create `scripts/v7/track_b/chemprop_parity.py`; output to `results/v7/track_b/chemprop_parity/`. Install `chemprop>=2` into the GPU venv only.

**Interfaces:** Two targets, `flash_C` and `ysi`, single-task, same outer folds, Chemprop defaults with matched depth 4 and hidden 200, a 5-model ensemble.

- [ ] Predeclare. Parity holds if the paired gain CI of in-house vs Chemprop contains 0 for both targets. If in-house is worse with the CI excluding 0, stop Phase 3 Track B tuning and open a debugging task. Run on the GPU, commit.

### Task 8: Sample-size learning curves

**Files:** Create `src/fuelblend/v7/sample_size.py` and `scripts/v7/sample_size_curves.py`; output to `results/v7/sample_size/<model>/` (`curve.csv` with `frac, repeat, fold, n_train_groups, mae`; `curve.png`). Test: `tests/test_v7_sample_size.py` (`test_subsample_is_group_level`, `test_fracs_monotone_nesting`: the groups at 0.2 are a subset of those at 0.4 for the same repeat).

**Interfaces:** `group_subsample(groups: np.ndarray, frac: float, rng) -> np.ndarray` (indices). Fractions are `0.1 0.2 0.4 0.6 0.8 1.0`, with 3 repeats; only the training groups of each outer fold are subsampled. Models: Track B D-MPNN (GPU, v1 config) and the Task 6 GBM, the v5 viscosity and density pair residuals (`viscosity_pair.LGB_PARAMS` and `density_pair.LGB_PARAMS` refit in v7 dirs on v5 dev rows), and the freezing LOGO ridge and residual (`scripts/v6/freezing/step3_models.py`).

- [ ] Test → implement → run → commit. In the reading, report the slope of MAE vs log(n) over the top two fractions with a bootstrap CI. A slope CI that excludes 0 means more data is expected to help; one that includes 0 means the curve has flattened. Give no forecast beyond the data.

**Phase 2 exit.** Run the full test suite. Write `results/v7/PHASE2_READINGS.md`, in which the Task 5 to 8 readings fix the Phase 3 trial budgets. The rule is predeclared: a model family whose Phase 2 result shows no headroom (Task 5 flat, Task 8 slope CI including 0, Task 6 comparator equal or better) gets the minimum budget (20 trials per outer fold); the others get 60.

---

## Phase 3: instrumented retrain with nested HPO

### Task 9: Optuna HPO harness

**Files:** Create `src/fuelblend/v7/hpo.py` and `configs/v7/search_spaces.yaml`. Test: `tests/test_v7_hpo.py`.

**Interfaces:**
- `nested_cv(df, groups, outer_folds, *, fit_score: Callable[[dict, idx_tr, idx_va, CurveLogger], float], space: dict, n_trials: int, timeout_s: int, study_prefix: str, storage: str) -> pd.DataFrame` gives one row per outer fold: `fold, best_params (json), inner_score, outer_pred_path`. The inner split is `GroupKFold(3)` on training groups only. `MedianPruner(n_warmup_steps=20)` prunes through the `on_epoch` hook (Task 2), reporting steps from inner fold 0 only. Each trial fits one model per inner fold; the 5-model ensemble is fitted only for the final outer-fold fit. The storage is `sqlite:///C:/fuelblend-gpu/runs/optuna.db`, and the study name is `{study_prefix}_fold{k}`, `load_if_exists=True`.
- Search spaces (YAML keys = model names):
  - `dmpnn`: hidden 100 to 600 (log), depth 2 to 6, dropout 0 to 0.4, lr 1e-4 to 3e-3 (log), wd 1e-6 to 1e-3 (log), bs one of 32, 64 or 128, scheduler one of none, plateau or cosine; epochs 1000, patience 100.
  - `molpool`: d 4 to 32, rank 2 to 16, hidden 16 to 128, wd 1e-3 to 1 (log), lr 3e-4 to 1e-2 (log), dropout 0 to 0.3; epochs 2000, patience 150.
  - `lgbm`: num_leaves 3 to 31, min_child_samples 5 to 60, learning_rate 0.01 to 0.2 (log), subsample 0.5 to 1, colsample 0.5 to 1, lambda_l2 1e-3 to 10 (log); n_estimators 5000 with 100-round early stopping.
- [ ] **Step 1: Write the failing tests.** `test_study_resume_no_duplicate_trials` (run 4 trials, run again with `n_trials=6`: storage has 6 trials, not 10); `test_outer_never_in_inner` (an assertion inside `fit_score` that no outer test index appears in the train or val indices); `test_outer_scored_once` (an outer predictions file already present makes the call skip that fold).
- [ ] Steps 2 to 5 as usual. Commit `feat(v7): nested Optuna harness with resumable studies`.

### Task 10: Track B D-MPNN retrain on the GPU

**Files:** Create `scripts/v7/track_b/retrain.py`; output to `results/v7/track_b/retrain/` (`cv_summary.csv` in the v6 columns plus `mae_gbm` from Task 6, `oof.csv`, `best_params.csv`, `curves.csv`, `curves.png`, `meta.json`). TensorBoard logs go to `$FUELBLEND_RUNS/track_b/`.

**Interfaces:** Consumes `hpo.nested_cv` with `space["dmpnn"]` and the Phase 2 budget. The final fit per outer fold is a 5-model ensemble with the chosen params. On CUDA OOM, halve `bs` once and log `oom_retry=1`. This is implemented as `fit_with_oom_retry(fit_fn, **kw)` in `runner.py`, with the test `test_oom_retry_halves_bs_once` (a stub that raises `torch.cuda.OutOfMemoryError` on the first call: the second call gets `bs/2`, a third OOM re-raises).

- [ ] Predeclare the comparators (v6 v1 OOF on the same rows, the Task 6 GBM, the constant; Holm over 10 targets). Run in the background:

```powershell
$env:FUELBLEND_DEVICE="cuda"; $env:PYTHONPATH="$PWD\src"; C:\fuelblend-gpu\.venv\Scripts\python.exe scripts/v7/track_b/retrain.py *> scratch\v7_b_retrain.log
```

  Monitor with `scripts/v7/monitor.ps1`. Commit the results; v1 stays the declared v6 result.

### Task 11: Mixture encoder, MolPool heads and flash oxygenate arm

**Files:** Create `scripts/v7/mixture/retrain.py` and `scripts/v7/flash_oxy/retrain.py`; output to `results/v7/mixture/` and `results/v7/flash_oxy/`. The encoder pretrain reuses `scripts/v6/mixture/pretrain_encoder.py` logic through an import (no copy), with `epochs=1000, patience=100, logger` set, and saves to `results/v7/mixture/encoder.pt`.

**Interfaces:** Consumes `molpool.fit_predict(device=..., val_idx=..., patience=150, logger=...)`, `hpo.nested_cv(space["molpool"])` and `guard.mask_locked(df, "system_id", locks=("mixture_lock", "mixture_holdout51", "flash_oxy_lock"))`. The heads are the v6 heads from `run_cv.py`, and the comparators are the v6 physics and descriptor GBM baselines on the same rows. Dev OOF only.

- [ ] Predeclare → run on the GPU → commit. In the README, state that the rows are the v6 dev rows (already seen in v6 selection), so any gain is exploratory and cannot reopen the Track C verdict.

### Task 12: v5 set models (deep sets and pairwise)

**Files:** Create `scripts/v7/mixture/retrain_set_models.py`; output to `results/v7/set_models/`.

**Interfaces:** Consumes `set_models.train_ensemble(curve=..., device="cuda", epochs=500)` with a 20 % group-held-out val split for early stopping (added in the script around `train_ensemble`; the library code is unchanged beyond Task 2). Comparators are the v5 numbers read from the BASE files the v5 scripts wrote (read only).

- [ ] Predeclare → run → commit.

### Task 13: CPU tree models with curves and Optuna

**Files:** Create `scripts/v7/cpu/retrain_trees.py`; output to `results/v7/retrain_v5/<model>/` for `viscosity_pair`, `density_pair`, `fusion_estimate`, `ysi_lgbm` (from `delta_ml_ysi.py`), `freezing_step3` and `mixture_gbr`.

**Interfaces:** Two steps per model.
1. Reproduction gate: refit with the original parameters (for example `viscosity_pair.LGB_PARAMS`) in the v7 directory, on the same rows and groups as the v5 or v6 script that produced its BASE or v6 number. The MAE must match the number in that file to 1e-6 relative, recorded in `repro_gate.json`. If the gate fails, stop for that model and report why.
2. Refit through `lgb.train` with `lgb_curve_callback` (or `staged_predict` for sklearn GBR) under `hpo.nested_cv(space["lgbm"])`. The paired gain CI is computed against the per-row predictions of the step 1 refit.

Outputs: `repro_gate.json`, `summary.csv` (v7 MAE, refit-original MAE, the BASE or v6 MAE read from file, paired gain CI), `curves.csv`, `curves.png`, `best_params.json`. The script asserts that no file under `src/` or `results/v5/BASE/` changed (by hashing before and after).

- [ ] Predeclare → run (CPU, background) → commit. Adoption into `predict_fuel()` is a separate decision for the user and is not part of this plan.

**Phase 3 exit.** Run the full suite, then write `results/v7/PHASE3_READINGS.md`: per model, the v7 vs comparator gain CI and the curve reading (where the best epoch sits, and the train/val gap at the best epoch).

---

## Phase 4: literature methods

### Task 14: Freezing limiting-species baseline (d)

**Files:** Create `scripts/v7/freezing/baseline_d.py`, importing helpers from `scripts/v6/freezing/common.py` if needed. Never modify `scripts/v6/*`. Output to `results/v7/freezing/baseline_d/`.

**Interfaces:** For each fuel, the predicted freezing point is the max over species of the ideal-solubility SvL temperature `T_i` solving `ln x_i = -(dHfus_i/R)(1/T - 1/Tm_i)`. Dev LOGO rows (52 groups) only. Report MAE, bias, and one-sided over-prediction beside baselines (a) to (c) from `results/v6/freezing/dev_logo_summary.csv`.
- [ ] Test `test_limiting_species_single_component` (a pure-compound "fuel" returns its Tm) and `test_limiting_species_is_max` → implement → run → commit.

### Task 15: Freezing by family, heavy-tail GP and novelty flag

**Files:** Create `results/v7/freezing/families.csv` (`group, family`, assigned from fuel metadata before any scoring and committed first), `src/fuelblend/v7/freeze_gp.py` and `scripts/v7/freezing/lofo_gp.py`. Test: `tests/test_v7_freeze_gp.py`.

**Interfaces:** `heavy_tail_features(comp_df) -> pd.DataFrame` gives `nC14plus_mass`, `naphthalene_mass`, `dicyclo_mass`, `tm_max_present` and `tm_p90_mass`. `fit_gp(X, y) -> GaussianProcessRegressor` uses an ARD Matern 2.5 kernel plus white noise. `novelty_score(X_tr, X_te) -> np.ndarray` is the Mahalanobis distance in the standardised feature space, with the flag threshold at the 95th percentile of the training LOFO scores. Evaluation: leave-one-family-out and LOGO on the 52 dev groups. Coverage of the live 13.84 C margin is reported by family.
- [ ] Tests: `test_heavy_tail_features_sum_bounds`, `test_novelty_flags_far_point`. Predeclare (comparators: constant, ridge, UNIFAC-SvL, baseline d; family-level reading) → implement → run → commit. In the README, freezing stays "status: no skill shown" unless the gate passes on a new lock (Task 20).

### Task 16: Joint pure plus fuel MolPool for freezing (GPU)

**Files:** Create `scripts/v7/freezing/molpool_joint.py`; output to `results/v7/freezing/molpool_joint/`.

**Interfaces:** Consumes the Task 11 encoder, `molpool.fit_predict` (with the Task 2 kwargs) over the GC x GC species list per fuel, and an auxiliary pure-Tm loss (multi-task: the fuel freezing label plus pure Tm rows as single-species "fuels"). MC-dropout (30 passes) gives the uncertainty. LOFO and LOGO on dev, with the same comparators as Task 15.
- [ ] Predeclare → run → commit.

### Task 17: HANNA gamma vs UNIFAC in the v5 flash route

**Files:** Create `src/fuelblend/thermo/gamma_hanna.py` (a wrapper with the same call signature as the UNIFAC gamma used by the v5 flash route) and `scripts/v7/flash/hanna_vs_unifac.py`. Test: `tests/test_v7_gamma_hanna.py`.

**Interfaces:**
- Step 0, licence: confirm the licence of the HANNA code and weights and their source URL. Ask the user before downloading. If the licence forbids research use or redistribution, stop and record that in the README.
- `gamma_hanna(smiles: list[str], x: np.ndarray, T_K: float) -> np.ndarray`, binary only; it raises `NotImplementedError` for more than 2 components (HANNA is binary-only, failure doc section 4b). Pretrained, no refit.
- Rows: Leenhouts hydrocarbon binary flash rows from `measurements_dedup.csv` (dev only, `guard` applied). Metrics: MAE and one-sided over-prediction of the v5 flash route.
- Three arms. Production flash is ideal (`predict_route.py:233`, `FLASH_METHOD = "ideal_le_chatelier_anchored"`, `use_unifac=False`). The arms are therefore: ideal (gamma = 1, the production comparator), UNIFAC (`flash_point_unifac_detail(use_unifac=True)`) and HANNA, the last through the same function with the gamma source swapped.
- [ ] Tests: `test_hanna_pure_limit_gamma_one` (x to 1 gives gamma 1 within 1e-3), `test_hanna_rejects_ternary`. Predeclare (primary comparison: HANNA vs ideal; secondary: HANNA vs UNIFAC; Holm over the two; gate: non-inferior on over-prediction with the upper CI bound below 0.5 C, and the MAE gain CI above 0) → implement → run → commit. Not used in freezing.

### Task 18: Equation-parameter heads (Antoine and Andrade)

**Files:** Create `src/fuelblend/ml/param_heads.py` and `scripts/v7/track_b/param_heads.py`; output to `results/v7/track_b/param_heads/`. Test: `tests/test_v7_param_heads.py`.

**Interfaces:** `class ParamHeadDMPNN(DMPNN)` outputs `(A, B, C)` for Antoine (log10 P/kPa = A - B/(T_K + C)) and `(a, b)` for Andrade (ln mu = a + b/T_K). The loss is computed on measured `(T, value)` pairs from `measurements_dedup.csv` pure rows (P, viscosity), so one molecule contributes many points. C is constrained by a softplus offset to keep `T_K + C > 0` over the data range. Evaluation: scaffold GroupKFold, MAE in log P and log mu at held-out T, and flash recomputed through the v5 route from the predicted Antoine constants for hydrocarbons with measured flash. Respect Trap 3: no evaluation of un-refit Antoine for C10+ alkanes in the 40 to 100 C window.
- [ ] Tests: `test_antoine_head_shape`, `test_andrade_monotone_in_T` (b > 0 gives mu falling with T). Predeclare → implement → run on the GPU → commit.

### Task 19: Family-conditional and shift-weighted conformal

**Files:** Create `src/fuelblend/uq/conformal_v7.py` and `scripts/v7/uq/coverage_by_family.py`; output to `results/v7/uq/`. Test: `tests/test_v7_conformal.py`.

**Interfaces:** `one_sided_margin(resid, alpha=0.10, weights=None) -> float` gives the weighted (1 - alpha) quantile with the finite-sample correction. `group_conditional_margins(resid, groups, alpha) -> dict[str, float]` uses Mondrian margins per family (with a fallback to the global margin when a family has fewer than 10 calibration rows). `shift_weights(X_cal, X_test) -> np.ndarray` uses a logistic-regression density ratio, clipped to [0.05, 20]. Applied on dev to freezing (Task 15 families), YSI (scaffold families) and flash. The report gives coverage by family and states the number of groups behind each coverage figure.
- [ ] Tests: `test_margin_exact_on_uniform` (residuals 1 to 100, alpha 0.10: the margin equals the 91st value), `test_mondrian_fallback_small_group`, `test_weights_clipped`. Implement → run → commit.

### Task 20: Fresh lock (user decision)

**Files:** Create `scripts/v7/lock/check_external_overlap.py` and `results/v7/lock/overlap.json`.

**Interfaces:** `data/interim/external_staging/canonical_leenhouts_flash_mix_measurements.csv` mapped to `system_id`, then checked for overlap with `measurements_dedup.csv` systems and `guard.locked_ids(name)` for every name in `guard.LOCKS`. Report the counts of systems and rows not seen anywhere.
- [ ] Run and commit the overlap report. **Stop and ask the user** whether to spend the unseen systems as the v7 flash lock and which single primary comparison it tests (for example Task 17 HANNA vs UNIFAC). If yes, write `results/v7/lock/LOCK.md` and `lock.json` with the comparison family and the Holm rule, commit them before any fit, and score once.

---

## Phase 5: documentation

### Task 21: Protocol and doc fixes

**Files:** Modify `docs/v6_summary.md` (header link to both reviews: fix 3.6), `docs/literature_review_v6.md` (Leenhouts and Hall years flagged as online vs print, pending a DOI check: fix 3.3), `docs/v6_negative_results.md` (add v7 dev-only rows where a v7 model shows no skill). Create `docs/v7_protocol.md` (multiplicity rule: fix 3.5; the lock guard; the CPU/GPU separation; dev-only status of all Phase 2 to 4 results) and `docs/v7_summary.md` (from the `PHASE*_READINGS.md` files, every number traced to a `results/v7/` file).
- [ ] Mark fixes 3.1 and 3.2 as already applied. Run the full suite and record the counts. Commit `docs(v7): protocol, summary, and v6 doc fixes`.

---

## Deferred (start only after Task 17 and with a new predeclaration)

- COSMO-SAC or ML-generated mixture pretraining for MolPool (F2 fix 2): needs a COSMO-SAC implementation and sigma profiles for all components.
- CheMeleon descriptor-pretrained encoder (gap 3): check the licence first.
- A tabular foundation model on descriptors as the Track C comparator (gap 7): the source is non-peer-reviewed.
- GE-GNN gamma plus predicted Antoine and LFL as the flash oxygenate physics comparator (F4, gap 6): depends on Tasks 17 and 18.
