# v9 predeclaration (fixed before any v9 data is gathered or read)

## User decisions (recorded 2026-11-23, before any v9 work)

* Targets: all four (real-fuel freezing, mixture systems, pure constants, real-fuel cold viscosity).
* Licences: open and government sources plus numeric facts transcribed with citation from journal supporting information (Vozka 2018 practice). Raw PDFs and licensed files are never committed.
* Web access: fetching pages and downloading papers, SI and the ThermoML archive is allowed without asking.
* Re-evaluation: lock a share of the new data before looking, re-read models on dev, score the new lock once.
* Git: commit after each task; push to `origin feat/external-datasets-ingest` once, after Task 9.
* Execution: native, no subagents; Gemini 3.8 Flash for gathering and drafting, Sonnet for verification, locks, readings and push (user decision 2026-11-23).

## Lock rule

A new system or fuel (a group with no v5 row, `v9_new`) is in the v9 lock when `int(md5(group_id).hexdigest(), 16) % 4 == 0` (`fuelblend.v9.lock.is_v9_locked`). `results/v9/lock.json` is written once by Task 8 and refused if it exists. Only data new in v9 can enter the v9 lock, so no v5, v6 or v7 lock is re-read.

## Readings and pass rules

Statistics: paired group bootstrap over fuels or systems, 90 % CI.

| reading | pass rule |
|---|---|
| R9-1 real-fuel freezing | skill shown only if the CI of the gain (constant MAE minus model MAE) excludes 0 on the v9 lock. A dev-only gain cannot reopen "status: no skill shown". If the v9 lock holds fewer than 8 freezing fuels with carbon split (8 = the v5 locked freezing set, `docs/final_numbers_v5.md` section A), the lock reading is descriptive only and the status stays "no skill shown" |
| R9-2 viscosity residual with more data | retrained model vs the v5 `ml_pair_residual`, both out-of-fold on the v5 dev systems: "more data helps" only if the CI of the MAE gain excludes 0; then confirmed on the v9 lock with the same rule |
| R9-3 real-fuel cold viscosity (253.15 K, 233.15 K) | as R9-1, per temperature, same lock floor of 8 fuels |
| R9-4 v5 routes on new mixture systems | descriptive: MAE and bias per property vs the v5 lab figures; no verdict change |
