# v9 readings

Rules: `docs/v9_predeclaration.md` (verbatim). Data added in v9 (`results/v9/inventory/delta.csv`): 3 mixture systems with density or vapour pressure, 4 fuel viscosity cells at 233 K and 6 at 253 K (Boehm 2022 neat fuels), 316 groups merged into existing v5 keys. The Edwards 2017 freezing drafts and the pure-constant drafts could not be verified and are excluded (`results/v9/verification.md`).

## Dev (data new in v9, not in the v9 lock, no id in any earlier lock)

| reading | n | verdict |
|---|---|---|
| R9-1 real-fuel freezing | 0 new fuels | status: no skill shown. No freezing value or composition was added; a dev-only gain could not reopen it anyway. |
| R9-2 viscosity residual with more data | 0 new composition-bearing systems | more data helps: not shown. The retrain was not run: the only new viscosity systems are real fuels without composition, so the training set equals the v5 one. |
| R9-3 real-fuel cold viscosity | 2 fuels with composition (fuel:143, fuel:24570) | status: no skill shown (floor of 8 fuels unmet; descriptive only). |
| R9-4 v5 routes on new mixture systems | 20 systems (3 of them new to v5), 407 new groups (ThermoML 2020 files) | descriptive; no model scored (density, vapour pressure and viscosity groups only; no flash). No verdict change. |

## v9 lock readings (scored once)

| reading | n | floor | verdict |
|---|---|---|---|
| R9-1 | 0 | 8 fuels | status: no skill shown (floor unmet, descriptive only); no verified freezing fuel was added |
| R9-2 | 0 | CI gain > 0 | more data helps: not shown (retrain not run); no new composition-bearing viscosity system |
| R9-3 | 0 | 8 fuels | status: no skill shown (floor unmet, descriptive only) |
| R9-4 | 3 | descriptive | descriptive, no verdict change (105 groups) |

Per-row files: `R9-*_v9_lock.csv`. The lock reading was run once.
