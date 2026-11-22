# v9 readings

Rules: `docs/v9_predeclaration.md` (verbatim). Data added in v9 (`results/v9/inventory/delta.csv`): 3 mixture systems with density or vapour pressure, 4 fuel viscosity cells at 233 K and 6 at 253 K (Boehm 2022 neat fuels), 316 groups merged into existing v5 keys. The Edwards 2017 freezing drafts and the pure-constant drafts could not be verified and are excluded (`results/v9/verification.md`).

## Dev (data new in v9, not in the v9 lock, no id in any earlier lock)

| reading | n | verdict |
|---|---|---|
| R9-1 real-fuel freezing | 0 new fuels | status: no skill shown. No freezing value or composition was added; a dev-only gain could not reopen it anyway. |
| R9-2 viscosity residual with more data | 0 new composition-bearing systems | more data helps: not shown. The retrain was not run: the only new viscosity systems are real fuels without composition, so the training set equals the v5 one. |
| R9-3 real-fuel cold viscosity | 2 fuels with composition (fuel:143, fuel:24570) | status: no skill shown (floor of 8 fuels unmet; descriptive only). |
| R9-4 v5 routes on new mixture systems | 20 new systems, 407 groups (ThermoML 2020 files) | descriptive; no model scored (density, vapour pressure and viscosity groups only; no flash). No verdict change. |
