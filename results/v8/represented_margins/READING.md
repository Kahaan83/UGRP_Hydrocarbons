# v8 G2 reading: freezing-only check (represented blendstock route, real fuels)

Scope: freezing only, as the go/no-go for Task 3. Route: designated UNIFAC Schroeder-van Laar (`unifac_desig`); constant: leave-one-group-out constant `C`. Fuels: `results/v5/BASE/delta_freeze_t2_unifac/rows.csv` (AJFTD, AV-33-22 and FuelLib entities with a finite route value; 58 fuels, each its own group; the committed pooled reading uses 50, the headline role only). Script: `scripts/v8/g2_freezing_check.py`; output `freezing_g2.csv`.

| quantity | value |
|---|---|
| MAE route / constant | 8.55 C / 6.24 C |
| bias of route (pred - meas) | +5.66 C (safe side) |
| group-mean gain of route over constant | -2.30 C, 90 % CI (-4.09, -0.47) |
| condition (a) beats constant, CI excludes 0 | not met (CI is entirely on the constant's side) |
| leave-source-out coverage AJFTD to others / others to AJFTD | 1.00 / 0.76 |

Reading: status: no skill shown. Condition (a) fails, so G2 fails for freezing; conditions (b) and (c) are not decisive (coverage of 1.00 is uninformative under safe-side bias). No v8 freezing margin is written; `configs/uncertainty_v8.yaml` is not created. Following the decision recorded for this run, Task 3 is not built for the other properties.
