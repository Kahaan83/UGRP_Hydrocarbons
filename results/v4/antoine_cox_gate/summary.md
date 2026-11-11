# Antoine reference gate: flash-point effect (2026-11-09)

Antoine sets changed by the gate: 6 (CCCCCCCCCCCC, CCCCCCCCCCCCCC, CCCCCCCCCCCCCCCC, CCCCCCCCCCCCCCCCCCCC, CCCCCC, CCCCCCC).
Error = predicted − measured (°C); positive = over-prediction (the unsafe direction). LFL fallback 1 % as in the R2 audit.

## 126 mixtures, all

| registry | n | MAE | bias | RMSE | over_5C_pct | max_over |
|---|---|---|---|---|---|---|
| before | 126.0 | 6.31 | 1.25 | 10.14 | 18.3 | 41.49 |
| after | 126.0 | 6.54 | 0.95 | 10.26 | 18.3 | 41.49 |
| anchored | 126.0 | 3.74 | -3.3 | 4.55 | 0.8 | 10.6 |

## 126 mixtures, only those containing a changed compound

| registry | n | MAE | bias | RMSE | over_5C_pct | max_over |
|---|---|---|---|---|---|---|
| before | 99.0 | 3.73 | -2.64 | 5.06 | 2.0 | 18.55 |
| after | 99.0 | 4.03 | -3.03 | 5.35 | 2.0 | 17.38 |
| anchored | 99.0 | 3.78 | -3.39 | 4.7 | 1.0 | 10.6 |

## pure flash records of the changed compounds

| registry | n | MAE | bias | RMSE | over_5C_pct | max_over |
|---|---|---|---|---|---|---|
| before | 6.0 | 3.64 | -1.91 | 4.19 | 0.0 | 3.77 |
| after | 6.0 | 3.27 | -2.22 | 4.65 | 0.0 | 1.91 |
| anchored | 6.0 | 5.58 | -4.51 | 7.53 | 0.0 | 2.11 |
