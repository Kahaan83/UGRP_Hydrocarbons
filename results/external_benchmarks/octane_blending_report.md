# KAUST octane blending benchmark (computed)

All numbers below are produced by `scripts/external_ingest/eval_octane_blending_benchmarks.py`.
Data: KAUST Table 4, 123 rows, of which 3 single-component anchor rows (TMB+0, CP+0, Hex+0) define the pure values and are excluded from the metrics; 120 blends scored. Validation: GroupKFold, 5 splits, groups = component set with ethanol removed (8 groups, so ethanol blends share a fold with their base fuel).

## Pure values used (from measured rows)

| Component | RON | MON | Source | Density kg/m3 | Density source |
|---|---|---|---|---|---|
| n-heptane | 0.0 | 0.0 | PRF scale definition | 683.8 | data/unified (median, 292-294.5 K) |
| i-octane | 100.0 | 100.0 | Table 3 row 21 | 691.9 | data/unified (median, 292-294.5 K) |
| toluene | 118.0 | 100.3 | Table 3 row 101 | 866.8 | data/unified (median, 292-294.5 K) |
| ethanol | 108.0 | 90.0 | Table 3 row 128 | 789.3 | AI4fuel CRC handbook (15-30 C) |
| 1,2,4-trimethylbenzene | 109.5 | 108.0 | Table 4 row TMB+0 | 875.7 | data/unified (median, 292-294.5 K) |
| cyclopentane | 100.0 | 85.6 | Table 4 row CP+0 | 744.5 | data/unified (median, 292-294.5 K) |
| 1-hexene | 73.6 | 64.5 | Table 4 row Hex+0 | 673.4 | data/unified (median, 292-294.5 K) |
| diisobutylene | 100.6 | 86.5 | Table 3 row 72 | 715.0 | AI4fuel CRC handbook (15-30 C) |

Diisobutylene uses the Table 3 row 72 isomer (2,4,4-trimethyl-1-pentene) as a proxy; no pure diisobutylene row exists. n-Heptane uses the PRF scale definition (not tabulated in Table 3).

## Overall metrics (ON)

| Target | Model | n | MAE | RMSE | Bias | Max |
|---|---|---|---|---|---|---|
| RON | Constant (train mean) | 120 | 12.476 | 14.825 | +1.541 | 49.38 |
| RON | Linear volume | 120 | 5.748 | 7.986 | -5.127 | 23.9 |
| RON | Linear molar | 120 | 2.791 | 3.783 | -1.141 | 13.87 |
| RON | Linear molar + GBM residual (GroupKFold) | 120 | 2.453 | 3.161 | +0.505 | 11.376 |
| MON | Constant (train mean) | 120 | 10.021 | 12.113 | +1.140 | 43.544 |
| MON | Linear volume | 120 | 5.644 | 7.751 | -4.767 | 23.6 |
| MON | Linear molar | 120 | 3.101 | 4.364 | -2.303 | 12.799 |
| MON | Linear molar + GBM residual (GroupKFold) | 120 | 2.092 | 2.815 | +0.133 | 10.156 |

## Skill checks (computed)

- RON: linear molar MAE 2.791 vs constant 12.476 -> beats the constant; GBM residual MAE 2.453 vs linear molar 2.791 -> improves on it; linear volume MAE 5.748 vs linear molar 2.791.
- MON: linear molar MAE 3.101 vs constant 10.021 -> beats the constant; GBM residual MAE 2.092 vs linear molar 3.101 -> improves on it; linear volume MAE 5.644 vs linear molar 3.101.

## Per base-fuel family (MAE, ON)

| Group | n | RON lin. mol | RON lin. vol | RON GBM | MON lin. mol | MON lin. vol | MON GBM |
|---|---|---|---|---|---|---|---|
| `1,2,4-trimethylbenzene` | 4 | 3.786 | 4.044 | 5.517 | 3.337 | 6.425 | 3.413 |
| `1-hexene` | 4 | 1.65 | 6.89 | 2.836 | 0.395 | 3.519 | 0.719 |
| `cyclopentane` | 4 | 0.606 | 0.15 | 2.336 | 1.001 | 0.585 | 1.301 |
| `diisobutylene+i-octane+n-heptane` | 10 | 5.085 | 4.524 | 4.736 | 4.195 | 3.657 | 4.236 |
| `diisobutylene+i-octane+n-heptane+toluene` | 16 | 3.804 | 4.246 | 3.558 | 4.001 | 4.241 | 3.576 |
| `i-octane+n-heptane` | 14 | 5.122 | 12.257 | 1.696 | 6.709 | 11.429 | 2.136 |
| `i-octane+n-heptane+toluene` | 61 | 1.938 | 5.194 | 1.968 | 2.351 | 5.195 | 1.584 |
| `n-heptane+toluene` | 7 | 1.296 | 6.265 | 0.502 | 1.408 | 7.684 | 0.446 |

## Density sensitivity (linear molar MAE, diisobutylene density scaled)

| DIB density scale | RON MAE | MON MAE |
|---|---|---|
| 0.95 | 2.804 | 3.11 |
| 1.05 | 2.778 | 3.091 |

## FACE gasoline + ethanol (Table 5), linear volume between measured FACE+0 and ethanol

n = 24 blends; RON MAE 7.044 (bias -7.044, max 13.100); MON MAE 4.960 (bias -4.960, max 11.050). No model is fitted here.

## Scope

Model-level benchmark on literature-compiled data. Oxygenate (ethanol) rows are in scope for octane only. No comparison to literature model errors is made because none is computed in this repository.
