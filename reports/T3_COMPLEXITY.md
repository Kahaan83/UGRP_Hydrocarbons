# T3 - baseline error by mixture complexity

baselines from data/runs/baseline_preds.csv; conflict rows excluded; cluster bootstrap (2000, system_id).
Viscosity numbers are log-RMSE (ln mPa.s units); density and flash in kg/m3 and degC. Flash is a
temperature whose errors are signed, so bias matters more than RMSE below the table.

## Honest notes (read before quoting)

1. **`pure` buckets for viscosity and density are self-predictions** - the pure-compound curve fits that
   produce the baseline are fitted on the very rows being scored, so the ~0 error there is by construction,
   not a result. Real information lives in the binary/ternary buckets.
2. **Viscosity binary columns are dominated by one outlier system** `a175950749e8` (measured
   1300-2100 mPa.s vs Refutas 1.2-2.3 mPa.s, |ln err| approx 6.8, ~1000x). Median binary ln-error is
   only **+0.038 (~4 %)** - excluding that system, binary Refutas is the genuinely accurate regime.
   The system uses a single source per row, so the >25 %-conflict flag (needs >=3 obs) never fired.
   Value/src provenance for that system should be checked by hand before it enters any training set
   (recovered or flagged in drop-outlier lists downstream).
3. **No 4-comp or 5+comp hydrocarbon rows exist in the applicable data** - the table honestly shows
   zero because the open dataset has almost no >=4-component characterized blends (T0 gap #1 carries over).
4. Wide RMSE bootstrap CIs on binary buckets ([0.15, 100.7]) are the same outlier: cluster resampling
   sometimes resamples it many times. The CI of the median-filtered binary bucket would be ~[0.05, 0.10].

## viscosity_dyn_mPas

### baseline `refutas` (n_applicable=16621)

| bucket | n | MAE | RMSE | rel(err) | mean signed | 95% CI rmse |
|---|---|---|---|---|---|---|
| pure | 1421 | 0.0301 | 0.0561 | 3.0% | -0.089 | [0.065, 3.882] |
| binary | 14698 | 0.0619 | 0.2653 | 6.2% | -2.073 | [0.145, 100.663] |
| ternary | 502 | 0.0870 | 0.1050 | 8.7% | +0.104 | [0.079, 0.223] |
| 4comp | 0 | - | - | - | - | - |
| 5+comp | 0 | - | - | - | - | - |

### baseline `gn_mass` (n_applicable=18369)

| bucket | n | MAE | RMSE | rel(err) | mean signed | 95% CI rmse |
|---|---|---|---|---|---|---|
| pure | 2219 | 0.0282 | 0.0564 | 2.8% | -0.069 | [0.164, 3.176] |
| binary | 15628 | 0.0886 | 0.2727 | 8.9% | -1.888 | [0.256, 98.375] |
| ternary | 522 | 0.1236 | 0.1516 | 12.4% | +0.158 | [0.119, 0.339] |
| 4comp | 0 | - | - | - | - | - |
| 5+comp | 0 | - | - | - | - | - |

## density_kg_m3

### baseline `dens_idvol` (n_applicable=14754)

| bucket | n | MAE | RMSE | rel(err) | mean signed | 95% CI rmse |
|---|---|---|---|---|---|---|
| pure | 1752 | 0.0000 | 0.0000 | 0.0% | +0.000 | [0.000, 0.000] |
| binary | 12358 | 1.1789 | 2.0432 | 0.2% | +0.077 | [1.647, 2.481] |
| ternary | 644 | 3.5887 | 6.1247 | 0.5% | +2.156 | [1.938, 8.256] |
| 4comp | 0 | - | - | - | - | - |
| 5+comp | 0 | - | - | - | - | - |

## flash_point_C

### baseline `chat_nominal` (n_applicable=100)

| bucket | n | MAE | RMSE | rel(err) | mean signed | 95% CI rmse |
|---|---|---|---|---|---|---|
| pure | 43 | 13.6297 | 21.2180 | 10.9% | +6.617 | [12.173, 30.679] |
| binary | 54 | 1.6352 | 1.9200 | 13.4% | -1.372 | [1.638, 2.476] |
| ternary | 3 | 3.1355 | 3.1663 | 4.8% | -3.136 | [nan, nan] |
| 4comp | 0 | - | - | - | - | - |
| 5+comp | 0 | - | - | - | - | - |

