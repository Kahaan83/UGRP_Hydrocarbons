# BASELINES.md - WP3 physics baselines on the unified data

## What ran - exact commands + seeds

```
python scripts/baselines_eval.py --eval
```

- script: `scripts/baselines_eval.py`; outputs `data/runs/baseline_preds.csv`, `data/runs/BASELINES_EVAL_LOG.md`, this file.
- seed 20261030 everywhere (splits, bootstrap, null-fold hashing); 2000 cluster-bootstrap resamples by `system_id` (`src/fuelblend/eval/metrics.py`); wall time 114 s.
- data: `data/processed/model_table.csv` (42,328 hydrocarbon-only de-duplicated records, built by `src/fuelblend/data/model_table.py`); component constants from `data/unified/compounds.csv`; splits `splits/s1s2s3.json`.
- nothing trained, nothing tuned; locked sets (`splits/locked.json`) never evaluated (asserted absent; 0 locked multi-component rows are demoted out of every split before scoring).
- **Mode A only**: a record is predicted iff every one of its components has pure-compound data of the needed property covering the record temperature (2 K pad). No novel-molecule extrapolation is attempted.
- pure-component curves are fitted from the pure rows of the same table: `ln eta = a + b ln T` (viscosity, in-range only) and `rho = a + b T` (density). A single-measurement compound is a constant. The viscosity fit starts from a Theil-Sen line and drops points more than a factor 3 off it (2 of 2199 pure viscosity points dropped over 433 compounds).
- null model: median of (`property`, `kind`, 20 K T-bin) from the **S1 training rows**; S1 rows are scored against a 5-fold system-disjoint median so the null is never in-sample. Same rows are used for the physics baseline and the null in every cell.
- default filter `conflict=False` for every number below unless the line says otherwise; the sensitivity run including `conflict=True` is in the log.

## Numbers table per property x kind x split: RMSE, MAE, 95 % CI from your bootstrap, frac_within_tolerance, over-prediction rate (flash)

### viscosity_dyn_mPas

log-space score = RMSE of ln(eta / mPa.s); tolerance 5 % relative. `gn_mass`/`gn_mole` = Grunberg-Nissan k_ij=0 (mass / mole weighting); `refutas` = VBN blending on kinematic viscosity, needs pure density (records with any nu_i < 0.2 cSt skipped).

records in table: 18873, evaluated: 16036 (85.0 %)

| baseline | split | kind | n | rmse 95% CI (ln(mPa.s)) | MAE | within tol (5%) |
|---|---|---|---|---|---|---|
| gn_mass | S1 | pure | 1725 | 0.0596 [0.0426, 0.0789] | 0.029 | 85.2% |
| gn_mass | S1 | binary | 9878 | 0.0951 [0.0808, 0.1148] | 0.071 | 42.5% |
| gn_mass | S1 | multi | 261 | 0.0767 [0.0608, 0.1370] | 0.065 | 37.9% |
| gn_mass | S2 | pure | 246 | 0.0396 [0.0269, 0.0535] | 0.024 | 86.2% |
| gn_mass | S2 | binary | 3527 | 0.1260 [0.1048, 0.1496] | 0.094 | 36.6% |
| gn_mass | S2 | multi | 57 | 0.1380 [0.1105, 0.1929] | 0.120 | 14.0% |
| gn_mass | S3 | pure | 248 | 0.0474 [0.0233, 0.0732] | 0.027 | 85.1% |
| gn_mass | S3 | binary | 21 | 6.7678 [6.7678, 6.7678] | 6.768 | 0.0% |
| gn_mole | S1 | pure | 1725 | 0.0596 [0.0426, 0.0789] | 0.029 | 85.2% |
| gn_mole | S1 | binary | 9878 | 0.0916 [0.0742, 0.1155] | 0.064 | 50.8% |
| gn_mole | S1 | multi | 261 | 0.0773 [0.0634, 0.1357] | 0.066 | 34.1% |
| gn_mole | S2 | pure | 246 | 0.0396 [0.0269, 0.0535] | 0.024 | 86.2% |
| gn_mole | S2 | binary | 3527 | 0.1025 [0.0843, 0.1245] | 0.074 | 47.7% |
| gn_mole | S2 | multi | 57 | 0.0895 [0.0657, 0.0974] | 0.080 | 19.3% |
| gn_mole | S3 | pure | 248 | 0.0474 [0.0233, 0.0732] | 0.027 | 85.1% |
| gn_mole | S3 | binary | 21 | 6.7447 [6.7447, 6.7447] | 6.745 | 0.0% |
| refutas | S1 | pure | 1145 | 0.0588 [0.0355, 0.0852] | 0.030 | 85.5% |
| refutas | S1 | binary | 9458 | 0.0584 [0.0536, 0.0633] | 0.046 | 62.1% |
| refutas | S1 | multi | 261 | 0.0613 [0.0494, 0.1055] | 0.050 | 52.5% |
| refutas | S2 | pure | 157 | 0.0403 [0.0225, 0.0623] | 0.026 | 87.3% |
| refutas | S2 | binary | 3236 | 0.0761 [0.0652, 0.0864] | 0.061 | 51.3% |
| refutas | S2 | multi | 49 | 0.0987 [0.0738, 0.1389] | 0.086 | 22.4% |
| refutas | S3 | pure | 119 | 0.0462 [0.0271, 0.0693] | 0.033 | 79.0% |
| refutas | S3 | binary | 21 | 6.7775 [6.7775, 6.7775] | 6.777 | 0.0% |
| null | S1 | pure | 1725 | 0.9500 [0.7960, 1.1446] | 0.712 | 5.5% |
| null | S1 | binary | 9878 | 0.6096 [0.5634, 0.6563] | 0.492 | 5.7% |
| null | S1 | multi | 261 | 0.4339 [0.3655, 0.7262] | 0.386 | 2.3% |
| null | S2 | pure | 246 | 0.7130 [0.5380, 0.9361] | 0.570 | 4.1% |
| null | S2 | binary | 3527 | 0.6174 [0.5578, 0.6861] | 0.520 | 4.0% |
| null | S2 | multi | 57 | 0.5159 [0.4037, 0.7244] | 0.438 | 5.3% |
| null | S3 | pure | 248 | 1.0899 [0.8111, 1.3227] | 0.923 | 2.0% |
| null | S3 | binary | 21 | 7.3906 [7.3906, 7.3906] | 7.389 | 0.0% |

† pure rows: for viscosity the value is the pure-compound fit evaluated at its own data (circular, included only to show the fit error); for density it is the identity.

### viscosity_kin_cSt

log-space score = RMSE of ln(nu / cSt); `refutas` predicts nu directly; `gn_mass`/`gn_mole` are converted with the ideal-volume mixture density.

records in table: 95, evaluated: 50 (52.6 %)

| baseline | split | kind | n | rmse 95% CI (ln(cSt)) | MAE | within tol (5%) |
|---|---|---|---|---|---|---|
| gn_mass | S1 | pure | 48 | 0.1165 [0.0346, 0.1557] | 0.086 | 54.2% |
| gn_mass | S2 | pure | 2 | 0.0182 [0.0182, 0.0182] | 0.018 | 100.0% |
| gn_mole | S1 | pure | 48 | 0.1165 [0.0346, 0.1557] | 0.086 | 54.2% |
| gn_mole | S2 | pure | 2 | 0.0182 [0.0182, 0.0182] | 0.018 | 100.0% |
| refutas | S1 | pure | 48 | 0.1165 [0.0346, 0.1557] | 0.086 | 54.2% |
| refutas | S2 | pure | 2 | 0.0182 [0.0182, 0.0182] | 0.018 | 100.0% |
| null | S1 | pure | 48 | 1.5416 [0.6900, 1.9763] | 1.192 | 4.2% |
| null | S2 | pure | 2 | 0.1618 [0.1618, 0.1618] | 0.159 | 0.0% |

### density_kg_m3

score = RMSE in kg/m3; `dens_idvol` = 1 / sum(w_i/rho_i). Pure rows are identity by construction (v_i = 1) - marked with † - and are not evidence about the physics.

records in table: 15200, evaluated: 12840 (84.5 %)

| baseline | split | kind | n | rmse 95% CI (kg/m3) | MAE | within tol (5%) |
|---|---|---|---|---|---|---|
| dens_idvol | S1 | pure | 1447 | 0.0000 [0.0000, 0.0000] | 0.000 | 100.0% |
| dens_idvol | S1 | binary | 7811 | 1.8269 [1.4639, 2.2426] | 1.089 | 100.0% |
| dens_idvol | S1 | multi | 349 | 7.7790 [1.1664, 9.2978] | 4.580 | 99.7% |
| dens_idvol | S2 | pure | 166 | 0.0000 [0.0000, 0.0000] | 0.000 | 100.0% |
| dens_idvol | S2 | binary | 2852 | 2.5686 [1.4690, 3.6416] | 1.498 | 100.0% |
| dens_idvol | S2 | multi | 76 | 1.4284 [0.7268, 2.7120] | 1.024 | 100.0% |
| dens_idvol | S3 | pure | 139 | 0.0000 [0.0000, 0.0000] | 0.000 | 100.0% |
| null | S1 | pure | 1447 | 98.6348 [84.9926, 113.3058] | 80.517 | 25.8% |
| null | S1 | binary | 7811 | 68.2094 [60.3053, 75.8399] | 52.532 | 48.0% |
| null | S1 | multi | 349 | 58.4606 [33.1174, 64.9970] | 51.596 | 31.2% |
| null | S2 | pure | 166 | 77.1544 [57.3378, 97.9074] | 67.747 | 22.9% |
| null | S2 | binary | 2852 | 71.3725 [63.6168, 80.3286] | 60.055 | 32.1% |
| null | S2 | multi | 76 | 64.1851 [40.0395, 69.8402] | 55.579 | 34.2% |
| null | S3 | pure | 139 | 129.2577 [81.1302, 168.2100] | 103.915 | 20.1% |

† pure rows: for viscosity the value is the pure-compound fit evaluated at its own data (circular, included only to show the fit error); for density it is the identity.

### flash_point_C

score = RMSE in degC; `chat_nominal` = Le Chatelier + Raoult (gamma=1) + Antoine with LFL25; `chat_lflT` adds Zabetakis LFL(T). over-prediction = pred - measured > 5 degC (the unsafe direction: the model claims a safer fuel than it is).

records in table: 915, evaluated: 148 (16.2 %)

| baseline | split | kind | n | rmse 95% CI (degC) | MAE | within tol | over-pred >5 C | mean +err (C) | |err|<=5 C |
|---|---|---|---|---|---|---|---|---|---|
| chat_nominal | S1 | pure | 35 | 22.1733 [11.0385, 32.8628] | 13.270 | 34.3% | 40.0% [23,57] | 15.60 | 37.1% |
| chat_nominal | S1 | binary | 31 | 2.0192 [1.5370, 2.7594] | 1.770 | 9.7% | 0.0% [0,0] | 1.77 | 96.8% |
| chat_nominal | S1 | multi | 3 | 3.1663 [3.1663, 3.1663] | 3.136 | 33.3% | 0.0% [0,0] | 0.00 | 100.0% |
| chat_nominal | S2 | pure | 3 | 19.7820 [4.1466, 25.3491] | 17.391 | 0.0% | 66.7% [0,100] | 24.01 | 33.3% |
| chat_nominal | S2 | binary | 22 | 1.6420 [1.6228, 2.0040] | 1.354 | 9.1% | 0.0% [0,0] | 0.00 | 100.0% |
| chat_nominal | S3 | pure | 5 | 13.9771 [12.5971, 15.2325] | 13.893 | 0.0% | 40.0% [0,80] | 14.19 | 0.0% |
| chat_lflT | S1 | pure | 35 | 21.0721 [10.0253, 31.4416] | 12.235 | 34.3% | 37.1% [20,54] | 12.99 | 40.0% |
| chat_lflT | S1 | binary | 31 | 1.9484 [1.4991, 2.7915] | 1.673 | 12.9% | 0.0% [0,0] | 1.95 | 96.8% |
| chat_lflT | S1 | multi | 3 | 3.6512 [3.6512, 3.6512] | 3.626 | 33.3% | 0.0% [0,0] | 0.00 | 100.0% |
| chat_lflT | S2 | pure | 3 | 15.9885 [3.7192, 20.5642] | 14.151 | 0.0% | 66.7% [0,100] | 19.37 | 33.3% |
| chat_lflT | S2 | binary | 22 | 1.4890 [1.4245, 2.4825] | 1.186 | 22.7% | 0.0% [0,0] | 0.14 | 100.0% |
| chat_lflT | S3 | pure | 5 | 14.2301 [12.2747, 16.1308] | 14.058 | 0.0% | 40.0% [0,80] | 12.30 | 0.0% |
| null | S1 | pure | 35 | 104.8438 [86.1229, 121.6430] | 86.883 | 5.7% | 20.0% [9,34] | 35.33 | 8.6% |
| null | S1 | binary | 31 | 57.2568 [38.8270, 60.3902] | 55.178 | 0.0% | 90.3% [57,98] | 58.25 | 0.0% |
| null | S1 | multi | 3 | 17.4031 [17.4031, 17.4031] | 16.283 | 0.0% | 0.0% [0,0] | 0.00 | 0.0% |
| null | S2 | pure | 3 | 153.0272 [56.5750, 183.1000] | 140.925 | 0.0% | 33.3% [0,100] | 56.57 | 0.0% |
| null | S2 | binary | 22 | 45.7016 [9.7500, 46.7287] | 44.705 | 0.0% | 95.5% [0,100] | 46.37 | 0.0% |
| null | S3 | pure | 5 | 125.3564 [90.0504, 152.7083] | 119.490 | 0.0% | 20.0% [0,60] | 105.50 | 0.0% |

### vapor_pressure_kPa

score = RMSE of log10(P/kPa); this is a sanity check that `compounds.csv` Antoine reproduces the pure vapour-pressure measurements it was (partly) refitted from - it is **not** an independent validation.

records in table: 6437, evaluated: 4106 (63.8 %)

| baseline | split | kind | n | rmse 95% CI (log10(kPa)) | MAE | within tol |
|---|---|---|---|---|---|---|
| psat_antoine | S1 | pure | 3348 | 0.0140 [0.0092, 0.0195] | 0.006 | 96.2% |
| psat_antoine | S2 | pure | 418 | 0.0067 [0.0048, 0.0083] | 0.004 | 99.0% |
| psat_antoine | S3 | pure | 248 | 0.0088 [0.0035, 0.0125] | 0.005 | 96.4% |
| null | S1 | pure | 3348 | 2.2222 [1.9829, 2.4803] | 1.635 | 1.2% |
| null | S2 | pure | 418 | 2.5791 [1.9004, 3.2931] | 1.849 | 1.0% |
| null | S3 | pure | 248 | 2.7111 [2.0006, 3.5955] | 2.137 | 0.8% |

### Mixture-only summary (binary + multi, conflict=False)

Physics baselines only, n-weighted across kinds within each split; this is the honest headline (pure rows excluded when they are identity/self-fit). The comparable null level, on exactly the same rows, is in the paired bootstrap table below. The S3 viscosity cells are dominated by impossible source labels (system `a175950749e8`, see *What failed*); use the S1/S2 columns for generalisation.

| property | baseline | split | n(mixtures) | n-weighted score | within 5% rel |
|---|---|---|---|---|---|
| viscosity_dyn_mPas | gn_mass | S1 | 10139 | 0.0947 | 42.4% |
| viscosity_dyn_mPas | gn_mass | S2 | 3584 | 0.1262 | 36.2% |
| viscosity_dyn_mPas | gn_mass | S3 | 21 | 6.7678 | 0.0% |
| viscosity_dyn_mPas | gn_mole | S1 | 10139 | 0.0912 | 50.4% |
| viscosity_dyn_mPas | gn_mole | S2 | 3584 | 0.1023 | 47.3% |
| viscosity_dyn_mPas | gn_mole | S3 | 21 | 6.7447 | 0.0% |
| viscosity_dyn_mPas | refutas | S1 | 9719 | 0.0585 | 61.9% |
| viscosity_dyn_mPas | refutas | S2 | 3285 | 0.0765 | 50.8% |
| viscosity_dyn_mPas | refutas | S3 | 21 | 6.7775 | 0.0% |
| density_kg_m3 | dens_idvol | S1 | 8160 | 2.0814 | 100.0% |
| density_kg_m3 | dens_idvol | S2 | 2928 | 2.5390 | 100.0% |
| flash_point_C | chat_nominal | S1 | 34 | 2.1204 | 11.8% |
| flash_point_C | chat_nominal | S2 | 22 | 1.6420 | 9.1% |
| flash_point_C | chat_lflT | S1 | 34 | 2.0987 | 14.7% |
| flash_point_C | chat_lflT | S2 | 22 | 1.4890 | 22.7% |

Success criterion (reports/NOISE_FLOOR.md): a model only counts as an improvement over physics if it beats the physics baseline by at least the noise floor with a non-overlapping 95 % CI - viscosity 0.016 log-units median spread (p90 0.12), density 0.73 kg/m3 (p90 5.44), flash 6 degC (p90 29.4). The same bar applies to physics vs the null model below.

### Does physics beat the null model?

Per split, n-weighted across kinds and across the physics cells listed for that property (each cell pairs the physics baseline with the null on the same rows; pure rows included, so density/viscosity pure cells are identity/self-fit and not evidence). The paired, mixture-only, same-rows comparison with bootstrap CIs is the second table.

| property | baseline | split | physics score | null score | delta | verdict |
|---|---|---|---|---|---|---|
| viscosity_dyn_mPas | gn_mass | S1 | 0.0896 | 0.6488 | -86.2% | physics better |
| viscosity_dyn_mPas | gn_mass | S2 | 0.1206 | 0.6225 | -80.6% | physics better |
| viscosity_dyn_mPas | gn_mass | S3 | 0.5720 | 1.6340 | -65.0% | physics better |
| viscosity_dyn_mPas | gn_mole | S1 | 0.0866 | 0.6488 | -86.7% | physics better |
| viscosity_dyn_mPas | gn_mole | S2 | 0.0982 | 0.6225 | -84.2% | physics better |
| viscosity_dyn_mPas | gn_mole | S3 | 0.5702 | 1.6340 | -65.1% | physics better |
| viscosity_dyn_mPas | refutas | S1 | 0.0585 | 0.6488 | -91.0% | physics better |
| viscosity_dyn_mPas | refutas | S2 | 0.0748 | 0.6225 | -88.0% | physics better |
| viscosity_dyn_mPas | refutas | S3 | 1.0559 | 1.6340 | -35.4% | physics better |
| viscosity_kin_cSt | gn_mass | S1 | 0.1165 | 1.5416 | -92.4% | physics better |
| viscosity_kin_cSt | gn_mass | S2 | 0.0182 | 0.1618 | -88.7% | physics better |
| viscosity_kin_cSt | gn_mole | S1 | 0.1165 | 1.5416 | -92.4% | physics better |
| viscosity_kin_cSt | gn_mole | S2 | 0.0182 | 0.1618 | -88.7% | physics better |
| viscosity_kin_cSt | refutas | S1 | 0.1165 | 1.5416 | -92.4% | physics better |
| viscosity_kin_cSt | refutas | S2 | 0.0182 | 0.1618 | -88.7% | physics better |
| density_kg_m3 | dens_idvol | S1 | 1.7679 | 72.4379 | -97.6% | physics better |
| density_kg_m3 | dens_idvol | S2 | 2.4028 | 71.5062 | -96.6% | physics better |
| density_kg_m3 | dens_idvol | S3 | 0.0000 | 129.2577 | -100.0% | physics better |
| flash_point_C | chat_nominal | S1 | 12.2922 | 79.6624 | -84.6% | physics better |
| flash_point_C | chat_nominal | S2 | 3.8188 | 58.5807 | -93.5% | physics better |
| flash_point_C | chat_nominal | S3 | 13.9771 | 125.3564 | -88.9% | physics better |
| flash_point_C | chat_lflT | S1 | 11.7229 | 79.6624 | -85.3% | physics better |
| flash_point_C | chat_lflT | S2 | 3.2289 | 58.5807 | -94.5% | physics better |
| flash_point_C | chat_lflT | S3 | 14.2301 | 125.3564 | -88.6% | physics better |
| vapor_pressure_kPa | psat_antoine | S1 | 0.0140 | 2.2222 | -99.4% | physics better |
| vapor_pressure_kPa | psat_antoine | S2 | 0.0067 | 2.5791 | -99.7% | physics better |
| vapor_pressure_kPa | psat_antoine | S3 | 0.0088 | 2.7111 | -99.7% | physics better |

Pooled binary+multi rows (conflict-free), same rows on both sides, RMSE with a 95 % cluster-bootstrap CI - 'separated' means the physics CI lies entirely below the null CI:

| property | baseline | split | n | physics RMSE [lo, hi] | null RMSE [lo, hi] | separated |
|---|---|---|---|---|---|---|
| viscosity_dyn_mPas | gn_mass | S1 | 10139 | 0.0947 [0.0809, 0.1137] | 0.6057 [0.5598, 0.6536] | yes |
| viscosity_dyn_mPas | gn_mass | S2 | 3584 | 0.1262 [0.1044, 0.1499] | 0.6160 [0.5540, 0.6800] | yes |
| viscosity_dyn_mPas | gn_mass | S3 | 21 | 6.7678 [6.7678, 6.7678] | 7.3906 [7.3906, 7.3906] | yes |
| viscosity_dyn_mPas | gn_mole | S1 | 10139 | 0.0912 [0.0747, 0.1152] | 0.6057 [0.5598, 0.6536] | yes |
| viscosity_dyn_mPas | gn_mole | S2 | 3584 | 0.1023 [0.0835, 0.1248] | 0.6160 [0.5540, 0.6800] | yes |
| viscosity_dyn_mPas | gn_mole | S3 | 21 | 6.7447 [6.7447, 6.7447] | 7.3906 [7.3906, 7.3906] | yes |
| viscosity_dyn_mPas | refutas | S1 | 9719 | 0.0585 [0.0535, 0.0638] | 0.6031 [0.5559, 0.6520] | yes |
| viscosity_dyn_mPas | refutas | S2 | 3285 | 0.0765 [0.0661, 0.0862] | 0.6212 [0.5550, 0.6895] | yes |
| viscosity_dyn_mPas | refutas | S3 | 21 | 6.7775 [6.7775, 6.7775] | 7.3906 [7.3906, 7.3906] | yes |
| density_kg_m3 | dens_idvol | S1 | 8160 | 2.4047 [1.5507, 3.3703] | 67.8211 [60.2972, 75.3462] | yes |
| density_kg_m3 | dens_idvol | S2 | 2928 | 2.5455 [1.4714, 3.6224] | 71.1951 [63.2573, 79.8857] | yes |
| flash_point_C | chat_nominal | S1 | 34 | 2.1452 [1.7484, 2.8726] | 54.9163 [27.0944, 59.4980] | yes |
| flash_point_C | chat_nominal | S2 | 22 | 1.6420 [1.6228, 2.0040] | 45.7016 [9.7500, 46.7287] | yes |
| flash_point_C | chat_lflT | S1 | 34 | 2.1535 [1.6282, 3.0272] | 54.9163 [27.0944, 59.4980] | yes |
| flash_point_C | chat_lflT | S2 | 22 | 1.4890 [1.4245, 2.4825] | 45.7016 [9.7500, 46.7287] | yes |

## real vs synthetic

Every number in this report is a **real-data** number computed from `data/processed/model_table.csv` (measured densities, viscosities, flash points and vapour pressures, de-duplicated medians). No synthetic label is used anywhere in the evaluation. The only synthetic values in this work package are the hand-made fixtures inside `tests/test_baselines_eval.py` (invented pure viscosities and Antoine constants, used to prove the code implements the formulas); they are clearly labelled in the test file and are **not** reported as results.

## What failed / unevaluable (counts + reasons)

| property | records | evaluable (any baseline) | share |
|---|---|---|---|
| viscosity_dyn_mPas | 18873 | 16036 | 85.0% |
| viscosity_kin_cSt | 95 | 50 | 52.6% |
| density_kg_m3 | 15200 | 12840 | 84.5% |
| flash_point_C | 915 | 148 | 16.2% |
| vapor_pressure_kPa | 6437 | 4106 | 63.8% |

Per-baseline skip reasons (a record can be skipped by one baseline and predicted by another):

| property | baseline | reason | records |
|---|---|---|---|
| vapor_pressure_kPa | psat_antoine | not_pure | 2022 |
| viscosity_dyn_mPas | refutas | no_pure_rho_at_T | 1759 |
| flash_point_C | chat_lflT | no_antoine | 766 |
| flash_point_C | chat_nominal | no_antoine | 766 |
| density_kg_m3 | dens_idvol | no_pure_rho_at_T | 443 |
| viscosity_dyn_mPas | gn_mass | no_pure_eta_at_T | 428 |
| viscosity_dyn_mPas | gn_mole | no_pure_eta_at_T | 428 |
| viscosity_dyn_mPas | refutas | no_pure_eta_at_T | 428 |
| vapor_pressure_kPa | psat_antoine | no_antoine | 309 |
| viscosity_kin_cSt | gn_mass | no_pure_eta_at_T | 45 |
| viscosity_kin_cSt | gn_mole | no_pure_eta_at_T | 45 |
| viscosity_kin_cSt | refutas | no_pure_eta_at_T | 45 |
| density_kg_m3 | dens_idvol | no_T | 3 |
| viscosity_dyn_mPas | gn_mass | no_T | 3 |
| viscosity_dyn_mPas | gn_mole | no_T | 3 |
| viscosity_dyn_mPas | refutas | no_T | 3 |
| viscosity_dyn_mPas | refutas | nu_lt_0p2 | 1 |

Evaluable records per property x kind (a record counts as evaluable if at least one baseline produced a number):

| property | kind | records | evaluable |
|---|---|---|---|
| viscosity_dyn_mPas | pure | 2270 | 2269 |
| viscosity_dyn_mPas | binary | 13754 | 13449 |
| viscosity_dyn_mPas | multi | 346 | 318 |
| viscosity_kin_cSt | pure | 67 | 50 |
| viscosity_kin_cSt | binary | 28 | 0 |
| density_kg_m3 | pure | 1752 | 1752 |
| density_kg_m3 | binary | 10982 | 10663 |
| density_kg_m3 | multi | 454 | 425 |
| flash_point_C | pure | 789 | 92 |
| flash_point_C | binary | 84 | 53 |
| flash_point_C | multi | 3 | 3 |
| vapor_pressure_kPa | pure | 4415 | 4106 |
| vapor_pressure_kPa | binary | 1685 | 0 |
| vapor_pressure_kPa | multi | 1 | 0 |

**Real fuels**: 0 real-fuel records for these five properties are in the table and 0 of them have a molecule list, so these baselines cannot score them at all. They are held for a distribution-input model (class x carbon-number, GCxGC) and for inverse-design validation - listed as future work, never silently dropped (every one is in `baseline_preds.csv` with NaN predictions; the counted reason `no_molecule_list` is in the log's skip-reason table).

**S3** (792 rows) was evaluated and is reported where a baseline reached it; S3 holds out whole molecules, but Mode A physics does not need those molecules to be seen, so only Antoine/LFL/pure-curve availability limits coverage. Cells with zero coverage are simply absent from the tables above.

### Worst records (audit trail, not excluded)

Largest errors of the best baseline per property (mixture rows, conflict-free). They are kept in every number above - nothing here is filtered out.

| property | baseline | system_id | split | kind | T_K | measured | predicted | |err| |
|---|---|---|---|---|---|---|---|---|
| viscosity_dyn_mPas | refutas | a175950749e8 | S3 | binary | 303.1 | 1546.0000 | 1.6571 | 6.838 |
| viscosity_dyn_mPas | refutas | a175950749e8 | S3 | binary | 323.1 | 1115.0000 | 1.1986 | 6.835 |
| viscosity_kin_cSt | refutas | 34c6d0ac127e | S1 | pure | 345.1 | 0.1903 | 0.2443 | 0.250 |
| viscosity_kin_cSt | refutas | 34c6d0ac127e | S1 | pure | 340.1 | 0.1971 | 0.2505 | 0.240 |
| density_kg_m3 | dens_idvol | 42c008af920c | S1 | multi | 293.1 | 728.6000 | 765.4263 | 36.826 |
| density_kg_m3 | dens_idvol | 42c008af920c | S1 | multi | 303.1 | 720.7000 | 756.0343 | 35.334 |
| flash_point_C | chat_lflT | 82f488beee5b | S1 | pure | - | 120.0000 | 209.5576 | 89.558 |
| flash_point_C | chat_lflT | d54904a86050 | S1 | pure | - | 160.3000 | 216.5964 | 56.296 |
| vapor_pressure_kPa | psat_antoine | f6c175883f73 | S1 | pure | 302.1 | 12.0000 | 7.4240 | 0.480 |
| vapor_pressure_kPa | psat_antoine | f9e1cfd55b52 | S1 | pure | 300.1 | 14.0000 | 22.0658 | 0.455 |

**Impossible-looking labels.** A physics baseline cannot be off by 3 orders of magnitude on a mixture of well-characterised components, so these rows are almost certainly unit/assignment errors in the source records (e.g. mPa.s vs Pa.s families) and they dominate the S3 viscosity cells:

- `viscosity_dyn_mPas`: 21 of 13767 mixture records are more than 50x away from the ideal/log-linear prediction. Example: system `a175950749e8` (binary, 313.1 K) measured 1304.000, predicted 1.412.
- `density_kg_m3`: 0 of 11088 mixture records are more than 50x away from the ideal/log-linear prediction.

## Decisions made

1. **Component curves from the pure rows of the same table**, `ln eta = a + b ln T` and `rho = a + b T`, evaluated only inside (2 K pad) each compound's own measured range. This replaces the earlier nearest-measurement rule (`baselines.py`, 7.5 K tolerance) with interpolation - coverage went up and no molecule is extrapolated.
2. **Mode A is enforced structurally**: no component curve -> no prediction, with the reason logged. We did not impute missing pure properties; imputation would make the baselines non-reproducible.
3. **Refutas needs kinematic viscosity**, so it also needs pure density; records where any component lacks a density curve are reported as `no_pure_rho_at_T` rather than silently switching to the log-linear baseline. Records with any nu_i < 0.2 cSt are excluded (`nu_lt_0p2`) because Refutas is undefined there.
4. **Pure rows of density (and of viscosity) are identity/circular**, so they are marked † and excluded from the mixture headline summary. Flash and vapour-pressure pure rows are genuine predictions and are kept.
5. **Null model is split-aware** (S1 medians for S2/S3, system-disjoint folds for S1) so the physics comparison is not flattered by an in-sample null. The null is binned by `kind` and a 20 K temperature bin; for flash (no condition T in the table) that degenerates to a per-kind constant, which is the right sanity reference - it is the 'always predict the median flash point' model.
6. **`conflict=False` is the headline filter**; the sensitivity run with conflicts included is in the log. `n_obs`-based `quality_weight` was not used inside the baselines (they are per-record, and weighting would hide the disagreement).
7. **No tuning of any kind**: baselines are formulas; the only fitted objects are the per-compound pure-property curves, which use pure rows only (never a blend row), so a blend record cannot leak into its own prediction.
8. Mole fractions for flash are recomputed from `w_mass` and the registry molar masses, as specified, not taken from the stored `x_mole` column.

## Next steps (what a learned residual model should target)

- Residual target per property, from the tables above (mixture cells, conflict-free):
  - `viscosity_dyn_mPas` / `gn_mass`: S1 0.0947 -> S2 0.1262 ln(mPa.s) (multiplicative correction in ln space; the S1->S2 gap is the molecule-generalisation headroom).
  - `viscosity_dyn_mPas` / `gn_mole`: S1 0.0912 -> S2 0.1023 ln(mPa.s) (multiplicative correction in ln space; the S1->S2 gap is the molecule-generalisation headroom).
  - `viscosity_dyn_mPas` / `refutas`: S1 0.0585 -> S2 0.0765 ln(mPa.s) (multiplicative correction in ln space; the S1->S2 gap is the molecule-generalisation headroom).
  - `density_kg_m3` / `dens_idvol`: S1 2.0814 -> S2 2.5390 kg/m3 (additive correction in the target unit; the S1->S2 gap is the molecule-generalisation headroom).
  - `flash_point_C` / `chat_nominal`: S1 2.1204 -> S2 1.6420 degC (additive correction in the target unit; the S1->S2 gap is the molecule-generalisation headroom).
  - `flash_point_C` / `chat_lflT`: S1 2.0987 -> S2 1.4890 degC (additive correction in the target unit; the S1->S2 gap is the molecule-generalisation headroom).
- Learn the **activity-coefficient-like** term for flash point: ideal Raoult (`chat_*`) cannot bend below the lightest component's own flash point, and the binary alkane mixtures in this table do; that is the residual with the clearest physical meaning.
- Learn the **excess viscosity** correction; keep Refutas as the floor (it is the strongest baseline on binaries) and only allow the learned term to act for n_comp >= 3.
- Density: the ideal-volume rule should already be close; any residual is a volume-contraction signal, so fit it as a function of aromatic content and T rather than as a free-form correction.
- Flash point must be trained with an **asymmetric loss** on over-prediction; report `over-pred >5 C` as a first-class metric with its bootstrap CI, exactly as tabulated here.
- Unevaluable records (no Antoine / no LFL / no pure curve) are the coverage ceiling for any residual model: fixing them is a data task (interpolate Antoine vs carbon number, extend LFL), not a modelling task.
- Real fuels (no molecule list) stay out of supervised training; use them as a distribution-input validation set and for inverse-design sanity checks.

