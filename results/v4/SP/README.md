# Smoke point / sooting evaluation (v4 Task 14)

All numbers below are read from `smoke_point_eval.csv` (regenerate: `scripts/v4/build_sooting_inputs.py` then `scripts/v4/eval_smoke_point.py`).

## Results (mm; lower is better)

| slice | model | N | MAE | RMSE | bias | MAE minus const MAE, 95% bootstrap CI |
|---|---|---|---|---|---|---|
| pure_ysi_route | const_geomean | 23 | 6.98 | 13.09 | -3.28 |  |
| pure_ysi_route | MW/SP=c0+c1*YSI | 23 | 5.07 | 10.43 | -3.69 | [-2.95, -0.88] |
| pure_ysi_route_no_isooctane | MW/SP=c0+c1*YSI | 22 | 4.14 | 9.19 | -2.71 |  |
| ajftd100_bulk | const_mean | 100 | 1.14 | 1.39 | +0.00 |  |
| ajftd100_bulk | arom+diarom_fixed | 100 | 0.90 | 1.14 | +0.01 | [-0.36, -0.12] |
| ajftd100_bulk | selected_in_fold | 100 | 0.90 | 1.14 | +0.01 | [-0.36, -0.12] |
| ajftd_broad_bulk | const_mean | 849 | 1.65 | 3.10 | -0.00 |  |
| ajftd_broad_bulk | arom+diarom_fixed | 849 | 1.68 | 2.30 | -0.00 | [-0.09, +0.15] |
| ajftd_broad_bulk | selected_in_fold | 849 | 1.68 | 2.31 | -0.00 | [-0.09, +0.15] |

## Censoring

* Pure compounds: 5 compounds have true SP > 50 mm (D1322 lamp limit); never fit as exact targets. Fraction the YSI model predicts >= 50 mm: 0.00 (constant baseline 0.00). Predictions: {'n-hexane': 14.8, 'n-heptane': 16.5, 'cyclohexane': 13.2, 'n-decane': 20.2, 'n-dodecane': 22.1}
* AJFTD: 0 censored rows in the 100-fuel file.

## Honest reading

* Bulk route (AJFTD 100 fuels, GroupKFold by fuel_id, feature set selected inside training folds): MAE 0.90 mm vs constant-mean 1.14 mm (y range 19.3-25.0 mm, sd 1.38). AJFTD smoke point resolution is ~0.5-1 mm and fuel_ids from one producer/batch are not grouped, so the CV is if anything optimistic. No fuel is below 19 mm or above 25 mm, so the 18/25 mm constraint boundaries are not probed.
* Broad AJFTD sensitivity (N=849, D1322, 5<SP<50 filter): selected MAE 1.68 vs const 1.65.
* YSI route, pure compounds (leave-one-out, N=23): MAE 5.07 vs geometric-mean baseline 6.98.
* The linear MW/SP-vs-YSI form fails on saturated hydrocarbons: every censored (SP > 50 mm) alkane/cycloalkane is predicted far below 50 mm (see Censoring), so the YSI route is not reliable for low-sooting components. On the broad AJFTD set the bulk model has no skill over the constant mean (CI spans 0).
* YSI route on real fuels: NOT EVALUATED (no fuel with both speciation and measured smoke point). It must not be cited as evidence for fuel smoke point.

## Provenance and caveats

* YSI: raw_file, verified equal to `UnifiedYSI` in Harvard Dataverse YSI Database Volume 2 (benzene 100.3, toluene 170.9, n-heptane 36.0).
* Pure smoke points (28 rows, 23 uncensored): `cited_typed` from extras/build_tsi_table.py (Mensch 2010, Li & Sunderland 2012, Boehm 2022 SI); SI numbers are not in this repo and could not be checked. Known discrepancy: iso-octane typed 40.0 mm, Boehm 2022 text quotes 42.8 mm.
* EXCLUDED as circular: 175 rows of extras `tsi_table.csv` (derived=True), smoke point computed from YSI.
* AJFTD smoke point/aromatics: raw_file (match `ajftd_fuels.csv`/`ajftd_composition.csv`); smoke point is not derived from aromatics.

## Claims-register ids that would update (not edited here)

Any claim listing smoke point as an unmodeled property (D1655 smoke point now has a bulk-composition baseline and a constraint entry); no pure/fuel sooting claim is supported beyond the rows above.
