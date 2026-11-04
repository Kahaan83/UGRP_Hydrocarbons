# T1 — Strict family-level holdout re-check of the pure flash claim

## What ran

```
python scripts/t1_leakage_audit.py
```

- Property `flash_point_C`, pure rows of `model_table.csv` in S1∪S2: **n=710 rows / 710 molecules / 122 scaffold families** (Murcko generus from `src/fuelblend/features/family.py`; acyclic chains bucketed by unsaturation, e.g. `open_u0` = all n-/iso-alkanes = 133 rows).
- 5-fold **leave-family-out** CV (seed 20261030): a whole scaffold family is held out of training in each fold; every row is scored exactly once across folds.
- Pipeline identical to WP4 (`collect/prepare/evaluate` reused from `scripts/train_pure_models.py`, HGB 5-seed ensembles [11, 23, 37, 41, 59]) — only the group unit changed (molecule → family).
- 2,000-resample cluster bootstrap by `system_id` for CIs.

## Numbers

| Split | n | MAE (°C) | RMSE (°C) | 95 % CI |
|---|---|---|---|---|
| WP4 heldback (molecule holdout, 20 % of mols) | 142 | 0.5615 | 1.5927 | RMSE [0.6801, 2.4739] |
| WP4 S3 (novel molecules) | 79 | 0.6458 | 1.1993 | RMSE [0.6888, 1.7340] |
| **T1 family holdout (5-fold, full pool)** | **710** | **0.7485** | **1.7332** | RMSE [1.3595, 2.1254], MAE [0.6441, 0.8650] |

Safety metrics at family strictness: bias +0.14 °C; over-prediction rate (>5 °C) 1.27 %; mean positive error 0.71 °C.

Biggest fold hits: fold 0 (multi-family holdout incl. `open_u1` + 20 scaffold families, MAE 1.04), fold 2 (MAE 1.12); smallest fold 3 (MAE 0.50, `open_u0` n-paraffins held out — easily interpolated).

## Verdict (plain language)

1. **The 0.56 °C claim does not survive strict family holdout as stated.** MAE rises 0.56 → 0.75 °C (+34 %): sibling molecules (e.g. close homologs sharing a scaffold) were in the WP4 training set, so part of the accuracy is molecular interpolation.
2. **But the effect is modest, not a collapse.** RMSE 1.59 → 1.73 °C stays **inside** the WP4 heldback 95 % CI [0.68, 2.47], and the family-holdout result (n=710, every row scored) is arguably a *fairer* number than the WP4 heldback snapshot (n=142 on a 20 % subset).
3. Recommended public claim: **pure flash ≈ 0.75 °C MAE / 1.73 °C RMSE under chemical-family holdout; 0.56–0.65 °C is the optimistic molecule-grouped bound.** Update PURE_MODELS.md sentiment accordingly; use the family number for any external claim.
4. Safety metrics degrade only slightly (bias +0.14 °C, ~1 % over-prediction at >5 °C) — the conformal flash boundary from WP6 remains valid in direction and magnitude.

## Real vs synthetic

All numbers computed on measured `model_table.csv` pure rows (real experimental data). The only synthetic data involved is `tests/test_t1_leakage_audit.py` fixtures (unit tests of the fold machinery itself; 3/3 pass).
