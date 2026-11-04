# PURE_MODELS.md - WP4 Pure-Compound Property Models

## What ran / exact commands + seeds

```bash
python scripts/train_pure_models.py
```

- **Script**: `scripts/train_pure_models.py`
- **Output files**: `data/runs/pure_preds.csv` (2,758 rows), `data/runs/pure_feature_cache.csv`, `data/runs/PURE_EVAL_LOG.md`
- **Seeds**: model seeds `[11, 23, 37, 41, 59]`, split seed `20261030`
- **Features** (`src/fuelblend/features/pure_features.py`): 22 tabular features derived from SMILES via RDKit (`mol_wt`, `heavy_atom_count`, `n_aromatic_rings`, `ring_count`, `formal_charge`, `hbond_donor`, `hbond_acceptor`, `rotatable_bonds`, `labute_asa`, `tpsa`, `n_carbon`, `log_T`, `T_K`, `inv_T`, one-hot class indicators `cls_*`).
- **Model** (`src/fuelblend/models/pure_gbm.py`): 5-seed ensemble of `HistGradientBoostingRegressor` with early stopping and sample weights (`quality_weight = 1 / sqrt(n_obs)`).

## Numbers table per property x split

Evaluated on two held-out partitions:
1. **Heldback**: 20 % random molecule holdout from S1∪S2 (GroupShuffleSplit by `inchikey`, seed `20261030`).
2. **S3**: Strict all-novel-molecule extrapolation partition from `splits/s1s2s3.json`.

| Property | Target Space | Split | N | MAE | RMSE | 95 % Bootstrap CI | Top Features (Permutation Importance) |
|---|---|---|---|---|---|---|---|
| `viscosity_dyn_mPas` | `ln(eta)` | Heldback | 425 | 0.1796 | 0.2498 | [0.1804, 0.3254] | `heavy_atom_count`, `n_aromatic_rings`, `mol_wt` |
| | | S3 | 259 | 0.2289 | 0.6323 | [0.1792, 1.1333] | |
| `density_kg_m3` | Linear (kg/m³) | Heldback | 357 | 10.227 | 14.262 | [9.067, 20.675] | `n_carbon`, `cls_aromatic`, `tpsa` |
| | | S3 | 139 | 25.055 | 38.989 | [10.974, 57.473] | |
| `flash_point_C` | Linear (°C) | Heldback | 142 | 0.5615 | 1.5927 | [0.6801, 2.4739] | `cls_aromatic`, `tpsa`, `n_carbon` |
| | | S3 | 79 | 0.6458 | 1.1993 | [0.6888, 1.7340] | |
| `vapor_pressure_kPa` | `log10(P)` | Heldback | 1,070 | 0.2325 | 0.4096 | [0.2216, 0.6572] | `hbond_donor`, `formal_charge`, `n_aromatic_rings` |
| | | S3 | 287 | 0.1377 | 0.1990 | [0.1322, 0.2714] | |

### Safety Metrics for Flash Point (°C)
- **Heldback**: Bias = +0.224 °C; Over-prediction rate (> 5 °C) = 1.41 %; Mean positive error = 0.627 °C.
- **S3**: Bias = +0.205 °C; Over-prediction rate (> 5 °C) = 1.27 %; Mean positive error = 0.730 °C.

## Real vs Synthetic
All metrics presented in this report were computed strictly on actual experimental measurements from `data/processed/model_table.csv` and `data/unified/compounds.csv`. The synthetic test cases in `tests/test_pure_features.py` and `tests/test_pure_gbm.py` are unit-test fixtures used solely for code verification.

## What failed / skipped
- `freezing_point_C` and `pour_point_C`: Skipped from training because pure-compound training rows were < 200 (cold-flow pure count is sparse in open hydrocarbon datasets, per plan constraints).

## Decisions made
1. Target transformation: Viscosity trained in natural log space (`ln(mPa.s)`), vapor pressure in `log10(kPa)`. Density and flash point kept linear.
2. S3 Novel-Molecule Evaluation: Reported directly as the primary generalization test.
3. Ensemble: Mean of 5 seeds provides the point prediction; standard deviation across seeds provides an initial epistemic uncertainty proxy.

## Next steps
Feed these pure-compound predictions (`data/runs/pure_preds.csv`) into the Mode B mixture residual architecture (WP5) when component pure properties are not experimentally observed.
