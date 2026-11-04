# RESIDUAL_MODELS.md - WP5 Mixture Residual Models (Phase A / Post-R1)

## What ran - exact commands + seeds

```bash
python scripts/train_mixture_residual.py
```

- **Script**: `scripts/train_mixture_residual.py`
- **Output files**: `data/runs/RESIDUAL_EVAL_LOG.md`
- **Seeds**: model seeds `[11, 23, 37, 41, 59]`, split seed `20261030`
- **Model Table**: `data/processed/model_table.csv` filtered strictly to Rule R1 compliant datasets (42,662 rows, 42,328 non-locked).
- **Architecture**: 5-seed ensemble of `HistGradientBoostingRegressor` (max_iter=150, learning_rate=0.05, early stopping) trained on permutation-invariant mixture features (`src/fuelblend/features/mixture_features.py`).
- **Target definitions**:
  - Viscosity (`viscosity_dyn_mPas`): Residual in log space: $\Delta \ln \eta = \ln(\eta_{\text{measured}}) - \ln(\text{gn\_mass})$.
  - Density (`density_kg_m3`): Residual in linear space: $\Delta \rho = \rho_{\text{measured}} - \rho_{\text{idvol}}$.

## Numbers Table: Physics Baseline vs. Residual Model

Evaluated on S1 (interpolation training pool) and held-out S2 (molecule holdout test set) partitions:

| Property | Metric Space | Split | N | Physics Baseline | Residual Model | Error Reduction (%) |
|---|---|---|---|---|---|---|
| **viscosity_dyn_mPas** | $\ln(\text{mPa}\cdot\text{s})$ RMSE | **S1** (Train) | 11,864 | 0.0904 | 0.0336 | 62.8 % |
| | $\ln(\text{mPa}\cdot\text{s})$ RMSE | **S2** (Test) | 3,830 | 0.1225 | **0.0781** | **36.2 %** |
| **density_kg_m3** | Linear RMSE ($\text{kg/m}^3$) | **S1** (Train) | 9,607 | 2.216 | 0.969 | 56.3 % |
| | Linear RMSE ($\text{kg/m}^3$) | **S2** (Test) | 3,094 | 2.476 | **2.459** | **0.7 %** |

### Observations:
1. **Viscosity**: Substantial generalisation gain on novel molecule mixtures (S2). Physics Grunberg-Nissan log-linear mixing baseline achieves ln-RMSE 0.1225; the learned residual model cuts errors by 36.2% to 0.0781.
2. **Density**: The ideal volume baseline ($\rho_{\text{idvol}}$) is already near experimental precision on hydrocarbon mixtures (~2.48 kg/m³ error across S2). The residual model fits S1 excess volume closely (RMSE 0.969 kg/m³) and matches/slightly edges baseline on S2 (2.459 kg/m³).

## Compliance with Rule R1
- Evaluated on verified L3 gravimetric/molar mixtures and pure hydrocarbons.
- Uncharacterized L1 and L0 blends were filtered to `data/processed/excluded_by_R1.csv` prior to split generation and model training.
