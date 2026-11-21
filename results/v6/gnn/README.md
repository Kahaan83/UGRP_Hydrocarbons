# v6 Track B: pure-compound D-MPNN

Pool: `pool.csv` (1,208 compounds, 10 targets; `scripts/v6/gnn/build_pool.py`). Sources: data/unified pure HC rows, AI4fuel
(HC and oxygenates), v5 ThermoML fusion, v4 YSI table. Oxygenates carry only ron/mon labels. Tb and LHV are not trained
(no Tb source; AI4fuel heat-of-combustion units unverified). Density is ambient (unified: 293 K +-5; AI4fuel: no T given).
Viscosity is log10 mPa.s at 293 K.

Model: `src/fuelblend/ml/dmpnn.py`, 5-model ensemble, early stopping on an inner group-held-out validation split,
2 configurations chosen per outer fold (`chosen_cfg.json`). CV: GroupKFold(5) on generic Murcko scaffold (acyclics bucketed
by carbon count), out-of-fold scoring (`oof.csv`, `cv_summary.csv`). No locked holdout and no conformal margins yet.

| target | n | MAE GNN | MAE constant | gain 90% CI (const - GNN) | existing estimator MAE |
|---|---|---|---|---|---|
| cetane | 259 | 9.04 | 18.57 | [7.84, 12.32] | |
| dcn | 89 | 5.56 | 18.37 | [9.77, 17.85] | |
| density kg/m3 | 558 | 21.90 | 107.31 | [71.27, 102.65] | |
| dHfus kJ/mol | 81 | 12.80 | 19.14 | [2.78, 10.74] | |
| flash C | 837 | 19.11 | 55.22 | [26.77, 51.24] | |
| log10 visc 293 K | 180 | 0.225 | 0.457 | [0.143, 0.347] | |
| mon | 363 | 8.24 | 14.53 | [5.12, 7.63] | |
| ron | 378 | 8.99 | 16.58 | [5.05, 9.50] | |
| Tm K | 547 | 30.42 | 71.82 | [31.93, 56.43] | 28.36 (v5 fusion_estimate) |
| YSI | 224 | 69.40 | 290.69 | [162.49, 330.10] | 106.46 (v4 ridge GC) |

Readings
- All heads beat the constant, but the constant is a pooled-median baseline over mixed classes, so it is a weak bar. Absolute errors
  are large for property use (flash MAE 19 C, density MAE 22 kg/m3, one-sided flash over-prediction mean 9.4 C): no head is
  accurate enough to replace physics or measured values for those properties.
- YSI: GNN MAE 69.4 vs ridge 106.5 on the same folds. Positive, but the v5 ridge failed its gate, and no holdout was used here.
- Tm: status: no improvement over the existing estimator (30.4 vs 28.4 K). Not used to fill Track A isomer constants.
- dHfus has no existing comparator and n = 81; the CI is wide.
- Not done: locked holdout, conformal margins, feeding missing Tm/dHfus into Track A tiles (Tm head does not beat the existing estimator;
  Track A result "no skill shown" is unchanged).
