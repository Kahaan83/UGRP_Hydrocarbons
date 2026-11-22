# Task 14 reading: limiting-species freezing baseline
Predicted freezing point = max over species of the ideal-solubility temperature T_i (1/T_i = 1/Tm_i - R ln x_i / dHfus_i), designated constants as in v6 baseline (b), mole fractions from the same mapped GC x GC composition. No fitting. Rows: the 52 dev entities (same order as `results/v6/freezing/dev_oof_predictions.csv`, asserted). Gain = constant-baseline MAE minus model MAE per entity, 90 % percentile bootstrap (4000 draws, entities as groups). Comparator rows a, b, c1 and c2 are recomputed from the v6 file and match `dev_logo_summary.csv` (b: MAE 8.830, bias +5.513; c1: 16.717; c2: 6.434).

| model | MAE | bias | mean over-pred | frac over | p90 signed | gain vs constant [90 % CI] |
|---|---|---|---|---|---|---|
| a constant | 6.635 | 0.000 | 3.318 | 0.42 | 10.0 | 0 |
| b UNIFAC-SvL | 8.830 | +5.513 | 7.171 | 0.98 | 11.9 | -2.194 [-4.387, -0.098] |
| c1 weighted Tm | 16.717 | +16.717 | 16.717 | 1.00 | 22.5 | -10.081 [-12.479, -7.565] |
| c2 ridge | 6.434 | +2.659 | 4.547 | 0.40 | 5.3 | 0.202 [-3.177, 2.452] |
| limiting species (new) | 9.229 | +5.916 | 7.573 | 0.98 | 12.8 | -2.594 [-4.783, -0.496] |

Reading: the limiting-species model over-predicts the freezing point for 98 % of dev fuels (mean +5.9 degC) and is significantly worse than the constant (interval below 0). It behaves like the ideal-solubility UNIFAC-SvL baseline (b), slightly worse, so replacing the full liquidus by its single limiting species does not remove the over-prediction bias. Status: no skill shown. Tests: `tests/test_v7_limiting_species.py`.
