# Task 15 reading: freezing, heavy-tail GP and novelty flag
Protocol: `../PREDECLARATION_task15.md`; families `../families.csv`. Dev rows only (52 groups; holdout not read). Outputs: `summary.csv`, `coverage_by_family.csv`, `novelty.csv`, `predictions.csv`, `features.csv`. Gain = comparator MAE minus model MAE per group, 90 % percentile bootstrap over groups.

| scheme | model | MAE | bias | gain vs constant [CI] | gain vs ridge [CI] |
|---|---|---|---|---|---|
| LOGO | constant | 6.635 | 0.00 | | |
| LOGO | ridge (same features) | 5.827 | +0.87 | +0.81 [-0.77, 2.32] | |
| LOGO | GP (Matern 2.5 ARD) | 6.007 | -0.18 | +0.63 [-0.40, 1.60] | -0.18 [-1.41, 1.05] |
| LOGO | UNIFAC-SvL (b) | 8.830 | +5.51 | -2.19 [-4.39, -0.10] | -3.00 [-4.33, -1.65] |
| LOGO | limiting species | 9.229 | +5.92 | -2.59 [-4.78, -0.50] | -3.40 [-4.75, -2.03] |
| LOGO | v6 LightGBM residual (read) | 3.315 | -0.16 | +3.32 [1.83, 4.91] | +2.51 [1.10, 4.26] |
| LOFO | constant | 6.994 | -0.54 | | |
| LOFO | ridge | 5.233 | +1.80 | +1.76 [0.16, 3.49] | |
| LOFO | GP | 6.701 | +1.15 | +0.29 [-0.53, 1.07] | -1.47 [-2.64, -0.44] |
| LOFO | UNIFAC-SvL (b) | 8.830 | +5.51 | -1.84 [-4.07, 0.30] | -3.60 [-5.54, -1.96] |
| LOFO | limiting species | 9.229 | +5.92 | -2.24 [-4.47, -0.12] | -4.00 [-5.96, -2.35] |

Reading:
- The GP has no skill over the constant (both intervals include 0) and is no better than ridge under LOGO and significantly worse under LOFO (gain vs ridge -1.47 [-2.64, -0.44]). The rule (intervals above 0 against constant and ridge under both schemes) is not met. A ridge on the same five heavy-tail features does beat the constant under LOFO (+1.76 [0.16, 3.49]); that is a cross-family result on 4 families (one with 3 groups) and is exploratory. The v6 LightGBM residual on the isomer-tile base remains the best dev LOGO model in this table (3.315) and is not changed by this task.
- Novelty flag: threshold = 95th percentile of training leave-one-out scores. LOGO flags 4 of 52 fuels (ajftd_posf and crc_jeta), GP MAE 27.5 on flagged versus 4.2 on unflagged; LOFO flags 5 (all ajftd_posf), 23.8 versus 4.9. The flag isolates the fuels with the largest errors (heavy-tail or unusual compositions) but is based on 4 to 5 fuels, so its value is a lead, not a verified rule.
- Coverage of the live 13.8427 C margin by family (`coverage_by_family.csv`, groups per cell stated there): ajftd_posf 0.87 (GP, LOGO, 15 groups; constant 0.93), crc_jeta 1.00 (23), crc_jeta1 1.00 (11), fuellib_reference 1.00 (3). The shortfall sits in the AJFTD family, where the largest errors also sit; coverage of the other families is 1.00 for every model and is uninformative (UNIFAC-SvL over-predicts, so the margin is not binding).
- Status: freezing stays "status: no skill shown" for real fuels in the README; nothing is adopted. Tests: `tests/test_v7_freeze_gp.py`.
