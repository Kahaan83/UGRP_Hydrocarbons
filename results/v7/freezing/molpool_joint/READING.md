# Task 16 reading: joint pure-Tm plus fuel MolPool for freezing (dev, 52 fuels; exploratory)

Run: `scripts/v7/freezing/molpool_joint.py`, 3 seeds, 30 MC-dropout passes. Numbers: `summary.csv`, `coverage_by_family.csv`. Reading rule: improves only if the 90 % group-bootstrap CI excludes 0 with positive gain.

- LOGO MAE (degC): joint 3.00, fuel-only 3.07, v6 LGBM residual 3.32, ridge 5.83, constant 6.64. Joint gain vs constant 3.63 (CI 2.08 to 5.26), vs ridge 2.83 (1.38 to 4.66): CI excludes 0.
- LOFO MAE: joint 2.80, fuel-only 2.97, ridge 5.23, constant 6.99. Gain vs constant 4.19 (2.65 to 5.82), vs ridge 2.43 (1.28 to 3.67).
- Joint vs fuel-only: gain 0.07 (LOGO, CI -0.21 to 0.34) and 0.17 (LOFO, CI -0.03 to 0.38). The CI includes 0, so the pure-Tm auxiliary data show no demonstrated benefit. The gain comes from the MolPool residual model itself, not from the extra pure rows.
- Joint beats the v6 LGBM residual numerically (3.00 vs 3.32) but this pairing was not bootstrapped; no claim.
- Sign of error under joint: mean over-prediction 1.43 and under-prediction 1.57 degC (LOGO); the unsafe side for freezing is under-prediction.
- MC-dropout SD vs absolute error Spearman: 0.40 (LOGO) and 0.11 (LOFO) for joint; 0.15 and 0.42 for fuel-only. Inconsistent across schemes, no reliable uncertainty signal.
- Coverage with the live margin 13.84 degC (joint): 0.933 for ajftd_posf (15 groups), 1.000 elsewhere (23, 11 and 3 groups).
- Leak check: the encoder was pretrained on pure Tm labels only, no fuel labels; fuel residual base is `pred_e`, as in v6. Hyperparameters were fixed before the run.

Status: dev-only and exploratory. Passes the predeclared dev reading rule against constant and ridge, but it cannot change the real-fuel freezing verdict ("status: no skill shown") without the fresh lock of Task 20.
