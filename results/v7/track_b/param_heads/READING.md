# Task 18 reading: D-MPNN Antoine and Andrade parameter heads (dev, scaffold GroupKFold; exploratory)

Run: `scripts/v7/track_b/param_heads.py` (GPU, 3 seeds, 5 folds). Numbers: `summary.csv`, `flash_summary.json`, `flash_rows.csv`, `head_params.csv`, `curves.csv`. Predeclaration: `PREDECLARATION.md`.

- log10 P (189 molecules, 4416 points): head MAE 0.474, pooled curve 1.793, ridge 0.220, per-molecule refit (uses own labels) 0.005. Gain vs pooled +1.49 (CI 1.35 to 1.63); gain vs ridge -0.34 (CI -0.40 to -0.28): the head is significantly worse than the descriptor ridge.
- ln mu (136 molecules, 1947 points): head 0.280, pooled 0.726, ridge 0.284, refit 0.033. Gain vs pooled +0.41 (CI 0.33 to 0.50); gain vs ridge -0.038 (CI -0.079 to 0.002): no difference from ridge.
- Flash through the v5 route (pure hydrocarbons, 53 rows after 26 C10+ paraffin rows excluded by Trap 3 and 34 solver failures): MAE head 22.1 degC, pooled 102.8, registry (refit constants, reference) 9.3. Mean over-prediction, the dangerous side: head 16.4, registry 6.7. Head gain vs pooled +80.7 (CI 62.8 to 99.5) but far from the refit reference.
- The ridge comparator predicts log P directly at each T (not constrained to the Antoine form) and wins on vapour pressure; the equation constraint does not pay off at this sample size.

Status: the heads beat only the pooled single curve. No skill shown against the descriptor ridge, and the flash recompute over-predicts by 16 degC on average, so predicted Antoine constants are not usable for flash. Dev only, nothing adopted.
