# Task 10 reading: Track B D-MPNN retrain (nested Optuna, GPU)

Run: 60 trials per outer fold, 5-model final ensembles, 8836 s total, 0 OOM retries. Rule from PREDECLARATION.md applied unchanged. Dev rows only; the search was tuned on these rows, so gains are exploratory.

- Beats v1 (Holm p < 0.10, positive gain): tm_K only (gain 3.20 K, Holm p 0.015). No evidence of improvement on the other nine targets. ysi has a large point gain (5.64) but Holm p 1.0, so no claim.
- GNN earns its place vs descriptor trees: dcn only (gain 1.70, Holm p 0.044).
- Trees better than the retrained GNN (Holm p < 0.10, negative gain): flash_C (RF 11.61 vs 17.22, p 0.010), dhfus_kJ_mol (RF 7.23 vs 10.98, p 0.032), logvisc_mPas_293 (RF 0.148 vs 0.223, p 0.056).
- Other targets (cetane, density, mon, ron, tm_K, ysi): no evidence either way against trees.

Numbers: `cv_summary.csv`; chosen parameters: `best_params.csv`; curves: `curves.csv`, `curves.png`.
