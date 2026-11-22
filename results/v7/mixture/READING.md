# Task 11 reading: MolPool mixture retrain (nested Optuna, long-cap encoder)
Rows are the v6 Track C dev rows (already seen in v6), so everything below is exploratory and cannot reopen the v6 verdicts. Locked groups were not scored. Protocol: `PREDECLARATION.md`; outputs `retrain/cv_summary.csv`, `retrain/oof.csv`. Gains are comparator MAE minus v7 MAE (positive = v7 better), paired group bootstrap 90 % CI, Holm over the five heads at 0.10.

| head / prop | v7 MAE | best comparator (MAE) | gain vs best [90 % CI] | Holm p | v6 MolPool v1 MAE | gain vs v1 [90 % CI] | Holm p |
|---|---|---|---|---|---|---|---|
| boehm density | 2.947 | physics (2.873) | -0.074 [-0.924, 0.848] | 1.000 | 3.865 | 0.918 [0.058, 1.734] | 0.110 |
| boehm ln_visc | 0.0593 | desc_gbm (0.0551) | -0.004 [-0.016, 0.005] | 1.000 | 0.0780 | 0.019 [0.007, 0.030] | 0.018 |
| kaust mon | 2.697 | lin_molar (2.912) | 0.216 [-0.456, 0.370] | 1.000 | 3.347 | 0.651 [0.131, 2.523] | 0.110 |
| kaust ron | 2.022 | desc_gbm (1.486) | -0.536 [-1.634, -0.173] | 0.140 | 2.603 | 0.581 [0.344, 1.394] | 0.000 |
| leenhouts mon | 3.695 | desc_gbm (3.084) | -0.612 [-1.174, -0.156] | 0.146 | 2.294 | -1.401 [-1.813, -0.729] | 0.008 |

Reading:
- No head beats its best comparator after Holm. `beats_best` is False for all five; on kaust ron and leenhouts mon the v7 model is worse than the descriptor GBM (intervals below 0 before correction).
- Against the untuned v6 MolPool, v7 is better on boehm ln_visc and kaust ron (Holm p 0.018 and 0.000); boehm density and kaust mon have positive point gains but Holm p 0.110 (not significant at 0.10). On leenhouts mon v7 is significantly worse than v1 (-1.40, Holm p 0.008), so the tuned search does not help there.
- So tuning recovers part of the v1 gap on some heads but does not make the mixture MolPool better than a physics baseline or a descriptor GBM. Status: no skill shown over the best comparator; no change to the v6 Track C verdict.
- Curves: `curves_encoder.csv` and per-head `retrain/<head>_<prop>/curves_*.csv`. Search budget 60 trials per outer fold, timeout 1800 s per fold, device cuda (see `retrain/meta.json`).
- Stated deviation: the v6 pretraining script runs on import, so `pretrain_encoder.py` refits the encoder with the same recipe on the Track B pool with the ysi_holdout masked.
