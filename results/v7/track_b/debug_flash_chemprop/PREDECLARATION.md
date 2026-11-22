# Debugging task: why is the in-house D-MPNN 0.36 degC worse than Chemprop on flash_C (Task 7)?

Facts (Task 7, committed 2f86d77): flash_C OOF MAE in-house 17.19, Chemprop 16.83, paired gain CI [-0.75, -0.06].

Differences read from the code (not yet tested):
- H1 learning-rate schedule: in-house uses a constant lr 1e-3 (AdamW, wd 1e-5); Chemprop uses linear warmup then decay to a low final lr within the 150-epoch budget.
- H2 epoch budget: Task 5 showed median best epoch 136 with the long cap; the v1 cap (150, patience 25) may stop the constant-lr model early.
- H3 featurisation: in-house atom features (element, degree, Hs, hybridisation, aromatic, ring, mass/100) differ from Chemprop's default atom and bond features. Not tested here (needs a feature ablation).

Experiment: in-house single-task flash_C, same rows, outer folds, inner splits, hidden 200, depth 4, dropout 0.1, 5 seeds, as Task 7. Variants: (b) scheduler cosine, epochs 150, patience 25; (c) scheduler plateau, epochs 150, patience 25; (d) constant lr, epochs 1000, patience 100. The Task 7 in-house run is variant (a) and is read from its `oof.csv`; Chemprop OOF is read from the same file.
Metric: OOF MAE and paired group-bootstrap 90 % CI of gain (Chemprop minus variant).
Reading rule: a variant whose CI contains 0 (or is above 0) is evidence that the gap comes from the training budget or schedule (H1 or H2); then the Phase 3 search space (which includes scheduler and long caps) already covers it and Task 10 may proceed with that variant named in its predeclaration. If no variant closes the gap, H3 or an implementation difference remains and Task 10 stays held; the next step is a feature ablation. One target, one comparison family; no multiplicity claim is made.
