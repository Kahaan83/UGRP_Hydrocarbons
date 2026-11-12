# delta_leenhouts_splits (Task 16b)

Declared 2026-11-12 before `delta_leenhouts_splits.py` exists; no result file in this commit. Frozen reference `5a62e20`. Skip-guarded (`extras/`).

The Leenhouts benchmark's own 5-fold split definitions (`splits.zip`: mixture, molecule, pure_to_mixture, random_rowwise) are read in memory (nothing extracted), their
flash-mixture test rows are matched to our 126 hydrocarbon lab rows by (components, mole fractions to 3 decimals, measured value to 0.1 C), and the frozen ideal
Le Chatelier errors (`results/v4/F/flash_lab_loso.csv`, zero fitted parameters) are summarised per split type and fold. Viscosity splits are not used: our lab viscosity cells come
from S1 residual files keyed by system, not by benchmark row (stated, not attempted).
Reading (fixed now): descriptive only. Because the physics model has no fitted parameter, per-fold numbers are a partition of the same 126 errors; they give the
numbers on the benchmark's own folds for comparison with published learned models, not a new test.
