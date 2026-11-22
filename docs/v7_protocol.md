# v7 protocol

1. **Predeclaration.** Each task writes a `PREDECLARATION.md` (comparators, metric, gate) and commits it before any result; results and `READING.md` follow in separate commits. Reading rule: a model improves only if its 90 % paired group-bootstrap interval (4000 draws) excludes 0 with positive gain.
2. **Multiplicity.** Holm correction within each declared comparison family. Exploratory dev gains are not corrected across tasks and carry no claim beyond the dev rows.
3. **Lock guard.** `src/fuelblend/v7/guard.py` lists every lock (YSI holdout, freezing holdout, mixture locks, flash oxygenate lock, `splits/locked.json`) and masks or rejects locked ids. v7 fits use dev rows only; the frozen v5 BASE (`results/v5/BASE/`) and `src/` are not modified by experiment runs (hash-asserted by the retrain driver).
4. **Fresh lock.** `results/v7/lock/` holds 61 Leenhouts flash-mixture systems unseen anywhere, committed before any fit, one primary comparison (HANNA vs UNIFAC), scored once.
5. **CPU and GPU separation.** Tree models and CPU tasks run in `.venv`; GNN and MolPool tasks run in `C:\fuelblend-gpu\.venv` with `FUELBLEND_DEVICE=cuda`. Parallel Optuna processes each use their own `FUELBLEND_RUNS` directory (shared SQLite databases lock).
6. **Status.** All Phase 2 to 4 results are dev-only and exploratory; no v7 model is adopted into `predict_fuel`, and real-fuel freezing stays "status: no skill shown".
