# Task 17 predeclaration: HANNA gamma vs UNIFAC and ideal in the v5 flash route (dev rows)

Date 2026-11-22. HANNA: MIT licence code and weights (tspecht93/hanna), ChemBERTa-77M-MTR embeddings, pretrained, no refit, binary only.
Rows: Leenhouts hydrocarbon binary flash-point rows of `data/unified/measurements_dedup.csv` (kind binary, is_hc, no conflict), excluding systems in any system-id lock (guard). Pure dev rows; the fresh lock (`results/v7/lock/`) is not touched.
Arms (same Le Chatelier route, Antoine and LFL, only gamma differs): ideal (gamma 1, production), UNIFAC (`make_unifac_gamma_fn`), HANNA. Rows where any arm fails are dropped for all arms and counted.
Metric: per-system MAE in degC (rows averaged within system), and mean one-sided over-prediction (pred - measured, positive part). Paired group bootstrap by system, 4000 draws, 90 % percentile CI.
Primary: HANNA vs ideal. Secondary: HANNA vs UNIFAC. Holm over these two.
Gate (dev reading): MAE gain CI above 0 and upper CI bound of (HANNA minus comparator) mean over-prediction below 0.5 degC.
Caveat: HANNA's training data (Dortmund Data Bank binaries) may include some of these systems; a dev gain would not show generalisation. The fresh lock is the only scoring that could.
Exploratory: result cannot change any v5 verdict or `predict_fuel` (production stays ideal).
