# T2 — predict(composition, T) wrapper (deliverable)

Status: complete. Subagent rate-limited; finished by direct execution. All 5 tests pass.
Model: `google/antigravity-claude-sonnet-4-6` was the assigned model for T2 (most complex coding task); subagent `ses_e64258fceffe2WBP7qyPUAexMZ` created file, tests ran 5/5.

## File
- `src/fuelblend/api/predict_blend.py` (959 lines, full implementation)
- `tests/test_predict_blend.py` (5/5 PASS)

## API
`predict(composition, T_K=298.15, T_C=None, mole_fractions=None)` → dict with
`viscosity_dyn_mPas`, `density_kg_m3`, `flash_point_C`, `T_K`, `warnings`.

## Paths implemented
- Viscosity: pure GBM ensemble per component → Refutas (or GN fallback) → + mixture-residual GBM (lazy-train on first call if missing)
- Density: pure GBM / fit curves → ideal-volume → + residual GBM
- Flash: Le Chatelier + Raoult (Antoine kPa + LFL, zabetakis T-correct optional); fallback to pure-GBM weighted mix; else `none`
- Conformal: split-conformal 90% quantile on S1 over-prediction (computed from baseline_preds.csv, module-level, ~+2.9 °C depending on data)

## Example outputs (verified, real runs)
(a) n-decane/n-dodecane 50/50 @ 298.15 K: viscosity 1.04 mPa·s (residual_gbm+refutas), density 737 kg/m³, flash 55.6 °C (chat, interval [41.7, 69.4]). Closest measured binary viscosity 1.092 at 298.1 K.
(b) 3-comp (tetradecane/decane/methane): viscosity 1.42 (gn_fallback — methane has no pure-data curve), density 747, flash 21.2 (gbm_mix; methane missing Antoine).
(c) Unknown + decane: global warnings (unknown_compound resolved; mole_fraction_approx); viscosity 0.68 (gn_fallback); flash 37.1 (gbm_mix); density 706.

## Honest limitations (per spec, not hidden)
- Residual GBM only trained when `run/mix_models/` missing; if data split changes, the lazy-train rebuilds from baseline_preds.csv + model_table using `system_id` + `T_K` rounded join (the `_ensure_mix_models()` logic noted in the user's prompt was the active fix; the file already uses `T_K_r4` matching).
- `flash_data_thin`: mixture flash in model_table = 123 binary + 3 multi — bias-only calibrator (T5) is what works at this n.
- Density residual is ~0.7% improvement — applied, noted.
- Pure extrapolation warnings always emit when T > pure-data range ±20 K (correct per WP4/5).
- Unknown/non-HC inputs: degrade gracefully (no crash; method falls back; warning emitted).

## Decisions recorded
- `system_id` / `property` / `T_K_r4` 4-decimal string join for mixture-residual rebuild (matches user's fix direction).
- Antoine convention: `log10(P/kPa) = A - B/(T/K + C)` (used consistently; refit constants in `compounds.csv`).
- LFL: Zabetakis temperature correction applied when `lfl_temp=True`; pure LFL values used when temperature unknown.
- Conformal interval computed once at module load from `baseline_preds.csv` (not per-call) for speed.
