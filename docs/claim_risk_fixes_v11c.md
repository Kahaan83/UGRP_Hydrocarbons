# Claim-risk addendum (post v11c)

No BASE file or number was changed. All figures below are read from committed results or recomputed from committed row files.

## 1. Density margin independence
Use the dev out-of-fold margin as the reported margin, not the locked margin.
- v5 ML: dev OOF margin 5.0671 kg/m3. Locked coverage under this dev margin: 38 of 41 systems (0.927), independent of the locked set (`results/v5/BASE/delta_ml_density/run_log.txt`, R-MLD-4).
- Tuned ML (v11b): dev OOF margin 4.5424 kg/m3 (`docs/final_numbers_v11.md`). Locked coverage under this margin: 38 of 41 systems (0.927), computed by `scripts/v11/density_locked_coverage.py` (margin reproduces 4.5424; output `results/v11/claim_risk/density_locked_coverage.json`). Margin from dev only, model fit on dev only.
- The locked margins (4.4248, 4.4241) are descriptive only and must not be presented as calibration.

## 2. Freezing coverage and bias
Real-fuel headline set, UNIFAC designated route (`results/v5/BASE/delta_freeze_t2_unifac/rows.csv`, 50 fuels with a measurement):
- Mean prediction error +5.88 C (predictions high); MAE 8.42 C against a constant baseline of 6.28 C on these rows. Status: no skill shown.
- 94 % of fuels fall within the 13.843 C margin. This reflects margin width, not skill, and must not be reported as coverage evidence.
- One-sided 90th percentile of over-prediction error: 11.96 C.
- Diagnostic only: subtracting the in-sample mean bias gives MAE 3.87 C. This uses the test fuels to estimate the bias, so it is not a model result and must not be cited as performance.

## 3. Locked-set reuse and post hoc hybrid
Not fixable with current data (no fresh holdout). Paper text: locked sets were scored a second time in v11 to v11c; the v11c warm-bin viscosity hybrid is labelled post hoc.
Nested evidence that uses dev systems only (`results/v11/nested/visc_summary.json`, nested 5-fold GroupKFold on system_id, 264 dev systems, warm bins T >= 283.15 K): macro MAE 0.0599 (v5) to 0.0578 (tuned), mean gain 0.00206 ln against a 0.002 threshold, 95 % CI [0.0002, 0.0041], Holm-adjusted p 0.024. The gain is therefore real on dev but only just above the pre-declared threshold. The whole-model variant failed the cold-bin margin guard (T < 263 K margin 0.2875 to 0.3864); restricting the tuned model to the warm bins is the post hoc step. A nested check of that restriction itself was not run: the warm-only choice was made after seeing the dev guard result, so its dev evidence is selection-affected.

## 4. D86 limits
Not fixable without the primary ASTM D86 / D1655 text (paywalled). Paper text: certified designs depend on D86 limits taken from secondary sources only; certification is model-level.
