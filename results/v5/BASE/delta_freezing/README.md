# delta_freezing (Task 14)

Declared 2026-11-12 before `delta_freezing.py` exists; this commit contains no result file. Frozen reference `5a62e20`.

(a) F1x (S5 F1 with the complete Coutinho n-alkane table, column `F1x` in `S5/s5_predictions.csv`) as a declared variant, scored with `group_stats`
against the leave-one-group-out constant, on AV-33-22, AJFTD and pooled (headline rows). Check: pooled MAE reproduces S5's 3.07 C (n = 49) within 0.02.
(b) `data/processed/v2/wax_dauphin_coutinho.csv` (57 wax-appearance temperatures) as a second T1 SLE source: ideal solid-liquid equilibrium
(max over components of 1/(1/Tm - R ln x / dHfus)) with the SLE2 `constants_extended.csv` values, at two cumulative tiers (v4_verified+two_compilations; +single_compilation).
A system is scored only if every component has Tm and dHfus at that tier; others are counted. Method label `wax_appearance`; reported apart from liquidus data, never in a liquidus cell.
(c) D86 lead: on the AV-33-22 fuels with a measured freezing point, ridge (leave-one-group-out, alpha by inner LOGO) on D86 T90, D86 FBP and measured n-paraffin mass share,
against the LOGO constant.
Reading (fixed now): "F1x beats the constant on AV-33-22 and on the pooled set; AJFTD stays not evaluable below 20 groups".

## Result 2026-11-12

Run with `scripts/v5/base/delta_freezing.py`; `delta.csv`, `rows_wax_appearance.csv`, `rows_d86_ridge.csv`, `readings.json`.
* (a) Reading PASSES. F1x pooled MAE 3.068 C (n = 49, reproduces S5's 3.07), beats the LOGO constant pooled and on AV-33-22 (bootstrap interval above 0);
  AJFTD: 12 groups, no skill shown, below 20 groups as declared.
* (b) Wax appearance (not liquidus): with v4_verified+two_compilations constants only 9 of 57 systems have constants for all components: MAE 3.87 K, all under-predicted, no skill
  shown against the constant. Adding single_compilation constants (34 systems) gives MAE 10.8 K, worse than the constant. Ideal SLE under-predicts wax appearance here.
* (c) D86 lead confirmed on 34 AV-33-22 fuels (one per S3 group, censored value excluded): ridge on D86 T90, FBP and n-paraffin share, LOGO, MAE 2.98 C against 6.49 C for the constant,
  interval +2.06 to +4.87 C (beats_trivial). Exploratory: one source, features chosen after the audit had seen the lead.
  (The composition file's `double_count_group` column is one constant for all AV-33-22 rows; the S3 groups from `t2_fuels.csv` were used.)
