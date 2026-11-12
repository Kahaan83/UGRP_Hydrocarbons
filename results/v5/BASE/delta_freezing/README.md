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
