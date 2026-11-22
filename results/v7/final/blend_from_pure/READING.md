# Reading: blend properties from pure-compound predictions

Rule applied as predeclared (Holm p below 0.10, gain above 0, interval above 0). Dev rows only, locked systems removed, no row dropped for failed predictions. Numbers: `summary.csv`.

| task | rows | systems | MAE model | MAE constant | gain (90% CI) | Holm p | reading |
|---|---|---|---|---|---|---|---|
| dcn_leenhouts | 251 | 31 | 6.90 | 21.57 | 14.68 (11.28 to 17.90) | 0.000 | usable, exploratory |
| ron_kaust | 90 | 5 | 6.38 | 13.56 | 7.18 (4.28 to 9.81) | 0.006 | usable, exploratory (5 systems only) |
| mon_kaust | 90 | 5 | 6.49 | 10.41 | 3.92 (2.36 to 5.88) | 0.006 | usable on this source (5 systems only) |
| mon_leenhouts | 135 | 32 | 10.60 | 10.05 | -0.55 (-4.34 to 4.07) | 0.744 | no skill shown |

Notes.
* Derived cetane: the blend error (6.90) equals the error of the same rule with measured pure labels (6.87 on 242 rows), so the linear mole-fraction rule, not the pure model, is the limit.
* MON has opposite readings on its two sources: usable on KAUST (up to 4 components, 5 systems), no skill on Leenhouts (up to 34 components, 32 systems) where even measured pure labels give 7.50 against the pure-model blend 10.89 on the same rows. The predeclaration scored tasks, not properties; for the property-level decision to add keys to `predict_fuel`, MON is treated as unsupported because its evidence does not agree across sources. This property-level choice was made after seeing the results.
* RON rests on 5 systems from one source: usable by the rule, thin by any other standard.
* Cetane has no mixture data: not tested, not added.
* Dev rows, no conformal margin, no lab validation. No v5 or v6 verdict is reopened.
