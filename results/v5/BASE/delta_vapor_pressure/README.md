# delta_vapor_pressure (Task 16a)

Declared 2026-11-12 before `delta_vapor_pressure.py` exists; no result file in this commit. Frozen reference `5a62e20`.

Data: AV-33-22 vapour pressure (266 rows, 38 fuels, 273-393 K, ASTM D6378). Model: Raoult + Antoine (the flash registry's constants) on the mapped measured composition,
P = sum x_i Psat_i(T), mole fractions from mass; error ln(pred/meas). Rows with mole coverage of Antoine constants below 0.90 or measured value <= 0 are excluded and counted.
Caveat stated in advance: D6378 is a triple-expansion vapour pressure at V/L = 4, not a bubble point, so Raoult is expected to over-predict for a wide-boiling fuel.
Use: per fuel, the ln VP error at 313.15 K (close to the flash temperatures) is set beside the frozen flash error, to separate Psat error from LFL error.
Reading (fixed now): "if the flash error is mainly Psat error, flash error and ln VP error at 313 K correlate positively (Spearman > 0.5 over the AV-33-22 groups)".
