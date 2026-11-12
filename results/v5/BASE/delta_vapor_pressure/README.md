# delta_vapor_pressure (Task 16a)

Declared 2026-11-12 before `delta_vapor_pressure.py` exists; no result file in this commit. Frozen reference `5a62e20`.

Data: AV-33-22 vapour pressure (266 rows, 38 fuels, 273-393 K, ASTM D6378). Model: Raoult + Antoine (the flash registry's constants) on the mapped measured composition,
P = sum x_i Psat_i(T), mole fractions from mass; error ln(pred/meas). Rows with mole coverage of Antoine constants below 0.90 or measured value <= 0 are excluded and counted.
Caveat stated in advance: D6378 is a triple-expansion vapour pressure at V/L = 4, not a bubble point, so Raoult is expected to over-predict for a wide-boiling fuel.
Use: per fuel, the ln VP error at 313.15 K (close to the flash temperatures) is set beside the frozen flash error, to separate Psat error from LFL error.
Reading (fixed now): "if the flash error is mainly Psat error, flash error and ln VP error at 313 K correlate positively (Spearman > 0.5 over the AV-33-22 groups)".

## Result 2026-11-12

246 rows scored, 20 excluded (measured 0 kPa or Antoine mole coverage < 0.90). Raoult + Antoine UNDER-predicts D6378 vapour pressure: ln bias -0.35 (273 K), -0.68 (293 K),
-0.23 (313 K), -0.21 (333 K), -0.20 (353 K), -0.15 (373 K), +0.10 (393 K). This is the opposite of the declared expectation (V/L = 4 should make Raoult over-predict).
Reading FAILS: Spearman(flash error, ln VP error at 313 K) = -0.13 over 35 AV-33-22 groups. The frozen AV-33-22 flash bias is -1.06 C (predicted flash too low, i.e. too volatile),
while the predicted vapour pressure is too low (too little volatile). A Psat error of this sign would push flash up, not down, so the AV-33-22 flash under-prediction is not
explained by Psat; the LFL side (or the light end of the mapped composition) is the remaining suspect.
