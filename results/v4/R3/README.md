# r3 re-score: CRC 647 compositions rebuilt from measured hydrocarbon types (Table G.1)

Date: 2026-11-11. Numbers are copied from `readings.csv` and `n_share_sensitivity.csv` in this folder; those files win over this text.

## Headline findings

* **Flash: raw ideal Le Chatelier dev bias +0.44 C (v4 -9.44, r2 +4.22); 13% of dev fuels over-predicted by more than 5 C (r2 53%); locked +2.46 C, 3 of 12.** The r2 sign flip to the unsafe side is mostly gone once the measured classes are used, and the bias stays between -0.26 and +1.25 C over the whole tested range of the one remaining prior.
* Density 15.56 C: dev bias -14.01 kg/m3 (v4 -11.11, r2 -4.69), outside +/-3. It is -13.9 to -14.1 for every tested n-share, so the remaining prior does not explain it. The better r2 number depended on r2's ring share (0.526 of the saturates against a measured mean of 0.26).
* Freezing (ideal SLE, mole): dev bias +4.14 C (r2 +0.32); beating the dev-mean constant is not established (bootstrap of the MAE difference [-0.75, +1.26]). The dev bias runs from -5.9 to +9.8 C as the n-share goes from 0.10 to 0.50, so freezing is controlled by a quantity Table G.1 does not measure.
* Viscosity (raw Refutas, dev, ln): bias -0.400 (v4 -0.301, r2 -0.244), not below 0.10. Recalibrated Refutas dev log-RMSE 0.114 (r2 0.123); structured ridge 0.078 (r2 0.078). The model ranking is unchanged.

The 12 locked fuels were already scored under the v4 and r2 reconstructions. r3 reads no target property and nothing was tuned; the readings below were fixed before any r3 score existed and are judged on dev rows.

## What r3 is

r2 split every fuel's saturates with one fixed prior (n, iso, cyclo) = (0.1366, 0.3374, 0.5260) from 4 GC fuels. r3 uses each fuel's own measured split from CRC Report 647 Appendix G, Table G.1 (PDF pages 125-127), which the original CRC parser skipped:

* ASTM D6379 total aromatics and saturates (vol%), and ASTM D2425 hydrocarbon types: paraffins, cyclo-, dicyclo-, tricycloparaffins, alkylbenzenes, indans/tetralins, indenes, naphthalene(s), acenaphthenes, acenaphthylenes, tricyclic aromatics.
* Parsed by `scripts/v4/parse_crc647_g1.py` into `data/raw/pour_point/crc647_hydrocarbon_types_g1.csv` (1026 rows = 57 fuels x 18 columns; 271 cells printed below the detection limit, kept with their "<" and counted as 0 by r3). Twelve fuels were compared with the page images; all values match. Closure: mono + di = total aromatics within 0.01; D2425 saturate types minus D6379 saturates between -0.2 and +0.7.
* Built by `src/fuelblend/data/reconstruction_r3.py` and `scripts/v4/rebuild_crc_reconstruction_r3.py` into `data/real_fuels/v4/compositions_crc647_r3.csv` (provenance `reconstructed_r3`). Details and the representative table: `reconstruction_report.md`.

Mean mass share over the 57 fuels:

| class | r2 | r3 |
|---|---|---|
| n-paraffin | 0.109 | 0.164 |
| iso-paraffin | 0.268 | 0.406 |
| one-ring cycloalkane | 0.418 | 0.104 |
| two-ring cycloalkane (decalin) | 0 | 0.085 |
| alkylbenzene | 0.083 | 0.145 |
| indan / tetralin | 0.094 | 0.052 |
| olefin | 0.009 | 0.009 |
| fallback lump | 0.020 | 0.035 |

Measured ring share of the saturates (cyclo + dicyclo + tricyclo over all saturate types): mean 0.2575, range 0.027 to 0.559; r2 used 0.526 for every fuel.

## Commands (venv python)

```
parse_crc647_g1.py                 -> data/raw/pour_point/crc647_hydrocarbon_types_g1.csv
rebuild_crc_reconstruction_r3.py   -> data/real_fuels/v4/compositions_crc647_r3.csv, reconstruction_report.md, fallback_detail.csv
eval_bulk_props_v4.py              --compositions data/real_fuels/v4/compositions_crc647_r3.csv --out results/v4/R3/B
eval_flash_v4.py                   --compositions data/real_fuels/v4/compositions_crc647_r3.csv --out results/v4/R3/F
build_viscosity_leaderboard_v4.py  --compositions data/real_fuels/v4/compositions_crc647_r3.csv --out results/v4/R3/V
r3_readings.py                     -> readings.csv (v4, r2, r3)
r3_n_share_sensitivity.py          -> n_share_sensitivity.csv, sens/ns_*/
```

Two evaluator lines were changed so that `reconstructed_r3` is handled exactly like `reconstructed_r2` (`eval_bulk_props_v4.py` allow-list; `build_viscosity_leaderboard_v4.py` direct mode). With default arguments and with the r2 compositions all three evaluators still reproduce `results/v4/{B,F,V}` and `results/v4/R2/{B,F,V}` byte for byte (checked 2026-11-11). The lump row is dropped and the rest renormalised, as in r2. All 57 fuels were scored for every property; viscosity unmappable mass 0.0000.

## Readings declared before any r3 score existed (dev, 45 fuels)

| property | rule | outcome |
|---|---|---|
| density | dev bias within +/-3 kg/m3 | **not met**: -14.01 [-15.26, -12.60] |
| viscosity | abs(raw Refutas dev bias) below 0.10 ln | **not met**: -0.400 |
| freezing | ideal-SLE bias; does the bootstrap interval of MAE(ideal SLE) - MAE(dev-mean constant) exclude 0 | bias +4.14 C; interval [-0.75, +1.26], does **not** exclude 0 |
| flash | sign and size of the raw ideal bias; fraction over-predicted by more than 5 C | +0.44 C; 6 of 45 (0.133) |

## Comparison (pred - meas)

| property | quantity | v4 | r2 | r3 |
|---|---|---|---|---|
| density 15.56 C (kg/m3) | bias dev / locked | -11.11 / -7.93 | -4.69 / -0.47 | -14.01 / -15.00 |
| | MAE dev / locked | 11.81 / 8.39 | 9.11 / 6.43 | 14.04 / 15.00 |
| freezing, ideal SLE mole (C) | bias dev / locked | +14.69 / +13.64 | +0.32 / -1.82 | +4.14 / +2.59 |
| | MAE dev / locked | 14.69 / 13.64 | 4.09 / 2.48 | 5.26 / 3.14 |
| freezing, WAT - dev-fit offset | MAE dev (in-sample) / locked | 4.40 / 3.61 | 4.15 / 2.68 | 3.79 / 2.26 |
| freezing, Coutinho library | bias dev / locked | +4.24 / +3.22 | -15.16 / -17.32 | -10.25 / -11.73 |
| freezing, dev-mean constant | MAE dev (leave-one-fuel-out) / locked | 5.02 / 3.32 | 5.02 / 3.32 | 5.02 / 3.32 |
| viscosity raw Refutas (ln) | bias overall dev / locked | -0.301 / -0.243 | -0.244 / -0.198 | -0.400 / -0.375 |
| | bias 20 C / -20 C / -40 C (dev) | +0.025 / -0.276 / -0.652 | +0.056 / -0.216 / -0.573 | -0.074 / -0.376 / -0.749 |
| recalibrated Refutas | log-RMSE dev / locked | 0.135 / 0.194 | 0.123 / 0.182 | 0.114 / 0.169 |
| structured ridge | log-RMSE dev / locked | 0.083 / 0.150 | 0.078 / 0.150 | 0.078 / 0.156 |
| flash, raw ideal (C) | bias dev / locked | -9.44 / -8.44 | +4.22 / +6.32 | +0.44 / +2.46 |
| | MAE dev / locked | 9.44 / 8.44 | 5.63 / 6.32 | 4.03 / 3.72 |
| | fraction over by > 5 C, dev / locked | 0 / 0 | 0.533 / 0.50 | 0.133 / 0.25 |
| flash, dev-fit offset | MAE dev (leave-one-fuel-out) / locked | 1.71 / 1.92 | 4.12 / 3.29 | 4.10 / 3.50 |
| flash, dev-mean constant | MAE dev / locked | 5.51 / 5.90 | 5.51 / 5.90 | 5.51 / 5.90 |

Flash residual correlations on dev under r3: +0.27 with aromatics vol%, -0.52 with D86 T10, -0.08 with the dropped lump mass (correlations only).

## Sensitivity to the remaining prior (n-paraffin share of the paraffins; nominal 0.288 from the same 4 GC fuels)

Dev rows; locked rows are in `n_share_sensitivity.csv`.

| n-share | density bias | freezing ideal-SLE bias | flash bias | flash fraction over by > 5 C | raw Refutas bias (ln) |
|---|---|---|---|---|---|
| 0.10 | -14.13 | -5.92 | -0.26 | 0.089 | -0.418 |
| 0.20 | -14.06 | +0.55 | +0.11 | 0.133 | -0.408 |
| 0.288 | -14.01 | +4.14 | +0.44 | 0.133 | -0.400 |
| 0.40 | -13.94 | +7.47 | +0.86 | 0.178 | -0.389 |
| 0.50 | -13.88 | +9.81 | +1.25 | 0.222 | -0.379 |

## What can and cannot be said

* The class split matters: three reconstructions give density biases of -11.1, -4.7 and -14.0 kg/m3 and freezing biases of +14.7, +0.3 and +4.1 C. A single reconstructed-composition number for density or freezing is not a stable accuracy claim.
* The r2 density and freezing figures should not be read as "the bias was mostly the old reconstruction": they were obtained with a ring share about twice the measured one.
* The cause of the -14 kg/m3 density bias under r3 is not isolated. It is not the n/iso prior. Not varied here: the carbon-number distribution (taken from D86 through n-alkane boiling points and applied to every class), the representative molecule inside a bin, and the lump handling.
* Flash is the property that improves and is stable to the remaining prior. A margin is still required: 13% (dev) and 25% (locked) of the fuels are over-predicted by more than 5 C.
* The viscosity conclusion is unchanged: raw Refutas fails at low temperature under every reconstruction, and the learned models keep their ranking.

## Caveats

* Compositions are still reconstructed, not measured: Table G.1 gives classes, not molecules or carbon numbers.
* Mean fallback mass 0.0349 (max 0.0817, CRC647_40): naphthalenes 0.0175 (no registry molecule has both pure eta/rho and an Antoine/LFL entry), tricycloparaffins 0.0089, olefin carbon numbers without a representative 0.0071, minor aromatics 0.0014.
* Every two-ring cycloalkane is C10 decalin and every indan/tetralin is indane or tetralin, whatever the carbon number.
* The D2425 columns carry no printed unit; they are treated as mass shares inside each D6379 fraction, and vol% to mass uses the class densities 0.87 / 0.78 / 0.79 g/cm3 as in r2.
* The G.1 D6379 aromatics differ from the D6379 column on page 92 of the same report (both transcribed correctly; correlation 0.97, G.1 higher by 0.84 vol% on average) and from the D1319 values r2 preferred (mean +2.49 vol%). r3 uses G.1, so the aromatic total also changes between r2 and r3.
* The r2 report line "Bicyclic cycloalkane candidates available: 0" is wrong (ring count taken from a temporary RDKit molecule); decalin was available. r2 files were not changed.
* r3 evaluators were run once; v4 and r2 reruns are byte-identical, r3 determinism was not separately tested.
* `scripts/v4/reproduce_all.py` and `SHA256SUMS_v4.txt` do not cover r3.
