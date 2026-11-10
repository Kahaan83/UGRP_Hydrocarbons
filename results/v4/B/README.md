# v4 Task 12: density, net heat of combustion, freezing point on the v4 corpus

Produced by `python scripts/v4/eval_bulk_props_v4.py` (about 13 s, deterministic: bootstrap seed 20261110).
Every number below was read from `density_summary_v4.csv`, `hoc_summary_v4.csv` and `freezing_summary_v4.csv`
(per-fuel detail in `density_v4.csv`, `hoc_v4.csv`, `freezing_v4.csv`). Test: `tests/test_v4_bulk_props.py`.

## What was evaluated

- **Fuels**: the 57 CRC647 aviation fuels that have a `raw_file` target in `data/real_fuels/v4/jet_properties_v4.csv` and a
  composition in `compositions_v4.csv`. **Dev = 45 fuels, locked = 12 fuels** (`splits/v4/real_fuels_v4_splits.json`
  `fuel_role`). Dev and locked are reported separately; nothing was selected, fitted or tuned on either slice.
- **Not evaluated**: the GC jet fuels in `included_jet_fuel_ids()` (JetA_API14avg, UDRI_HTA_repfuel, POSF-10325, POSF-13511,
  POSF-13516) have class x carbon-number compositions with no molecule list and no bulk raw_file target these models can use.
  POSF-10325 is the n=1 headline set and is excluded from every statistic (the test asserts this). No diesel or gasoline fuel
  qualifies (SD-02; the gasoline blends have no density/HoC/freezing targets).
- **Targets** (all `raw_file`, CRC 647 tables): density at 15.56 C (ASTM D4052, 57 fuels) plus Solartron densities at other
  temperatures (398 rows in total); net heat of combustion (57); freezing point, D5972 (57).
- **Compositions are reconstructed** (ruling R10: derived from measured aromatics plus D86, not measured by GC). Every error
  below is therefore **reconstruction error plus physics error**; the two cannot be separated with these data.
  The non-lumped molecules are renormalised, which assumes the dropped unresolved tail (mean 0.8 % of mass, max 2.3 %) behaves
  like the rest of the fuel.
- **Perdeuterated alkanes**: the reconstructed compositions contain perdeuterated n-alkane SMILES (on average 11.5 % of the mass,
  max 20.5 %). They are treated as the corresponding protium alkanes (isotope labels stripped before property lookup and before
  molecular weights are computed). This is a judgement call; it is not validated.

## Models (same as v3, same metric set: MAE, RMSE, MAPE, bias = pred - meas, fuel-resampled bootstrap 95 % CI)

| property | model | pure-component inputs (no fuel data used) |
|---|---|---|
| density | ideal volume mixing, `rho = 1/sum(w_i/rho_i(T))` | tier 1: linear-in-T fit of ThermoML liquid pure densities (`data/processed/v2/thermoml_clean_hc.csv`, robust MAD trim); tier 2: class + 1/nC regression on the tier-1 fits. **No constant fallback.** On the 15.56 C rows, 80.1 % of mass (mean) has a tier-1 density, 19.9 % tier 2 |
| HoC | mass-weighted group-additivity GCV (v2 table `gcv_v2_MJ_kg.csv`), converted to NCV with `mixture_heat_of_combustion_v2` | v2 table (Tier A measured formation enthalpies / Tier B group additivity); on average 53.4 % of mass is Tier A |
| freezing | `src/fuelblend/thermo/coutinho_wax.py` as shipped (Tm / dHfus for n-alkanes C6+ and four aromatics) | library table |

Freezing variants reported: library ideal WAT (the library puts **mass** fractions inside the Schroeder-van Laar log, as in
v3 R5), `WAT - 14.58` offset, Coutinho eutectic, and an ideal SLE with true **mole** fractions computed in this script
(not in the library).

## Results

### Density, kg/m3 (ideal volume)

| scope | slice | N | MAE [95 % CI] | RMSE | MAPE % | bias [95 % CI] |
|---|---|---|---|---|---|---|
| 15.56 C D4052 (v3-comparable) | dev | 45 | 11.81 [9.65, 14.30] | 14.36 | 1.46 | -11.11 [-13.98, -8.54] |
| 15.56 C D4052 (v3-comparable) | locked | 12 | 8.39 [5.37, 11.92] | 10.14 | 1.04 | -7.93 [-11.72, -4.51] |
| 15.56 C D4052 (v3-comparable) | dev+locked (not a headline) | 57 | 11.09 [9.11, 13.17] | 13.58 | 1.37 | -10.44 [-12.78, -8.23] |
| all temperatures (rows, not fuels) | dev | 315 | 11.87 [10.97, 12.78] | 14.42 | 1.47 | -11.12 [-12.13, -10.10] |
| all temperatures (rows, not fuels) | locked | 83 | 8.84 [7.53, 10.35] | 10.89 | 1.10 | -8.11 [-9.80, -6.65] |
| 15.56 C, >=80 % of mass with measured pure density | dev | 21 | 11.37 [7.53, 16.06] | 15.08 | 1.41 | -10.05 [-15.36, -5.62] |
| 15.56 C, >=80 % of mass with measured pure density | locked | 9 | 6.47 [4.07, 9.06] | 7.57 | 0.81 | -5.85 [-8.90, -2.74] |

Ideal volume mixing under-predicts density by about 8 to 11 kg/m3 on average (the CI excludes zero in every slice). The
sign is systematic, so the error is not noise. The cause is not isolated: candidates are the reconstructed class split, the
tier-2 class densities (20 % of mass), the dropped tail, and the deuterated-stand-in treatment; the subset whose mass is mostly
measured pure densities does not remove the bias (dev -10.05, locked -5.85). Do not read the locked MAE as better than dev:
the two CIs overlap.

### Net heat of combustion, MJ/kg (group additivity)

| model | slice | N | MAE [95 % CI] | RMSE | MAPE % | bias [95 % CI] |
|---|---|---|---|---|---|---|
| v2-comparable (constant 13.8 wt % H) | dev | 45 | 0.191 [0.151, 0.230] | 0.237 | 0.44 | -0.106 [-0.168, -0.043] |
| v2-comparable (constant 13.8 wt % H) | locked | 12 | 0.154 [0.099, 0.214] | 0.186 | 0.36 | -0.134 [-0.207, -0.065] |
| v2-comparable (constant 13.8 wt % H) | dev+locked | 57 | 0.183 [0.148, 0.218] | 0.227 | 0.43 | -0.112 [-0.161, -0.059] |
| composition-derived H | dev | 45 | 0.153 [0.123, 0.185] | 0.188 | 0.36 | -0.057 [-0.109, -0.003] |
| composition-derived H | locked | 12 | 0.111 [0.071, 0.155] | 0.135 | 0.26 | -0.085 [-0.143, -0.029] |
| measured H (diagnostic only; uses a measured descriptor as an input) | dev | 45 | 0.238 [0.195, 0.281] | 0.282 | 0.55 | -0.186 [-0.247, -0.123] |
| measured H (diagnostic only) | locked | 12 | 0.209 [0.149, 0.267] | 0.235 | 0.48 | -0.209 [-0.267, -0.149] |

The v2-comparable variant is the like-for-like row. Using measured hydrogen makes the error worse, not better, which
is consistent with part of the apparent v2 accuracy being error cancellation between the constant-H assumption (13.8 wt %
against a measured mean of 14.17) and the GCV bias, rather than correct hydrogen physics (not isolated here). The
composition-H variant (mean 13.57 wt % H) has the lowest error, which fits the same reading and is not evidence it is the
better physics.

### Freezing point, C (negative error = prediction colder than measured = optimistic, the unsafe direction)

| model | slice | N | MAE [95 % CI] | RMSE | bias [95 % CI] | n optimistic > 3 C | min error |
|---|---|---|---|---|---|---|---|
| library ideal WAT (mass fraction in ln) | dev | 45 | 18.04 [16.19, 20.02] | 19.25 | +18.04 [+16.19, +20.02] | 0/45 | +4.26 |
| WAT - 14.58 | dev | 45 | 5.33 [3.87, 7.13] | 7.56 | +3.46 [+1.61, +5.44] | 5/45 | -10.32 |
| Coutinho (library) | dev | 45 | 5.57 [4.16, 7.25] | 7.72 | +4.24 [+2.50, +6.12] | 5/45 | -8.05 |
| ideal SLE, mole fractions | dev | 45 | 14.69 [12.93, 16.57] | 16.04 | +14.69 [+12.93, +16.57] | 0/45 | +0.94 |
| library ideal WAT | locked | 12 | 16.94 [14.32, 19.72] | 17.60 | +16.94 [+14.32, +19.72] | 0/12 | +8.75 |
| WAT - 14.58 | locked | 12 | 4.12 [2.49, 6.39] | 5.33 | +2.36 [-0.26, +5.14] | 1/12 | -5.83 |
| Coutinho (library) | locked | 12 | 4.25 [2.57, 6.54] | 5.53 | +3.22 [+0.81, +5.89] | 1/12 | -4.29 |
| ideal SLE, mole fractions | locked | 12 | 13.64 [11.22, 16.15] | 14.32 | +13.64 [+11.22, +16.15] | 0/12 | +5.75 |

Both ideal SLE variants predict the freezing point far too warm (they are not freezing-point models: the ideal WAT is the
first wax appearance). The two corrected models land at 4 to 6 C MAE.

## v3 vs v4, side by side (v3 values from the committed v3 CSVs; the sets differ, see the last column)

| property | v3 result | v4 result (all 57) | why not identical |
|---|---|---|---|
| density, jet | `results/v3/G10`: MAE 16.86 kg/m3, 2.12 %, N=82 jet fuels | 11.09 kg/m3, 1.37 %, N=57 | v3 used 77 of 82 placeholder-constant jet densities as targets and an 800 kg/m3 fallback for any component without a pure density (ledger, CLM-25); v4 uses only raw CRC647 densities and no constant fallback. Not a model improvement |
| HoC (group additivity) | `results/v2/C1`: MAE 0.1792 [0.1462, 0.2154], bias -0.1039, N=57 | 0.183 [0.148, 0.218], bias -0.112, N=57 | same 57 CRC fuels and the same constant-H model; difference reflects the v4 composition table (isotope-normalised, reconstructed) vs the v2 bulk-descriptor route |
| freezing, real fuels | `results/v3/R5` (multi-component rows): offset MAE 4.99 / bias +3.14, Coutinho MAE 4.22 / bias +1.94, N=57 | offset 5.07 / +3.23, Coutinho 5.29 / +4.03, N=57 | same fuels; Coutinho gets about 1 C worse and 2 C warmer than v3 on the v4 compositions. The v3 R5 CSV pools these with 56 simple blends that are not in the v4 corpus and are not re-run here |

## Caveats

1. **Reconstructed compositions** (above): the density and HoC errors measure the composition reconstruction as much as the
   mixing rule.
2. **Freezing: the models do not beat a constant.** Predicting the mean measured freezing point of the slice gives an MAE of
   4.90 C (dev, SD 6.5 C) and 3.32 C (locked, SD 4.1 C), computed from `freezing_v4.csv` (in-sample constant, so a generous
   benchmark). Coutinho is 5.57 / 4.25 C. With N=45/12 there is no evidence that either corrected model carries composition
   information on these fuels. Six of 57 fuels are over-optimistic by more than 3 C (worst -10.3 C for the offset model).
   This bears on any design margin: do not use these MAEs as a margin.
3. **The library `coutinho_wax.py` has a basis problem.** It puts mass fractions inside the Schroeder-van Laar logarithm
   (parameter `basis` is accepted but never used). The mole-fraction variant is computed in this script only for comparison;
   the library was not modified (Task 12 scope). Coutinho's three parameters are the v3 empirical constants, not refit.
4. **Density fit domain**: tier-1 fits are linear in T over 250 to 380 K data; extrapolation up to 30 K beyond a compound's data
   range is allowed; the tier-2 class regression has no uncertainty estimate. Pyrene-like and other solid-at-ambient components
   in the reconstruction have no sensible liquid density at 15.56 C.
5. **Perdeuterated stand-ins** (11.5 % of mass on average) are treated as protium alkanes. Across the 57 fuels the stand-in mass
   fraction correlates +0.50 with the density error; this is an observation, not a diagnosed cause.
6. **HoC**: Tier A (measured formation enthalpy) covers only 53 % of mass; the rest is the v2 group-additivity table.
7. **Small N**: locked N=12 gives wide CIs (for example HoC bias CI width 0.14 MJ/kg). Dev and locked are never pooled for a
   headline. The "dev+locked" rows exist only to match the v3 N=57 comparison.
8. n=1 GC headline fuel (POSF-10325) and all GC fuels are absent from every statistic here.

## Claims-register rows these results would update (register not edited)

- **CLM-25** (G10 density: REVISED-PENDING): re-scored on raw CRC647 only, MAE 11.09 kg/m3 (dev 11.81, locked 8.39), systematic
  bias about -10 kg/m3.
- **CLM-28** (density margin 2.50 kg/m3 and the other margins): the density margin of 2.50 is far below the v4 MAE (8 to 12
  kg/m3) and the bias is one-sided, so it cannot stand as a calibrated margin on this evidence.
- **CLM-18** (group additivity as the recommended HoC route; MAE 0.1792): reproduced at 0.183 on v4 inputs; the measured-H
  diagnostic (0.232) weakens the "first-principles" reading.
- **CLM-07** (ridge residual on HoC): not re-run here (not in the Task 12 model list).
- **CLM-16** (offset / Coutinho for freezing): v4 real-fuel MAE 5.07 / 5.29 C, but neither beats a constant (caveat 2).
