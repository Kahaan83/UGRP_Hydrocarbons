# delta_d86_rayleigh: ASTM D86 as an ideal single-stage (Rayleigh) batch distillation

Declared 2026-11-15 BEFORE any declared reading in this folder was computed. Base commit 1e1c28c (v5).
Gates, spec limits, spec `verification:` flags, the judge and `configs/uncertainty_v5.yaml` are not changed. The WS7 margins file
`results/v5/BASE/delta_ws7_design_realism/d86_margins.json` (the default of `blendstock_v5.load_d86_margins`) is not touched; new margins go to
this folder only. Adopting any of them (designating the Rayleigh D86 model) is a **user decision**.

## Why (orientation, no declared reading; what shaped the declarations)

A probe (not a declared reading; same Antoine, MW and rho as below) pushed the stored pure-pool recipes through both models:

| recipe | Engler T10 / T90-T10 (C) | Rayleigh T10 / T90-T10 (C) |
|---|---|---|
| WS7 pool C Jet A (lam 0.25) | 175.7 / -5.7 (non-monotone) | 175.0 / 7.3 |
| WS7 P0 Jet A | 125.7 / 49.6 | 169.4 / 16.8 |
| P0x-a Jet A | 134.6 / 54.3 | 175.4 / 18.5 |
| P0x-b Jet A | 136.6 / 55.7 | 175.0 / 17.1 |

The Engler route interpolates a TBP step curve linearly, so 4-5 wt % isooctane (Tb 99 C) drags TBP T10 to about 110 C and D86 T10 to about 125-135 C,
which inflates T90-T10 to about 50 C. In an equilibrium batch distillation a few % isooctane is gone well before 10 vol % is recovered. For narrow
pure-compound recipes with a light minor component the Engler slope quantities are therefore biased on the unsafe side (T50-T10 and T90-T10 are
minimum limits). This is the motivation; the declared readings below decide whether the Rayleigh model is better on measured fuels.

## Declared model (fixed now; zero fitted parameters) - `src/fuelblend/baselines/distillation_rayleigh.py`

* Ideal Raoult (gamma = 1), P = 101.325 kPa, Antoine `log10(P/kPa) = A - B/(T/K + C)` from the default flash registry
  (`flash_unifac._lookup(smiles, None)`, the same sets `blendstock_v5.tb_C` inverts for the Engler route).
* Liquid moles n_i = w_i / M_i, M_i from RDKit `MolWt` of the registry SMILES.
* Single equilibrium stage, vapour removed at once: dn_i/dN = y_i = x_i P_i(T_b)/P, T_b the bubble point of the residual liquid (Newton).
* Recovered volume V_d = sum (n_i0 - n_i) M_i / rho_i(288.15 K), rho from the pure grid (as the Engler route's volume fractions).
* T_q = bubble point of the residual liquid at V_d/V_0 = q/100; q = 10, 50, 90 and the FBP cut 99 (`FBP_CUT_PCT`, as Engler). T50-T10, T90-T10 as differences.
* `solve_ivp` RK45, rtol 1e-9, exact event location at the cuts (smooth in the recipe).
* Domain rule unchanged: TBP T90-T10 inside the fit fuels' TBP range and a monotone D86 curve (Rayleigh is monotone by construction).
* No offset, stage count or other correction is fitted after seeing residuals. A fitted correction, if ever tried, is a separate reading with
  leave-one-fuel-out residuals.

## Declared readings

* **R-DR-1 (unit tests).** `tests/test_delta_d86_rayleigh.py`: (a) one compound -> every cut equals its Antoine Tb; (b) a binary with equal
  molar volumes and constant relative volatility (Antoine sets built so P2/P1 is constant) matches the analytic Rayleigh equation to 0.05 C;
  (c) `blendstock_v5.add_d86_models` with the new `d86_model` argument at its default (`"engler"`) gives identical quantities, margins and
  reasons to the current code on the WS7 pool C stored recipe; (d) Rayleigh cuts are monotone.
* **R-DR-2 (margins on measured fuels).** The 65 fuels of `delta_ws7_design_realism/d86_fit_rows.csv` (same data rules, same exclusions, same
  mapping and coverage), predicted with Rayleigh; margins with the unchanged `blendstock_v5.fit_d86_margins` (alpha 0.10, pooled and per-source,
  certifying = max(pooled, per-source with >= 10 fuels), gate >= 30 fuels with >= 10 per source). Written to `d86_rayleigh_margins.json` and
  `d86_rayleigh_fit_rows.csv`. **Adoption rule, per constrained quantity (T10, FBP, T50-T10, T90-T10):** Rayleigh is preferred for that quantity
  only if its certifying margin is smaller than Engler's AND its MAE is not larger than Engler's. Otherwise Engler stays for that quantity.
  Also reported: bias, MAE, MAE of the constant predictor, per source.
* **R-DR-3 (stored recipes).** Engler and Rayleigh D86 quantities for every stored WS7 / delta_smoke_branch_opt recipe (P0, P0x-a, P0x-b, P1,
  P2, P3, pool C), with the slack of each D86 constraint under each model's own certifying margin.
* **R-DR-4 (designs).** Only if R-DR-2 prefers Rayleigh for at least one quantity: pools P2, P3, P0x-b x ASTM Jet A and BIS Jet A-1, the
  delta_smoke_branch_opt protocol unchanged (D86 on, 10 starts, seed 42, lambda 0.25, synthetic cap 0.50, `alt_standin` on), with the D86 models
  per the R-DR-2 adoption rule. Reported as "under D86 variant Rayleigh (pending user adoption)". **Stated expectation:** P0x-b fails T90-T10
  (pure recipes are narrow); P2/P3 stay uncertified (represented-blendstock rule); 0 certified.

## Results

(appended 2026-11-15 after the runs; nothing above was edited)

### R-DR-1: pass
`tests/test_delta_d86_rayleigh.py` (single compound = Tb; analytic binary Rayleigh to 0.05 C; monotone cuts; default `d86_model="engler"`)
plus the WS7 and delta_smoke_branch_opt feature-off tests: 20 passed.

### R-DR-2 (65 fuels, `compare_margins.csv`, `d86_rayleigh_fit_rows.csv`, `d86_rayleigh_margins.json`)
| quantity | direction | Engler margin | Rayleigh margin | Engler MAE | Rayleigh MAE | Engler bias | Rayleigh bias | constant MAE | Rayleigh preferred |
|---|---|---|---|---|---|---|---|---|---|
| T10 (max 205) | under | 7.93 | 0.00 | 4.52 | 11.63 | -0.69 | +11.61 | 8.19 | no |
| T50 | (no limit) | 14.96 | 8.39 | 8.75 | **6.18** | -2.27 | +1.09 | 11.45 | - |
| T90 | (no limit) | 26.66 | 9.76 | 13.06 | **7.82** | -9.68 | -0.34 | 15.31 | - |
| FBP (max 300) | under | 3.26 | 6.10 | 15.75 | 9.54 | +15.11 | +7.61 | 16.52 | no |
| T50-T10 (min 15) | over | 18.68 | 1.24 | 7.67 | 10.88 | -1.58 | -10.52 | 9.32 | no |
| T90-T10 (min 40) | over | 18.25 | 5.31 | 12.56 | 13.79 | -8.98 | -11.95 | 16.53 | no |

Per-source Rayleigh MAE: AV-33-22 T10 9.97, T50 3.82, T90 4.51, FBP 6.89; AJFTD T10 13.55, T50 8.94, T90 11.68, FBP 12.64.

Reading:
* **By the declared rule Rayleigh is not adopted for any constrained quantity** (each fails "smaller margin AND MAE not larger"). R-DR-4 is
  therefore not run. This is a negative result for designation.
* Rayleigh is the better model of the curve's middle and back end: T50 MAE 6.2 vs 8.7 C, T90 7.8 vs 13.1 C (bias -0.3 vs -9.7), FBP 9.5 vs 15.8.
* Its weak point is the front end: T10 reads 11.6 C high on average (measured D86 vapour temperatures near the start run below the equilibrium
  bubble point of the liquid). That one bias makes both slopes about 11-12 C too small: safe-side (hence the small slope margins) but less accurate.

### R-DR-3 (stored recipes, `stored_recipes.csv`)
| recipe | Engler T10 / T50-T10 / T90-T10 | Rayleigh T10 / T50-T10 / T90-T10 |
|---|---|---|
| pool C Jet A (lam 0.25) | 175.8 / 1.8 / -5.7 | 175.0 / 4.2 / 7.3 |
| P0 Jet A | 125.7 / 49.6 / 49.6 | 169.4 / 10.1 / 16.8 |
| P0x-a Jet A | 134.6 / 47.2 / 54.3 | 175.4 / 8.7 / 18.5 |
| P0x-b Jet A | 136.6 / 42.2 / 54.9 | 175.0 / 7.6 / 17.1 |
| P1 Jet A | 168.8 / 32.8 / 55.9 | 181.9 / 21.3 / 53.4 |
| P2 Jet A (smoke stand-in) | 173.2 / 33.9 / 102.5 | 181.6 / 20.9 / 93.2 |
| P3 Jet A (smoke stand-in) | 173.5 / 33.9 / 91.5 | 183.0 / 19.9 / 87.5 |

Reading:
* For **pure-compound recipes with a light minor component** (P0, P0x) the two models disagree by about 35-45 C on the slopes. Engler's ~50 C
  T90-T10 comes from linear interpolation of a TBP step curve through 4-5 wt % isooctane; Rayleigh (which also under-predicts real-fuel slopes
  by about 12 C, so even +12 C does not close the gap) puts them at 17-19 C. **Conclusion: pure pools P0/P0x are far from T90-T10 >= 40 C
  under a physical model; the WS7/R-SB "fails by ~4-9 C" understates the gap.** The earlier result (0 certified) does not change.
* For real-fuel blends (P1, P2, P3) both models agree within 0-23 C and every D86 limit is met with both models' own margins.
