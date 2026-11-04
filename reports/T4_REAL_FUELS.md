# T4 — Real-composition test set scoring (CRC 647 surrogates, physics-only paths)

> **This is the honest scorecard of the composition→property pipeline on real fuels.**
> Every real-fuel data row sits behind the R1 wall (L0, in `measurements_dedup.csv`,
> never in `model_table.csv`), so nothing here was trainable — these are genuine
> out-of-distribution generalization numbers.

## Setup

- **Data**: `data/unified/measurements_dedup.csv`, `kind=real_fuel`, `sources=crc647`
  (CRC Handbook Report 647). 57 aviation fuels with **measured class composition**
  (saturates / olefins / aromatics vol% + naphthalenes) and measured density (398 rows),
  kinematic viscosity at −20/−40 °C (171 rows), closed-cup flash (57 fuels).
  AJFTD fuels were **excluded** from composition-respecting scoring: only ~82 of 2,672
  carry an aromatics number and it was never ingested into unified.
- **Surrogates** (4 representative molecules, one per class): n-undecane (n-paraffin),
  2,2,4-trimethylpentane (iso + olefins), methylcyclohexane (cycloalkane),
  ethylbenzene (all aromatics incl. naphthalenes). Saturates split by default
  0.60/0.25/0.15 n/iso/cyclo; vol% → mass% via equal-density approximation;
  mole fractions via component MWs.
- **Prediction paths (physics-only)**: density = ideal-volume mixing of pure-GBM
  component densities at T; viscosity = Refutas VBN blending of pure-GBM ν (GBM η/ρ);
  flash = Le Chatelier + Raoult with refit-measured Antoine (kPa convention) + LFL.
  No mixture-residual model is applied here (trained on molecular mixtures, not surrogates).
- Seeds [11, 23, 37, 41, 59]; sanity self-checks pass before scoring:
  GBM n-decane η(293 K)=1.196 mPa·s, ρ=740.5 kg/m³; chat cyclohexane −17.1 °C (lit −18).
- Known data bug found while scoring: CRC density rows carry `T_K=15.000` — the loader
  stored **°C in the Kelvin column** (all 398 rows). Fixed at scoring time only
  (+273.15); `TODO` logged for the CRC loader + rebuild in the next QA pass.

## Numbers (57 fuels)

| Property | n | MAE | RMSE | MAPE | Within 10 % | Note |
|---|---|---|---|---|---|---|
| density_kg_m3 (15 °C) | 398 | 56.3 kg/m³ | 62.0 | 6.9 % | 83.2 % | class-average surrogate can't see the fuel's C-number spread |
| viscosity_kin_cSt (−20 °C & −40 °C) | 171 | 4.12 cSt | 5.21 | 76 % | 0 % | cold-T pure-GBM extrapolation + C-number flattening |
| flash_point_C (Le Chatelier surrogate) | 57 | 53.6 °C | 54.1 | — | 0 % | **bias −53.6 °C: systematic hard under-prediction** |

### Flash: negative result, and why it matters

- Le Chatelier on the class surrogate scores ~−10 °C against measured 30–65 °C.
  The surrogate puts materially volatile molecules (C7–C10, incl. a C7 ring and a C8 aroma)
  at *class-scale* mass fractions, whereas a real jet fuel's flash behavior is set by
  **trace light ends** (<10 % by mole) its class breakdown doesn't resolve. Le Chatelier
  integrates vapour pressure over the *whole* composition — concentrated C7/C8 surrogates
  overwhelm it. The pure-sanity path works (cyclohexane −17.1 vs lit −18 °C), so the math
  is not the failure mode; the granularity is.
- Direction of the bias is the *safe* side (under-prediction; 0 % over-prediction >5 °C),
  but a −50 °C error makes the surrogate useless for certification — not just inaccurate,
  *unusable*.
- GBM mass-weighted mix path (fallback): MAE 24.3 °C — better, still a hard under-predictor
  (mean positive error 0). Same root cause: barycenter-of-classes ≠ volatility envelope.

## Verdict — when to trust what

| Path | Real-fuel flash | Molecular model datasets | Guidance |
|---|---|---|---|
| Le Chatelier + surrogate classes | ✗ (−54 °C bias) | ✓ (S2 binary MAE 1.64 °C) | only where true molecular composition exists |
| Le Chatelier + measured molecule list | — | ✓ | real surrogates for L3 systems |
| Pure-GBM → Refutas / id-volume | ✗ at cold T, ~✓ at RT | ✓ (binary ln-RMSE 0.27) | density usable for screening (MAPE 6.9 %) |
| D86-curve regression (WP7, trained on real fuels) | ✓ (CV MAE 3.36 °C) | n/a | the only real-fuel flash performer — requires distillation data |

## Honest limitations

1. Saturates 0.6/0.25/0.15 split and the equal-density vol%→mass% conversion are
   approximations; both are documented above and can be refined with GC×GC data (not open).
2. naphthalenes folded into ethylbenzene (no refit Antoine for methylnaphthalene) —
   worsens the volatility-envelope distortion further.
3. AJFTD real fuels (worst-case data thinness) remain unscoreable for flash by any
   composition path; this is a data-collection gap (T7), not a modeling gap.
4. CRC density `T_K` unit bug fixed at scoring time only — loader repair still pending.

All numbers computed from measured unified rows with the fixed seeds; no synthetic data
involved in evaluation.
