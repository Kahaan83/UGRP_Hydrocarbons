# v6 Track C: mixture ML retry (declared before any mixture fit)

Declared 2026-11-21, before any candidate or comparator was fitted or scored. Lock: `lock.json`, rows: `rows.csv`
(`scripts/v6/mixture/build_tables.py`). Nothing here is exposed in `predict_fuel` or certified; design certification stays model-level only.

## Candidate
Frozen Track B D-MPNN encoder (`scripts/v6/mixture/pretrain_encoder.py`, fitted once on the full pure pool, config hidden 300 / depth 3 /
dropout 0.2) gives one embedding per component. Pretraining on all pure labels is fair for the mixture task: the baselines also receive
measured pure values, the pool holds pure rows only, and no mixture label enters it. The encoder is not fine-tuned (tens of systems
cannot support it). On top: a per-component projection phi, mole-fraction pooling z = sum_i x_i phi_i (MolPool), a bilinear pair term
sum_{i<j} x_i x_j phi_i^T W phi_j, and a small MLP (plus 1000/T for the Boehm heads). Output is a **residual over the physics baseline**.
Hyperparameters (projection width, weight decay, epochs) are chosen per outer fold by inner GroupKFold(3) on the training groups;
final prediction is a 5-seed ensemble.

## Heads and data (after dedup; counts in `build_log.json`)
| head | source | baseline | groups (dev / locked) | rows (dev / locked) |
|---|---|---|---|---|
| kaust_ron, kaust_mon | KAUST Table 4 blends, vol% converted to mole with `kaust_registry.densities()` | linear molar | 5 / 3 families | 90 / 30 each |
| leenhouts_mon | Leenhouts/Kuzhagaliyeva MON mixtures | linear molar | 34 / 18 systems | 181 / 60 |
| boehm_ln_visc, boehm_density | Boehm 2022 cold-T (production loader output), mixtures only | Refutas (mass fractions, pure at same T); ideal volume | 12 / 11 systems | 303 / 95 each |

Octane heads ignore the `is_hc` filter on purpose (ethanol is in scope). Groups: KAUST family = component set without ethanol (nests
`system_id` and keeps ethanol blends with their base fuel); Leenhouts and Boehm `system_id`. Lock rule (outcome-blind): groups in
md5(`v6c:`+group) order are locked while locked rows stay <= 25 % of the head's rows. The first rule tried (first 20 % of groups) put
the largest KAUST family into the lock, leaving 49 dev rows, so it was replaced before any fit.

Exclusions, logged: Co-Optima (the only repo file, `data/processed/v2/octane_co_optima.csv`, is pure compounds and the blend file
`fuel_engine_co_optimization.csv` is absent; BOB-01 has no molecular composition); KAUST Table 5 FACE (no composition); 51 Leenhouts
rows whose pure-MON coverage (mole fraction with a measured pure MON in the pool) is < 0.9, dropped for baseline and candidate alike;
8 Boehm rows with no same-T pure value. Leenhouts linear molar is renormalised over the covered fraction (>= 0.9).

## Comparators (all fitted inside the same outer folds)
- Octane: (a) linear molar; (b) "corrected-pure" (our definition: linear molar with per-component effective blending values fitted
  in-fold by ridge shrinkage toward the measured pure value, lambda by inner CV); (c) linear molar + in-fold GBM on mole-weighted
  descriptors. Linear volume is reported for KAUST only.
- Boehm: (a) physics; (b) physics + in-fold per-temperature mean residual; (c) physics + in-fold GBM on mole-weighted descriptors and T.
  The v5 LightGBM pair models are not used: Boehm rows are part of the v5 data tables (`delta_ml_viscosity.py` reads them), so
  using them would overlap train and evaluation.
- Gate comparator = the comparator with the lowest **dev OOF** MAE (chosen on dev before the locked set is scored).

## Gate (per head, same logic as v5 WS4)
Pass only if all three hold: (1) dev out-of-fold paired group bootstrap, 90 % CI of (MAE_best_comparator - MAE_candidate) has lower
bound > 0; (2) the same on the locked groups (scored once); (3) MAE on rows with n_comp >= 3 is no worse than the comparator on dev
and locked (reported as not evaluable where a stratum is empty). Bootstrap resamples groups, 2000 draws. 5 heads, no multiplicity
correction. Conformal margins are fitted only for a head that passes.

## Expected outcome (stated now)
A tie or a loss. Dev has 5 KAUST families over 8 components, 34 Leenhouts systems and 12 Boehm systems; the GBM residual already
beats linear molar on KAUST, and v5 WS4 found no neural mixture encoder passing against LightGBM pair residuals. A negative result is
a valid result and goes to `docs/v6_negative_results.md`.

## Results (scored 2026-11-21: dev OOF first, then locked once; `summary.csv`, `readings.json`)
MAE of the gate comparator (best on dev OOF) vs MolPool; gain = comparator minus MolPool, 90 % paired group-bootstrap CI.

| head | gate comparator | dev MAE comp / MolPool | dev gain CI | locked MAE comp / MolPool | locked gain CI | pass |
|---|---|---|---|---|---|---|
| kaust_ron | desc GBM | 1.49 / 2.60 | [-2.95, -0.55] | 3.80 / 3.85 | [-0.23, 0.08] | no |
| kaust_mon | linear molar | 2.91 / 3.35 | [-2.78, 0.01] | 3.67 / 3.55 | [-0.46, 0.46] | no |
| leenhouts_mon | desc GBM | 3.08 / 2.29 | [0.04, 1.33] | 2.42 / 1.97 | [0.07, 0.98] | **yes** |
| boehm_ln_visc | desc GBM | 0.055 / 0.078 | [-0.043, -0.004] | 0.026 / 0.036 | [-0.017, -0.001] | no |
| boehm_density | physics | 2.87 / 3.87 | [-2.01, -0.09] | 2.42 / 2.42 (head chose no correction) | [0, 0] | no |

- Leenhouts MON (n-comp >= 3 rows): dev 2.74 -> 1.89, locked 1.97 -> 1.50, so the third gate holds. Linear molar there is 5.78 dev MAE and
  corrected-pure 6.50, so MolPool beats both by a wide margin; the gate comparator is the descriptor GBM. Out-of-fold 90 % margins
  (`leenhouts_mon_margin.json`): absolute 4.92 MON, one-sided over-prediction 3.43; locked coverage 93 % and 97 %.
- Caveats: one head of five passes with no multiplicity correction (expected false pass rate across 5 heads is not negligible); 34 dev
  and 18 locked systems; rows with pure-MON coverage < 0.9 were dropped; Leenhouts dev rows within a system are correlated (bootstrap
  resamples systems). No lab measurement exists, so this is a model-level result only.
- KAUST: no skill over the in-fold GBM or linear molar (8 components, 5 dev families). Boehm: MolPool is significantly worse on dev;
  for density it fell back to no correction in most folds.
- Corrected-pure (our definition) was worse than plain linear molar on every octane head.
