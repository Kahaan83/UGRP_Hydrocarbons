# delta_ws3_pure_gaps (WS3: Tm / dHfus and heavy n-alkane Antoine gap filling)

Declared 2026-11-14 before any script of this delta exists; this commit contains no result file. Frozen reference: branch v5 commit
`50f4333`. Plan: `docs/superpowers/plans/2026-11-14-v5-ws2-ws3-freezing.md` (WS3 of the 2026-11-13 roadmap; the YSI/TSI and LFL parts are
out of scope). The declared text below is never edited after results; results are appended under "Result".

## Rules that hold whatever the result
* A filled value never replaces a measured one. Tiers, in order: `verified` (`data/raw/sle_constants/tm_dhfus_verified.csv`, read-only),
  `thermoml_measured` (pure-component rows of the ThermoML 2020 archive, this delta), `ml_estimated` (this delta). Each written value carries
  its tier in `source`.
* Nothing in `src/fuelblend/api/`, `src/fuelblend/design/` or `configs/` reads any file written here. Wiring a filled value into the API,
  the designs or the flash registry is a **user decision**, proposed in the result, never done here. No yaml margin is changed.
* Train / evaluate only on measured values (no rule-labelled targets). Joback (thermo 0.6.1 `Joback`) is a comparator, never a label.

## Part A: Antoine for heavy n-alkanes (C11 low T, C13, C15, C17, C18, C19)
Data: measured pure vapour pressures in `data/unified/measurements_dedup.csv` (`property == vapor_pressure_kPa`, `kind == pure`; all are
ThermoML rows); refit constants `data/processed/antoine_refit.csv` (kPa form, `log10(P/kPa) = A - B/(T/K + C)`).
* **Knot rule (fixed now).** An n-alkane is a knot only if it has >= 6 measured ThermoML points and its `antoine_refit.csv` (A,B,C) reproduces
  its own measured points with rms(log10 P_pred - log10 P_meas) <= 0.03. Expected knots: C8-C10, C12, C14, C16, C20-C22 (C11 has 9 points
  only at 393-459 K and is tested as a knot too). The C13, C15, C17, C18, C19 refits have no measured points in `measurements_dedup.csv`
  (13 points each, exactly 313.15-373.15 K, provenance not established), so they can never be knots.
* **Designated method (fixed now).** For every T on the grid 293.15, 298.15, ..., 453.15 K, log10 P of each knot (its refit) is interpolated
  in carbon number with a monotone PCHIP; the target compound's (A,B,C) is then refit (least squares in log10 P) to its interpolated grid.
  Validity range written = 293.15-453.15 K. Comparator: literal per-parameter PCHIP interpolation of A, B and C separately.
* **R-W3-1 (descriptive):** knot table (n points, T range, rms of the refit on its own measured points, knot yes/no).
* **R-W3-2 (gate):** leave-one-knot-out over interior knots (each knot with a knot on both sides is removed and re-predicted from the rest);
  score = rms log10 error on the held-out compound's measured points inside 293.15-453.15 K. **Pass iff the pooled rms <= 0.03 AND each of
  the 4-carbon-gap cases (C12 held out with C10 and C14 present; C14 held out with C12 and C16 present) has rms <= 0.03.** The C16-C20 gap
  that C17-C19 sit in is 4 carbons wide, which is why those two cases are named.
* **Rule.** If R-W3-2 passes, `data/processed/v5/antoine_interpolated_nalkanes.csv` is written for C13, C15, C17, C18, C19 with
  `antoine_source = interpolated`, `gate = passed`. If it fails, the file is written with `gate = failed` and the note "not for use".
  In both cases the existing registry constants are not changed.
* **Descriptive extras:** C11 interpolated from its neighbours vs its 9 measured points (393-459 K); interpolated vs the existing
  `antoine_refit.csv` and v2 flash-registry (`compounds_anchored_vp.csv`, `anchored_thermo_ambient`) constants for C13-C19 at 313.15 and
  373.15 K. Known before running: the v2 registry has implausible constants from about C19 upward (C20: A = 17.05, C = +150; C27+ A = -5);
  reported, not fixed.

## Part B: Tm / dHfus
* **ThermoML harvest (fixed now).** Pure-component (`n_comp == 1`), all-hydrocarbon datasets of the 2020 archive with property
  `Normal melting temperature, K` or `Triple point temperature, K` (-> Tm) and `Molar enthalpy of transition or fusion, kJ/mol` counted as
  fusion only when the dataset's phases are exactly one Crystal phase plus Liquid. Pressure > 1000 kPa dropped. Per compound: median, min,
  max, n rows, n files; spread > 5 K (Tm) or > 15 % (dHfus) is flagged `spread`, not dropped. Output `data/processed/v5/thermoml_fusion_pure_hc.csv`.
  Caveat stated now: for compounds with a rotator (solid-solid) transition near the melt the ThermoML value is fusion only; the verified table's
  `svl_dhfus_kJ_mol` ("total" for odd n-alkanes C11/C13/C15) takes precedence wherever it exists.
* **Training set.** Hydrocarbons with a measured Tm (verified first, else thermoml_measured median). dHfus target = verified `svl_dhfus_kJ_mol`,
  else the thermoml_measured median. Model target for the enthalpy is dSfus = dHfus / Tm (J/mol/K).
* **Compound family (grouped CV key, fixed now)**, from the SMILES by `family_of()`: `n_alkane` (acyclic, saturated, unbranched),
  `iso_alkane` (acyclic saturated branched), `acyclic_unsaturated`, `cycloalkane_mono` (saturated, 1 ring), `cycloalkane_poly` (saturated,
  >= 2 rings), `cycloalkene` (non-aromatic unsaturated ring), `benzene_mono` (exactly one aromatic ring, no saturated ring),
  `aromatic_naphthenic` (aromatic ring plus a saturated ring, e.g. tetralin, indan), `aromatic_poly` (>= 2 aromatic rings, fused or not).
  CV = leave-one-family-out (every family is one fold). A compound is never in train and test together; no family is.
* **Model (fixed now).** LightGBM regressor (n_estimators 400, learning_rate 0.03, num_leaves 15, min_child_samples 5, seed 0) on RDKit
  descriptors: carbon count, MW, ring count, aromatic ring count, fraction sp3, number of rotatable bonds, branching (CH3 count, quaternary
  C count), Wiener index, Kier kappa1-3, Balaban J, molecular symmetry number (count of symmetry classes / atoms from canonical ranks),
  Joback Tm and Joback dHfus as features. Two models: Tm (K) and dSfus.
* **R-W3-3 (gate): leave-one-family-out OOF MAE(Tm) <= 5 K on hydrocarbons with 8 <= n_carbon <= 20.** Joback MAE on the same rows is
  reported. **Expected outcome stated now: fail.** Strict leave-family-out means the n-alkanes are predicted with no n-alkane in training;
  published group-contribution Tm errors are tens of K.
* **R-W3-4 (descriptive):** OOF MAE of dHfus (= predicted dSfus x predicted Tm) in kJ/mol and in %, same rows; Joback comparator.
* **R-W3-5 (descriptive):** per-family OOF MAE / RMSE and, as a non-decision reading, random 5-fold-by-compound OOF MAE (to show how much of
  the error is family transfer).
* **Rule.** Whatever R-W3-3 says, `data/processed/v5/fusion_ml_estimated.csv` is written for the gap-table compounds that have no verified
  and no thermoml_measured value, with `source = ml_estimated`, `oof_rmse_Tm_K`, `oof_rmse_dHfus_kJ_mol`, `gate_R_W3_3 = passed|failed`. A
  failed gate means "not for use in any margin"; even a pass gives only non-verified constants (freezing stays uncertified for such a member).

## Gap table (R-W3-6, descriptive)
Compounds: every member of the design pools `v5pool`, `ext` (`scripts/v5/run_design_userpool_v5.py extended_pool()`), pool A (+ toluene,
isooctane) and pool B (+ hexylbenzene, isooctane) of `delta_smoke_cert`, and every SubRepsV5 representative selected when the 113 mappable
AJFTD GCxGC fuels are mapped. Per compound: tier of Tm and of dHfus (`verified` / `thermoml_measured` / `needs_estimate`).

## Propagation (R-W3-7, descriptive)
On the final recipes of the pool A and pool B designs (`results/v5/BASE/delta_smoke_cert/designs/pool*/`), recompute ideal SvL (mole
fraction, `ideal_sle_C`) with constants verified > thermoml_measured > ml_estimated, then shift every non-measured component's Tm by +-1 OOF
RMSE(Tm) and dHfus by +-1 OOF RMSE(dHfus) (all four sign combinations) and report the largest change of the freezing value. Recipes are not
re-optimised.

## Outputs
`scripts/v5/base/delta_ws3_antoine.py`, `scripts/v5/base/delta_ws3_fusion.py`; files in this folder: `antoine_knots.csv`, `antoine_loo.csv`,
`antoine_targets.csv`, `fusion_oof.csv`, `fusion_family_summary.csv`, `gap_table.csv`, `propagation.csv`, `readings.json`.

## Result (appended 2026-11-14; declared text above unchanged)
Scripts `scripts/v5/base/delta_ws3_antoine.py`, `scripts/v5/base/delta_ws3_fusion.py`, extractor `src/fuelblend/data/thermoml_fusion.py`;
readings in `readings_antoine.json` and `readings_fusion.json` (two files instead of the one `readings.json` named above). Tests
`tests/test_ws3_antoine_interp.py`, `tests/test_thermoml_fusion.py` (trimmed real record 10.1016/j.fluid.2017.03.012), `tests/test_ws3_fusion_estimate.py`.

### Part A: Antoine
| reading | value | verdict |
|---|---|---|
| R-W3-1 knots | C8, C9, C10, C11, C21, C22. NOT knots: C12 (stored refit rms 0.087 on its own 74 measured points), C14 (0.0315), C16 (0.082), C20 (0.071), and C13/C15/C17/C18/C19 (no measured points) | descriptive |
| R-W3-2 | pooled LOO rms 0.016 (comparator per-parameter 0.017), but the two declared 4-carbon-gap cases cannot be formed (C12, C14 are not knots); the only interior knots are C9, C10 (gap 2) and C11, C21 (gap 11) | **fail** |

* Rule applied: `data/processed/v5/antoine_interpolated_nalkanes.csv` written for C13/C15/C17/C18/C19 with `gate_R_W3_2 = failed`, "not for use".
  No registry changed.
* Finding (not fixed): `data/processed/antoine_refit.csv` is stale against `data/unified/measurements_dedup.csv`. For C14 and C16 a fresh refit
  on the same measured points has rms 0.009 / 0.010, so the data are fine and the stored constants are not; for C12 and C20 the measured references
  disagree (fresh-refit rms 0.056 / 0.051). The stored C13/C15/C17/C18/C19 "refit_measured" sets (13 points each, exactly 313.15-373.15 K) have no
  measured source in the current unified data; their provenance is not established.
* Exploratory, post-hoc, NOT a gate (`antoine_exploratory_fresh_refit.csv`): with knots refit freshly on measured points (C8-C11, C14, C16, C21,
  C22) the leave-one-out errors are C11 0.002 (gap 4, 7 points only), C14 0.049 (gap 5), C16 0.072 (gap 7), C21 0.049 (gap 6). Interpolation across
  5-7 carbons misses the 0.03 target, so even fixed knots would not certify interpolated C17-C19 constants.
* Descriptive: at 313 K the interpolated C13-C19 log10 P agree with the stored refit within 0.03, while the v2 flash registry
  (`anchored_thermo_ambient`) is 0.03-0.12 log10 higher for C15-C19 (C19: +0.12, i.e. about 30 % higher Psat). The v2 registry constants from about
  C19 upward are implausible (C20 A = 17.05, C = +150; C27+ A = -5), as stated in advance; not fixed.

### Part B: Tm / dHfus
* ThermoML harvest: 321 rows, 102 hydrocarbons (`data/processed/v5/thermoml_fusion_pure_hc.csv`, `_by_compound.csv`). Implementation detail
  found on real records before any fit: ThermoML lists the surrounding atmosphere ("Air at 1 atmosphere", "Gas") as a dataset phase, and files
  crystal-crystal transitions under "Triple point temperature" and "Molar enthalpy of transition"; the phase rule therefore ignores the ambient
  phase and requires exactly one Crystal + Liquid for all three properties. Against the 29 verified compounds it also covers, Tm agrees within
  2.1 K (median |dT| 0.2 K); dHfus differs by more than 2.5 kJ/mol for n-C14 (+2.9), n-C16 (-3.3), n-C12 (-2.5), C36 (+3.2), 2,3-DMN (+4.6),
  1,8-DMN (+2.8) (verified values take precedence).
* Training set: 109 hydrocarbons with measured Tm (38 verified, 71 thermoml_measured); 94 with dHfus. Families: aromatic_poly 48, n_alkane 27,
  benzene_mono 10, cycloalkane_poly 7, cycloalkane_mono 5, iso_alkane 5, cycloalkene 3, aromatic_naphthenic 2, acyclic_unsaturated 2.

| reading | value | verdict |
|---|---|---|
| R-W3-3 leave-one-family-out OOF MAE(Tm), C8-C20 (n = 70) | 47.7 K (RMSE 60.4); Joback 40.7 K | **fail** (threshold 5 K; failure expected and stated in advance) |
| R-W3-4 OOF dHfus, C8-C20 (n = 60) | MAE 8.5 kJ/mol (42 %), RMSE 12.5; Joback 6.3 kJ/mol | descriptive (the model is worse than Joback) |
| R-W3-5 random 5-fold by compound, C8-C20 | MAE 47.7 K (n-alkanes alone 10.2 K) | descriptive: even without family transfer the data are too few and too aromatic-heavy |
| R-W3-6 gap table (55 compounds) | Tm: 15 verified, 10 thermoml_measured, 30 need an estimate; dHfus: 15 verified, 6 thermoml_measured, 34 need an estimate | descriptive |
| R-W3-7 propagation (pool A / pool B final recipes) | 0.0 C: no final recipe contains an estimated component (toluene and hexylbenzene end at 0 %); tiered SvL -40.83 C (Jet A) / -47.96 C (Jet A-1) | descriptive |

* Rule applied: `data/processed/v5/fusion_ml_estimated.csv` written for the 35 gap compounds missing a measured Tm and/or dHfus (30 Tm and
  34 dHfus estimates; where Tm is measured only dHfus is estimated, e.g. toluene 6.83 kJ/mol on the measured 178.1 K) (`source = ml_estimated`,
  `oof_rmse_Tm_K = 60.4`, `oof_rmse_dHfus_kJ_mol = 12.5`, `gate_R_W3_3 = failed`, "not for use in any margin"). Nothing reads it. The estimates are
  poor (e.g. n-nonane 194.8 K against about 219.7 K in the literature).
* Pool members: all v5pool / ext members and isooctane have verified Tm and dHfus. Toluene (pool A) has a ThermoML-measured Tm but no dHfus;
  hexylbenzene (pool B) has neither in these tiers. Observation (not changed): the pool A designs record freezing as uncertified because of
  "toluene.Tm=coutinho_table" although the final recipes contain 0 % toluene; the blocker comes from pool membership, not from the recipe.

### Post-result notes (added after review; no tier, gate or value changed)
* Gap-table scope as run differs from the declared wording: the AJFTD representatives come from `get_ctx().reps` (class `SubReps`) mapped over
  all 116 fuels in the GCxGC composition file (not only the 113 mappable ones), plus the 9 `SubRepsV5.CANDIDATES`; not a SubRepsV5 selection run.
* The declared tier rule did not consult `data/processed/v5/pure_props_v5.csv`, which lists Tm / dHfus with status `verified` (source mostly
  `inverse_design_v2.CANDIDATE_REGISTRY`, not primary calorimetry citations). `scripts/v5/base/delta_ws3_annotate_pure_props.py` adds columns
  `tm_pure_props_v5` / `dh_pure_props_v5` to `gap_table.csv` and `fusion_ml_estimated.csv` (22 of 55 gap compounds have an entry). 8 rows of
  `fusion_ml_estimated.csv` hold an ML estimate where pure_props_v5 lists a `verified` value: toluene (dHfus 6.64), hexylbenzene (Tm -61 C),
  tetralin (Tm -35.8 C), o-xylene (Tm -25.2 C), n-octane (Tm -56.8 C, dHfus 20.73), mesitylene (Tm -44.7 C, dHfus 9.51), n-hexane (Tm -95.3 C,
  dHfus 13.08), 2-methylpentane (Tm -153.7 C, dHfus 6.27). Those estimates must not be used where the listed value exists; the file stays unread.
* "Nothing reads it" checked by grep over `src/` and `scripts/` for `thermoml_fusion`, `fusion_ml_estimated`, `antoine_interpolated` and for globs
  over `data/processed`: the only readers are the delta scripts of this folder and `scripts/v5/base/delta_ws2_freeze_multicomp.py` (ThermoML tier only).

### Needs user decision / follow-ups
* Whether to refresh `antoine_refit.csv` from `measurements_dedup.csv` (C14, C16 fit their measured data to 0.01 log10 when refit). Not done here.
* Whether to curate primary-source Tm/dHfus (WebBook calorimetry, the `tm_dhfus_verified.csv` method) for toluene, hexylbenzene, n-octane, n-nonane,
  n-hexane, p-xylene and the gap-table representatives instead of estimating them: the gated estimator cannot fill them (47.7 K error).
* Acree-Chickos compendia were not ingested (large PDF tables); they remain the main open source of dHfus for alkylbenzenes, alkylcyclohexanes and
  branched alkanes.
