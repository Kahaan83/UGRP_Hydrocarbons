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
